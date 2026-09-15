//! The `find-isoforms` subcommand: resolve isoforms for an existing cluster assignment.
//!
//! The reads file does not need to be sorted or grouped by cluster. Clusters are packed into
//! passes of at most `--max-buffer-reads` reads; each pass rereads the file, keeps only its
//! own clusters' reads, and resolves them.

use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::fs;

use crate::io::SeqReader;

use super::options::{parse_flag, RESOLVE_HELP};
use super::output::{render, OutputFiles};
use super::{resolve_cluster, Cfg, Read};

pub struct Config {
    pub clusters: String,
    pub reads: String,
    pub out_dir: String,
    pub only_clusters: Vec<u32>,
    pub min_size: usize,
    pub max_size: usize,
    /// Reads held in memory per pass.
    pub max_buffer_reads: usize,
    pub cfg: Cfg,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            clusters: String::new(),
            reads: String::new(),
            out_dir: "lcpsynciso_isoform_out".into(),
            only_clusters: Vec::new(),
            min_size: 1,
            max_size: usize::MAX,
            max_buffer_reads: 5_000_000,
            cfg: Cfg::default(),
        }
    }
}

pub fn parse_args(argv: &[String]) -> Config {
    let mut a = Config::default();
    let mut i = 0;
    while i < argv.len() {
        let arg = argv[i].clone();
        let mut next = || {
            i += 1;
            argv.get(i).cloned().unwrap_or_else(|| {
                eprintln!("missing value for {arg}");
                std::process::exit(2);
            })
        };
        match arg.as_str() {
            "--clusters" => a.clusters = next(),
            "--reads" => a.reads = next(),
            "-o" | "--out" => a.out_dir = next(),
            "--only-clusters" => {
                a.only_clusters = next().split(',').map(|x| x.trim().parse().unwrap()).collect()
            }
            "--min-size" => a.min_size = next().parse().unwrap(),
            "--max-size" => a.max_size = next().parse().unwrap(),
            "--max-buffer-reads" => a.max_buffer_reads = next().parse().unwrap(),
            "--k" => a.cfg.k = next().parse().unwrap(),
            "--w" => a.cfg.w = next().parse().unwrap(),
            "-h" | "--help" => {
                usage();
                std::process::exit(0);
            }
            other => {
                if !parse_flag(&mut a.cfg, other, &mut next) {
                    eprintln!("unknown arg {other}");
                    usage();
                    std::process::exit(2);
                }
            }
        }
        i += 1;
    }
    if a.clusters.is_empty() || a.reads.is_empty() {
        usage();
        std::process::exit(2);
    }
    a
}

pub fn usage() {
    eprintln!(
        "lcpsynciso find-isoforms — group reads into isoforms by internal structure (in-process, no external aligner)

USAGE:
    lcpsynciso find-isoforms --clusters assignments.tsv --reads reads.fa[.gz|.bgz] -o DIR [options]

The reads need NOT be sorted or grouped — grouping is done internally in a few
bounded-memory streaming passes over the file (no external tool required).

FILES:
    --clusters FILE         TSV: read_name <tab> cluster_id (header allowed)
    --reads FILE            reads FASTA/FASTQ (plain/.gz/.bgz), in any order
    -o, --out DIR           output directory                             [lcpsynciso_isoform_out]

CLUSTER SELECTION:
    --only-clusters a,b,c   restrict to these cluster ids
    --min-size N            skip clusters with fewer reads               [1]
    --max-size N            skip clusters with more reads                [all]

MEMORY / SPEED:
    --max-buffer-reads N    reads held in RAM per pass; smaller = less   [5000000]
                            memory, more passes over the file

ISOFORMS:
    --k N / --w N           isoform minimizer k-mer length / window      [15 / 10]

{RESOLVE_HELP}

OUTPUT (in DIR):
    isoform_assignments.tsv   read_name <tab> cluster_id <tab> isoform_id
    isoform_summary.tsv       cluster_id <tab> isoform_id <tab> n_reads <tab> length
    isoforms.fasta            refined consensus per isoform"
    );
}

pub fn run(a: Config) {
    // ── Cluster assignment ──
    eprintln!("[find-isoforms] loading clusters from {} ...", a.clusters);
    let mut cl_of: HashMap<String, u32> = HashMap::new();
    let mut size: HashMap<u32, usize> = HashMap::new();
    for line in fs::read_to_string(&a.clusters).expect("read clusters").lines() {
        let mut it = line.split('\t');
        let name = it.next().unwrap_or("");
        // Lines without a numeric cluster id (such as the header) are skipped.
        match it.next().and_then(|c| c.parse::<u32>().ok()) {
            Some(cid) => {
                cl_of.insert(name.to_string(), cid);
                *size.entry(cid).or_insert(0) += 1;
            }
            None => continue,
        }
    }
    eprintln!("[find-isoforms] {} reads across {} clusters", cl_of.len(), size.len());

    let only: HashSet<u32> = a.only_clusters.iter().copied().collect();
    let is_target = |cid: u32| {
        (only.is_empty() || only.contains(&cid))
            && size.get(&cid).map_or(false, |&n| n >= a.min_size && n <= a.max_size)
    };

    // ── Plan passes ──
    // First-fit decreasing on cluster size: a cluster never spans two passes, and each pass
    // holds at most `budget` reads unless a single cluster is bigger than that.
    let budget = a.max_buffer_reads.max(1);
    let mut targets: Vec<u32> = size.keys().copied().filter(|&c| is_target(c)).collect();
    targets.sort_by_key(|c| Reverse(size[c]));
    let mut pass_of: HashMap<u32, usize> = HashMap::new();
    let mut pass_load: Vec<usize> = Vec::new();
    for &c in &targets {
        let sz = size[&c];
        let pi = pass_load.iter().position(|&load| load + sz <= budget).unwrap_or_else(|| {
            pass_load.push(0);
            pass_load.len() - 1
        });
        pass_load[pi] += sz;
        pass_of.insert(c, pi);
    }
    let n_passes = pass_load.len().max(1);
    eprintln!(
        "[find-isoforms] {} target clusters -> {} streaming pass(es) over {} (<= {} reads/pass) ...",
        targets.len(),
        n_passes,
        a.reads,
        budget
    );

    // ── Resolve, one pass at a time ──
    let mut out = OutputFiles::create(&a.out_dir);
    let (mut total_iso, mut n_clusters) = (0usize, 0usize);
    for p in 0..n_passes {
        let mut groups: HashMap<u32, Vec<Read>> = HashMap::new();
        let mut nread = 0u64;
        for (name, seq) in SeqReader::open(&a.reads) {
            match cl_of.get(&name) {
                Some(&cid) if pass_of.get(&cid) == Some(&p) => {
                    groups.entry(cid).or_default().push(Read { name, seq });
                    nread += 1;
                }
                _ => {}
            }
        }
        eprintln!(
            "[find-isoforms]   pass {}/{}: buffered {} reads in {} clusters; resolving ...",
            p + 1,
            n_passes,
            nread,
            groups.len()
        );
        // Ascending cluster id, so output order does not depend on HashMap order.
        let mut cids: Vec<u32> = groups.keys().copied().collect();
        cids.sort_unstable();
        for cid in cids {
            let reads = groups.remove(&cid).unwrap();
            let isos = resolve_cluster(&reads, &a.cfg);
            total_iso += isos.len();
            n_clusters += 1;
            out.append(&render(cid, &reads, &isos));
        }
    }
    out.flush();
    eprintln!(
        "[find-isoforms] wrote {} isoforms over {} clusters to {}/",
        total_iso, n_clusters, a.out_dir
    );
}
