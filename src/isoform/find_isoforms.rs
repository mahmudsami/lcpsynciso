//! The `find-isoforms` subcommand: detect isoforms for an existing cluster assignment.
//!
//! The reads file does not need to be sorted or grouped by cluster. Clusters are packed into
//! passes of at most `--max-buffer-reads` reads; each pass rereads the file, keeps only its
//! own clusters' reads, and detects their isoforms.

use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::fs;

use crate::seqio::SeqReader;

use super::options::{apply_isoform_flag, ISOFORM_HELP};
use super::output::{format_cluster_records, OutputFiles};
use super::{detect_isoforms, IsoformOptions, Read};

pub struct Config {
    pub clusters_path: String,
    pub reads_path: String,
    pub out_dir: String,
    pub only_clusters: Vec<u32>,
    pub min_cluster_size: usize,
    pub max_cluster_size: usize,
    /// Reads held in memory per pass.
    pub max_buffer_reads: usize,
    pub isoform_options: IsoformOptions,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            clusters_path: String::new(),
            reads_path: String::new(),
            out_dir: "lcpsynciso_isoform_out".into(),
            only_clusters: Vec::new(),
            min_cluster_size: 1,
            max_cluster_size: usize::MAX,
            max_buffer_reads: 5_000_000,
            isoform_options: IsoformOptions::default(),
        }
    }
}

pub fn parse_args(argv: &[String]) -> Config {
    let mut config = Config::default();
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
            "--clusters" => config.clusters_path = next(),
            "--reads" => config.reads_path = next(),
            "-o" | "--out" => config.out_dir = next(),
            "--only-clusters" => {
                config.only_clusters =
                    next().split(',').map(|x| x.trim().parse().unwrap()).collect()
            }
            "--min-size" => config.min_cluster_size = next().parse().unwrap(),
            "--max-size" => config.max_cluster_size = next().parse().unwrap(),
            "--max-buffer-reads" => config.max_buffer_reads = next().parse().unwrap(),
            "--k" => config.isoform_options.minimizer_k = next().parse().unwrap(),
            "--w" => config.isoform_options.minimizer_w = next().parse().unwrap(),
            "-h" | "--help" => {
                print_usage();
                std::process::exit(0);
            }
            other => {
                if !apply_isoform_flag(&mut config.isoform_options, other, &mut next) {
                    eprintln!("unknown arg {other}");
                    print_usage();
                    std::process::exit(2);
                }
            }
        }
        i += 1;
    }
    if config.clusters_path.is_empty() || config.reads_path.is_empty() {
        print_usage();
        std::process::exit(2);
    }
    config
}

pub fn print_usage() {
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

{ISOFORM_HELP}

OUTPUT (in DIR):
    isoform_assignments.tsv   read_name <tab> cluster_id <tab> isoform_id
    isoform_summary.tsv       cluster_id <tab> isoform_id <tab> n_reads <tab> length
    isoforms.fasta            refined consensus per isoform"
    );
}

pub fn run(config: Config) {
    // ── Cluster assignment ──
    eprintln!("[find-isoforms] loading clusters from {} ...", config.clusters_path);
    let mut cluster_of_read: HashMap<String, u32> = HashMap::new();
    let mut cluster_size: HashMap<u32, usize> = HashMap::new();
    for line in fs::read_to_string(&config.clusters_path).expect("read clusters").lines() {
        let mut it = line.split('\t');
        let name = it.next().unwrap_or("");
        // Lines without a numeric cluster id (such as the header) are skipped.
        match it.next().and_then(|c| c.parse::<u32>().ok()) {
            Some(cid) => {
                cluster_of_read.insert(name.to_string(), cid);
                *cluster_size.entry(cid).or_insert(0) += 1;
            }
            None => continue,
        }
    }
    eprintln!(
        "[find-isoforms] {} reads across {} clusters",
        cluster_of_read.len(),
        cluster_size.len()
    );

    let only: HashSet<u32> = config.only_clusters.iter().copied().collect();
    let is_eligible = |cid: u32| {
        (only.is_empty() || only.contains(&cid))
            && cluster_size
                .get(&cid)
                .map_or(false, |&n| n >= config.min_cluster_size && n <= config.max_cluster_size)
    };

    // ── Plan passes ──
    // First-fit decreasing on cluster size: a cluster never spans two passes, and each pass
    // holds at most `budget` reads unless a single cluster is bigger than that.
    let budget = config.max_buffer_reads.max(1);
    let mut eligible_clusters: Vec<u32> =
        cluster_size.keys().copied().filter(|&c| is_eligible(c)).collect();
    eligible_clusters.sort_by_key(|c| Reverse(cluster_size[c]));
    let mut pass_of: HashMap<u32, usize> = HashMap::new();
    let mut pass_load: Vec<usize> = Vec::new();
    for &c in &eligible_clusters {
        let sz = cluster_size[&c];
        let pass = pass_load.iter().position(|&load| load + sz <= budget).unwrap_or_else(|| {
            pass_load.push(0);
            pass_load.len() - 1
        });
        pass_load[pass] += sz;
        pass_of.insert(c, pass);
    }
    let n_passes = pass_load.len().max(1);
    eprintln!(
        "[find-isoforms] {} eligible clusters -> {} streaming pass(es) over {} (<= {} reads/pass) ...",
        eligible_clusters.len(),
        n_passes,
        config.reads_path,
        budget
    );

    // ── Detect isoforms, one pass at a time ──
    let mut output_files = OutputFiles::create(&config.out_dir);
    let (mut total_isoforms, mut n_processed_clusters) = (0usize, 0usize);
    for pass in 0..n_passes {
        let mut reads_by_cluster: HashMap<u32, Vec<Read>> = HashMap::new();
        let mut n_buffered_reads = 0u64;
        for (name, seq) in SeqReader::open(&config.reads_path) {
            match cluster_of_read.get(&name) {
                Some(&cid) if pass_of.get(&cid) == Some(&pass) => {
                    reads_by_cluster.entry(cid).or_default().push(Read { name, seq });
                    n_buffered_reads += 1;
                }
                _ => {}
            }
        }
        eprintln!(
            "[find-isoforms]   pass {}/{}: buffered {} reads in {} clusters; detecting ...",
            pass + 1,
            n_passes,
            n_buffered_reads,
            reads_by_cluster.len()
        );
        // Ascending cluster id, so output order does not depend on HashMap order.
        let mut cids: Vec<u32> = reads_by_cluster.keys().copied().collect();
        cids.sort_unstable();
        for cid in cids {
            let reads = reads_by_cluster.remove(&cid).unwrap();
            let isoforms = detect_isoforms(&reads, &config.isoform_options);
            total_isoforms += isoforms.len();
            n_processed_clusters += 1;
            output_files.append(&format_cluster_records(cid, &reads, &isoforms));
        }
    }
    output_files.flush();
    eprintln!(
        "[find-isoforms] wrote {} isoforms over {} clusters to {}/",
        total_isoforms, n_processed_clusters, config.out_dir
    );
}
