//! The `predict` subcommand: reads -> clusters -> isoforms in one process.
//!
//! The input file is read exactly once. Every read is kept 2-bit packed in memory
//! ([`ReadStore`]), so clustering hands sequences straight to isoform detection. Three phases:
//!
//!   A. Read and pack every read, counting how many reads contain each seed.
//!   B. Cluster the stored reads greedily, longest first ([`GreedyClusterer`]).
//!   C. Detect the isoforms of each cluster ([`isoform::detect_isoforms`]), in parallel across
//!      clusters, largest first, and write the output files.
//!
//! Phases A and B are shared with the `cluster` subcommand and live in the `cluster` module
//! (`load_reads`, `cluster_reads`).

use std::cmp::Reverse;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering::Relaxed};
use std::sync::Mutex;
use std::time::Instant;

use rayon::prelude::*;

use crate::cluster::{
    check_max_occ_frac, cluster_reads, effective_max_occ, informative_seeds, load_reads,
    GreedyClusterer, SeedParams,
};
use crate::isoform::{
    self, apply_isoform_flag, ClusterRecords, IsoformOptions, OutputFiles, PlacedRead, Read,
    StartSplit, StartStats, ISOFORM_HELP,
};
use crate::read_store::ReadStore;

pub struct Config {
    pub reads_path: String,
    pub out_dir: String,
    pub threads: usize,
    pub batch_size: usize,
    pub max_reads: usize,
    // Clustering.
    pub seed_k: usize,
    pub seed_s: usize,
    pub seed_t: usize,
    pub levels: usize,
    pub min_occ: u32,
    pub max_occ: u32,
    pub max_occ_frac: f64,
    pub min_shared: u32,
    pub min_shared_frac: f64,
    pub level_weights: Vec<u32>,
    // Which clusters to process, by read count.
    pub min_cluster_size: usize,
    pub max_cluster_size: usize,
    // Write each read's phase B cluster to `gene_clusters.tsv`.
    pub dump_clusters: bool,
    // Write what start validation sees to this file (`--dump-starts`; diagnostics only).
    pub dump_starts: Option<String>,
    // Isoform detection.
    pub isoform_options: IsoformOptions,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            reads_path: String::new(),
            out_dir: "lcpsynciso_out".into(),
            threads: 0,
            batch_size: 500_000,
            max_reads: usize::MAX,
            seed_k: 15,
            seed_s: 9,
            seed_t: usize::MAX,
            levels: 3,
            min_occ: 2,
            max_occ: 100_000,
            max_occ_frac: 0.01,
            min_shared: 3,
            min_shared_frac: 0.10,
            level_weights: vec![1, 2, 5],
            min_cluster_size: 1,
            max_cluster_size: usize::MAX,
            dump_clusters: false,
            dump_starts: None,
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
            "-o" | "--out" => config.out_dir = next(),
            "--threads" => config.threads = next().parse().unwrap(),
            "--batch" => config.batch_size = next().parse().unwrap(),
            "--max-reads" => config.max_reads = next().parse().unwrap(),
            "--k" => config.seed_k = next().parse().unwrap(),
            "--s" => config.seed_s = next().parse().unwrap(),
            "--t" => config.seed_t = next().parse().unwrap(),
            "--levels" => config.levels = next().parse().unwrap(),
            "--min-occ" => config.min_occ = next().parse().unwrap(),
            "--max-occ" => config.max_occ = next().parse().unwrap(),
            "--max-occ-frac" => config.max_occ_frac = next().parse().unwrap(),
            "--min-shared" => config.min_shared = next().parse().unwrap(),
            "--min-shared-frac" => config.min_shared_frac = next().parse().unwrap(),
            "--level-weights" => {
                config.level_weights =
                    next().split(',').map(|x| x.trim().parse().unwrap()).collect()
            }
            "--min-size" => config.min_cluster_size = next().parse().unwrap(),
            "--max-size" => config.max_cluster_size = next().parse().unwrap(),
            "--iso-k" => config.isoform_options.minimizer_k = next().parse().unwrap(),
            "--iso-w" => config.isoform_options.minimizer_w = next().parse().unwrap(),
            "--dump-clusters" => config.dump_clusters = true,
            "--dump-starts" => config.dump_starts = Some(next()),
            "-h" | "--help" => {
                print_usage();
                std::process::exit(0);
            }
            other => {
                if !apply_isoform_flag(&mut config.isoform_options, other, &mut next) {
                    if other.starts_with('-') {
                        eprintln!("unknown flag {other}");
                        std::process::exit(2);
                    }
                    config.reads_path = other.into();
                }
            }
        }
        i += 1;
    }
    if config.reads_path.is_empty() {
        print_usage();
        std::process::exit(2);
    }
    if config.seed_t == usize::MAX {
        config.seed_t = (config.seed_k - config.seed_s) / 2;
    }
    check_max_occ_frac(config.max_occ_frac);
    config
}

pub fn print_usage() {
    eprintln!(
        "lcpsynciso predict — fused reads-to-isoforms pipeline (file read once, all in RAM)

USAGE:
    lcpsynciso predict <reads.fq[.gz]|.fa[.gz|.bgz]> [options]

GENERAL:
    -o, --out DIR           output directory                             [lcpsynciso_out]
    --threads N             worker threads (0 = all)                     [0]
    --batch N               reads per batch                              [500000]
    --max-reads N           stop after N reads                           [all]
    --dump-clusters         also write gene_clusters.tsv (phase B clusters)
    --dump-starts FILE      write what start validation sees (start peaks, labelled
                            and capped reads per structure group) to FILE

CLUSTERING (seeds):
    --k N / --s N           k-mer / s-mer length                         [15 / 9]
    --t N                   open-syncmer offset of the minimal s-mer     [(k-s)/2]
    --levels N              use seeds from levels L1..N                  [3]
    --min-occ N             keep seed if it occurs in >= N reads         [2]
    --max-occ N             ... and in <= N reads                        [100000]
    --max-occ-frac F        ... or in <= F x reads, if that is more      [0.01]
    --level-weights a,b,c   per-level vote weights                       [1,2,5]
    --min-shared N          join a cluster if weighted shared >= N       [3]
    --min-shared-frac F     ... and >= F of the read's seed weight       [0.1]

ISOFORMS:
    --min-size N            skip clusters with fewer reads               [1]
    --max-size N            skip clusters with more reads                [all]
    --iso-k N / --iso-w N   isoform minimizer k-mer length / window      [15 / 10]

{ISOFORM_HELP}

OUTPUT (in DIR):
    isoform_assignments.tsv   read_name <tab> cluster_id <tab> isoform_id
    isoform_summary.tsv       cluster_id <tab> isoform_id <tab> n_reads <tab> length
    isoforms.fasta            refined consensus per isoform
    gene_clusters.tsv         read_name <tab> cluster_id <tab> cluster_size (--dump-clusters)"
    );
}

pub fn run(config: Config) {
    if config.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(config.threads).build_global().ok();
    }
    let seed_params = SeedParams {
        k: config.seed_k,
        s: config.seed_s,
        t: config.seed_t,
        levels: config.levels,
    };
    let run_start = Instant::now();

    // ── Phase A: read, pack, count seeds ──
    eprintln!(
        "[predict] phase A: reading {} once (pack + seed frequencies) ...",
        config.reads_path
    );
    let (store, counts, [d_read, d_seed, d_pack]) =
        load_reads(&config.reads_path, config.batch_size, config.max_reads, seed_params, "predict");
    let n_reads = store.len();
    let n_distinct_seeds = counts.len();
    let max_occ = effective_max_occ(config.max_occ, config.max_occ_frac, n_reads);
    let informative = informative_seeds(&counts, config.min_occ, max_occ);
    drop(counts);
    let phase_a_time = run_start.elapsed();
    eprintln!(
        "[predict] {} reads, distinct seeds={}, informative (in [{},{}])={}",
        n_reads,
        n_distinct_seeds,
        config.min_occ,
        max_occ,
        informative.len()
    );
    eprintln!(
        "[predict] === phase A took {:.1}s  [read/decompress={:.1}s  seed-freq={:.1}s  2bit-pack={:.1}s] ===",
        phase_a_time.as_secs_f64(),
        d_read.as_secs_f64(),
        d_seed.as_secs_f64(),
        d_pack.as_secs_f64()
    );

    // ── Phase B: cluster ──
    let phase_b_start = Instant::now();
    eprintln!(
        "[predict] phase B: greedy clustering (min-shared={}, min-shared-frac={}) ...",
        config.min_shared, config.min_shared_frac
    );
    let clusterer = GreedyClusterer::new(config.min_shared, config.min_shared_frac);
    let cluster_members = cluster_reads(
        &store,
        seed_params,
        informative,
        &config.level_weights,
        clusterer,
        config.batch_size,
        "predict",
    );
    let n_clusters = cluster_members.len();
    let phase_b_time = phase_b_start.elapsed();
    eprintln!("[predict] === phase B (clustering) took {:.1}s ===", phase_b_time.as_secs_f64());
    if config.dump_clusters {
        write_clusters(&config.out_dir, &cluster_members, &store);
    }

    // ── Phase C: detect isoforms ──
    if let Some(path) = &config.dump_starts {
        isoform::stats::open_start_dump(path);
    }
    let phase_c_start = Instant::now();
    let mut eligible_clusters: Vec<u32> = (0..n_clusters as u32)
        .filter(|&c| {
            let size = cluster_members[c as usize].len();
            size >= config.min_cluster_size && size <= config.max_cluster_size
        })
        .collect();
    eligible_clusters.sort_unstable_by_key(|&c| Reverse(cluster_members[c as usize].len()));
    eprintln!(
        "[predict] phase C: detecting isoforms for {} clusters (parallel) ...",
        eligible_clusters.len()
    );

    let ns_unpack = AtomicU64::new(0);
    let ns_format = AtomicU64::new(0);

    // With --start-split auto, the cap signal and the share of a start's reads that must carry
    // it are properties of the library, so phase C runs in two passes. The first finds every
    // cluster's structure groups and what start validation sees in them; the second detects
    // the isoforms from the same groups with the library's decision.
    let mut options = config.isoform_options.clone();
    let groups_by_position: Option<Vec<Vec<Vec<PlacedRead>>>> =
        (options.start_split == StartSplit::Auto).then(|| {
            let pass: Vec<(Vec<Vec<PlacedRead>>, StartStats)> = eligible_clusters
                .par_iter()
                .with_max_len(1)
                .map(|&cid| {
                    let reads = unpack_reads(&cluster_members[cid as usize], &store, &ns_unpack);
                    let groups = isoform::structure_groups(&reads, &options);
                    let stats = isoform::start_stats(&reads, &groups, &options);
                    (groups, stats)
                })
                .collect();
            let mut library_stats = StartStats::default();
            let mut groups_by_position = Vec::with_capacity(pass.len());
            for (groups, stats) in pass {
                library_stats.merge(stats);
                groups_by_position.push(groups);
            }
            let library = library_stats.library(&options);
            if library_stats.candidates.len() >= 100 {
                let (z, auc) = library_stats.cap_sharpness_test(library.threshold);
                eprintln!(
                    "[predict] start validation: capped start candidates against the others, \
                     sharpness AUC {auc:.3} (Mann-Whitney z {z:.1})"
                );
            }
            let rule = if library.cap_signal {
                format!(
                    "on: a start needs {:.0}% of its labelled reads capped ({} start candidates)",
                    100.0 * library.threshold,
                    library_stats.candidates.len()
                )
            } else {
                format!(
                    "off: a start needs {:.0}% of the reads starting within 300 bp to start within \
                     25 bp of it ({} start candidates)",
                    100.0 * library.sharpness,
                    library_stats.sharpness.len()
                )
            };
            eprintln!(
                "[predict] start validation: {} of {} labelled reads capped ({:.1}%) -> cap signal {}",
                library_stats.capped,
                library_stats.labelled,
                100.0 * library_stats.capped as f64 / library_stats.labelled.max(1) as f64,
                rule
            );
            isoform::stats::start_dump(format!(
                "library\t{}\t{}\t{}\t{}\t{}\t{}",
                library_stats.labelled,
                library_stats.capped,
                library.cap_signal as u8,
                library.threshold,
                library_stats.candidates.len(),
                library.sharpness
            ));
            options.start_library = Some(library);
            groups_by_position
        });

    let writer = Mutex::new(OrderedWriter::new(OutputFiles::create(&config.out_dir)));
    let finished = AtomicUsize::new(0);
    // Largest cluster first, each to whichever thread is free: `par_bridge` hands out one
    // cluster at a time, so the largest clusters spread over all threads. (An indexed
    // `par_iter` splits the list into contiguous ranges, and one thread would get the first
    // range: every one of the largest clusters.)
    eligible_clusters.iter().enumerate().par_bridge().for_each(|(position, &cid)| {
        let (records, n_isoforms) = detect_and_format_isoforms(
            cid,
            &cluster_members[cid as usize],
            groups_by_position.as_ref().map(|g| g[position].as_slice()),
            &store,
            &options,
            &ns_unpack,
            &ns_format,
        );
        writer.lock().unwrap().push(position, records, n_isoforms);
        let done = finished.fetch_add(1, Relaxed) + 1;
        if done % 10_000 == 0 {
            eprintln!("[predict]   phase C {} of {} clusters ...", done, eligible_clusters.len());
        }
    });
    let (total_isoforms, n_processed_clusters) = writer.into_inner().unwrap().finish();
    isoform::stats::close_start_dump();
    let phase_c_time = phase_c_start.elapsed();

    // ── Report ──
    let (ns_maps, ns_partition, ns_cons) = isoform::stats::timings_ns();
    let secs = |ns: u64| ns as f64 / 1e9;
    eprintln!(
        "[predict] === phase C took {:.1}s wall  [aggregate CPU across threads: unpack={:.0}s  minimizer-maps={:.0}s  split={:.0}s  consensus={:.0}s  format={:.0}s] ===",
        phase_c_time.as_secs_f64(),
        secs(ns_unpack.load(Relaxed)),
        secs(ns_maps),
        secs(ns_partition),
        secs(ns_cons),
        secs(ns_format.load(Relaxed)),
    );
    isoform::stats::print_report();
    eprintln!(
        "[predict] done: {} reads -> {} clusters -> {} isoforms ({} clusters processed). out: {}/",
        n_reads, n_clusters, total_isoforms, n_processed_clusters, config.out_dir
    );
    eprintln!(
        "[predict] TIMING  A(read+freq)={:.1}s  B(cluster)={:.1}s  C(isoforms)={:.1}s  total={:.1}s",
        phase_a_time.as_secs_f64(),
        phase_b_time.as_secs_f64(),
        phase_c_time.as_secs_f64(),
        run_start.elapsed().as_secs_f64()
    );
}

/// Write every read's phase B cluster, including reads that phase C later leaves out of
/// `isoform_assignments.tsv`.
fn write_clusters(dir: &str, cluster_members: &[Vec<u32>], store: &ReadStore) {
    fs::create_dir_all(dir).expect("mkdir out");
    let mut w = BufWriter::new(File::create(format!("{dir}/gene_clusters.tsv")).unwrap());
    writeln!(w, "read_name\tcluster_id\tcluster_size").unwrap();
    for (cid, member_indices) in cluster_members.iter().enumerate() {
        for &read_idx in member_indices {
            writeln!(w, "{}\t{cid}\t{}", store.read_name(read_idx as usize), member_indices.len())
                .unwrap();
        }
    }
    w.flush().unwrap();
}

/// Phase C output. Clusters finish in any order but are written in `eligible_clusters` order,
/// so the output files do not depend on thread timing: a finished cluster waits here until
/// every cluster before it has been written.
struct OrderedWriter {
    out: OutputFiles,
    /// Finished clusters by position in `eligible_clusters`, with their isoform counts.
    waiting: BTreeMap<usize, (ClusterRecords, usize)>,
    /// Position of the next cluster to write.
    next: usize,
    total_isoforms: usize,
}

impl OrderedWriter {
    fn new(out: OutputFiles) -> Self {
        OrderedWriter { out, waiting: BTreeMap::new(), next: 0, total_isoforms: 0 }
    }

    /// Take the finished cluster at `position`, then write every cluster now next in line.
    fn push(&mut self, position: usize, records: ClusterRecords, n_isoforms: usize) {
        self.waiting.insert(position, (records, n_isoforms));
        while let Some((records, n_isoforms)) = self.waiting.remove(&self.next) {
            self.out.append(&records);
            self.total_isoforms += n_isoforms;
            self.next += 1;
        }
    }

    /// Flush the files. Returns the isoform count and the number of clusters written.
    fn finish(mut self) -> (usize, usize) {
        assert!(self.waiting.is_empty(), "a cluster before position {} never finished", self.next);
        self.out.flush();
        (self.total_isoforms, self.next)
    }
}

/// One cluster's reads, unpacked from the store.
fn unpack_reads(member_indices: &[u32], store: &ReadStore, ns_unpack: &AtomicU64) -> Vec<Read> {
    let lap_start = Instant::now();
    let reads: Vec<Read> = member_indices
        .iter()
        .map(|&read_idx| Read {
            name: store.read_name(read_idx as usize).to_string(),
            seq: store.read_seq(read_idx as usize),
        })
        .collect();
    ns_unpack.fetch_add(lap_start.elapsed().as_nanos() as u64, Relaxed);
    reads
}

/// Phase C for one cluster, on a worker thread: unpack its reads, detect its isoforms (from
/// `groups`, its structure groups, when the first pass found them), and format its output
/// records. Returns the records and the isoform count.
fn detect_and_format_isoforms(
    cid: u32,
    member_indices: &[u32],
    groups: Option<&[Vec<PlacedRead>]>,
    store: &ReadStore,
    options: &IsoformOptions,
    ns_unpack: &AtomicU64,
    ns_format: &AtomicU64,
) -> (ClusterRecords, usize) {
    let reads = unpack_reads(member_indices, store, ns_unpack);
    let isoforms = match groups {
        Some(groups) => isoform::detect_isoforms_in_groups(&reads, groups, options),
        None => isoform::detect_isoforms(&reads, options),
    };

    let lap_start = Instant::now();
    let records = isoform::format_cluster_records(cid, &reads, &isoforms);
    ns_format.fetch_add(lap_start.elapsed().as_nanos() as u64, Relaxed);
    // the isoforms written to isoforms.fasta (the others to the .unvalidated_end files)
    (records, isoforms.iter().filter(|isoform| isoform.end_validated).count())
}
