//! The `predict` subcommand: reads -> clusters -> isoforms in one process.
//!
//! The input file is read exactly once. Every read is kept 2-bit packed in memory
//! ([`ReadStore`]), so clustering hands sequences straight to isoform detection. Three phases:
//!
//!   A. Read and pack every read, counting how many reads contain each seed.
//!   B. Cluster the stored reads greedily, longest first ([`GreedyClusterer`]).
//!   C. Detect the isoforms of each cluster ([`isoform::detect_isoforms`]), in parallel across
//!      clusters, largest first, and write the output files.

use std::cmp::Reverse;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering::Relaxed};
use std::sync::Mutex;
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::cluster::{
    count_seeds, informative_seeds, merge_seed_counts, weighted_seeds, GreedyClusterer, SeedCounts,
    SeedParams,
};
use crate::isoform::{
    self, apply_isoform_flag, ClusterRecords, OutputFiles, Read, IsoformOptions, ISOFORM_HELP,
};
use crate::read_store::ReadStore;
use crate::seedmap::{new_seed_map, SeedMap};
use crate::seqio::SeqReader;

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
    pub min_shared: u32,
    pub min_shared_frac: f64,
    pub level_weights: Vec<u32>,
    // Which clusters to process, by read count.
    pub min_cluster_size: usize,
    pub max_cluster_size: usize,
    // Write each read's phase B cluster to `gene_clusters.tsv`.
    pub dump_clusters: bool,
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
            min_shared: 3,
            min_shared_frac: 0.10,
            level_weights: vec![1, 2, 5],
            min_cluster_size: 1,
            max_cluster_size: usize::MAX,
            dump_clusters: false,
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

CLUSTERING (seeds):
    --k N / --s N           k-mer / s-mer length                         [15 / 9]
    --t N                   open-syncmer offset of the minimal s-mer     [(k-s)/2]
    --levels N              use seeds from levels L1..N                  [3]
    --min-occ N             keep seed if it occurs in >= N reads         [2]
    --max-occ N             ... and in <= N reads                        [100000]
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
    let (store, counts, [d_read, d_seed, d_pack]) = load_reads(&config, seed_params);
    let n_reads = store.len();
    let n_distinct_seeds = counts.len();
    let informative = informative_seeds(&counts, config.min_occ, config.max_occ);
    drop(counts);
    let phase_a_time = run_start.elapsed();
    eprintln!(
        "[predict] {} reads, distinct seeds={}, informative={}",
        n_reads,
        n_distinct_seeds,
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
    let cluster_members = cluster_reads(&config, seed_params, &store, informative);
    let n_clusters = cluster_members.len();
    let phase_b_time = phase_b_start.elapsed();
    eprintln!("[predict] === phase B (clustering) took {:.1}s ===", phase_b_time.as_secs_f64());
    if config.dump_clusters {
        write_clusters(&config.out_dir, &cluster_members, &store);
    }

    // ── Phase C: detect isoforms ──
    let phase_c_start = Instant::now();
    let mut targets: Vec<u32> = (0..n_clusters as u32)
        .filter(|&c| {
            let size = cluster_members[c as usize].len();
            size >= config.min_cluster_size && size <= config.max_cluster_size
        })
        .collect();
    targets.sort_unstable_by_key(|&c| Reverse(cluster_members[c as usize].len()));
    eprintln!(
        "[predict] phase C: detecting isoforms for {} clusters (parallel) ...",
        targets.len()
    );

    let writer = Mutex::new(OrderedWriter::new(OutputFiles::create(&config.out_dir)));
    let ns_unpack = AtomicU64::new(0);
    let ns_format = AtomicU64::new(0);
    let finished = AtomicUsize::new(0);
    // Largest cluster first, each to whichever thread is free: `par_bridge` hands out one
    // cluster at a time, so the largest clusters spread over all threads. (An indexed
    // `par_iter` splits the list into contiguous ranges, and one thread would get the first
    // range: every one of the largest clusters.)
    targets.iter().enumerate().par_bridge().for_each(|(position, &cid)| {
        let (records, n_isoforms) = detect_and_format_isoforms(
            cid,
            &cluster_members[cid as usize],
            &store,
            &config.isoform_options,
            &ns_unpack,
            &ns_format,
        );
        writer.lock().unwrap().push(position, records, n_isoforms);
        let done = finished.fetch_add(1, Relaxed) + 1;
        if done % 10_000 == 0 {
            eprintln!("[predict]   phase C {} of {} clusters ...", done, targets.len());
        }
    });
    let (total_isoforms, n_processed_clusters) = writer.into_inner().unwrap().finish();
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

/// Phase A: pack every read into a [`ReadStore`] while counting seeds, a batch at a time.
/// Also returns the time spent reading, counting, and packing.
fn load_reads(config: &Config, seed_params: SeedParams) -> (ReadStore, SeedCounts, [Duration; 3]) {
    let mut store = ReadStore::default();
    let mut counts: SeedCounts = new_seed_map();
    let (mut d_read, mut d_seed, mut d_pack) = (Duration::ZERO, Duration::ZERO, Duration::ZERO);
    let mut reader = SeqReader::open(&config.reads_path);
    let mut batch: Vec<(String, Vec<u8>)> = Vec::with_capacity(config.batch_size);
    loop {
        batch.clear();
        let remaining = config.max_reads.saturating_sub(store.len());
        if remaining == 0 {
            break;
        }
        let want = config.batch_size.min(remaining);
        let lap_start = Instant::now();
        while batch.len() < want {
            match reader.next() {
                Some(record) => batch.push(record),
                None => break,
            }
        }
        d_read += lap_start.elapsed();
        if batch.is_empty() {
            break;
        }

        let lap_start = Instant::now();
        let part = count_seeds(batch.par_iter().map(|(_, seq)| seq.as_slice()), seed_params);
        merge_seed_counts(&mut counts, part);
        d_seed += lap_start.elapsed();

        let lap_start = Instant::now();
        for (name, seq) in batch.drain(..) {
            store.push(name, &seq);
        }
        d_pack += lap_start.elapsed();
        eprintln!("[predict]   phase A {} reads ...", store.len());
    }
    (store, counts, [d_read, d_seed, d_pack])
}

/// The clustering order: longest read first, ties in file order. A counting sort over read
/// lengths, so visiting the reads in file order yields each one's rank without storing a
/// read -> rank table.
struct LengthRanks {
    /// Read length -> rank of the next read of that length.
    next: BTreeMap<usize, u32>,
}

impl LengthRanks {
    fn new(store: &ReadStore) -> Self {
        let mut next: BTreeMap<usize, u32> = BTreeMap::new();
        for read_idx in 0..store.len() {
            *next.entry(store.read_len(read_idx)).or_insert(0) += 1;
        }
        // Counts -> first rank of each length, longest first.
        let mut rank = 0u32;
        for c in next.values_mut().rev() {
            let count = *c;
            *c = rank;
            rank += count;
        }
        LengthRanks { next }
    }

    /// Rank of the next read of length `len`, visiting reads in file order.
    fn rank(&mut self, len: usize) -> usize {
        let r = self.next.get_mut(&len).expect("length counted in LengthRanks::new");
        *r += 1;
        (*r - 1) as usize
    }
}

/// Phase B: cluster the stored reads, longest first, and return each cluster's member read
/// indices in file order. Seeds are computed in parallel per batch; assignment is sequential.
///
/// The reads stay where they are in the store. One `u32` per read, `slot`, first holds the
/// read at each rank, then, once that read is assigned, its cluster id. So sorting costs
/// no memory beyond the per-read cluster id that clustering needs anyway.
fn cluster_reads(
    config: &Config,
    seed_params: SeedParams,
    store: &ReadStore,
    informative: SeedMap<u8>,
) -> Vec<Vec<u32>> {
    let n_reads = store.len();
    let mut slot: Vec<u32> = vec![0; n_reads];
    let mut ranks = LengthRanks::new(store);
    for read_idx in 0..n_reads {
        slot[ranks.rank(store.read_len(read_idx))] = read_idx as u32;
    }

    let mut clusterer = GreedyClusterer::new(config.min_shared, config.min_shared_frac);
    let mut start = 0;
    while start < n_reads {
        let end = (start + config.batch_size).min(n_reads);
        let seeds: Vec<Vec<(u64, u32)>> = slot[start..end]
            .par_iter()
            .map(|&read_idx| {
                weighted_seeds(
                    &store.read_seq(read_idx as usize),
                    seed_params,
                    &informative,
                    &config.level_weights,
                )
            })
            .collect();
        for (offset, read_seeds) in seeds.iter().enumerate() {
            slot[start + offset] = clusterer.assign(read_seeds);
        }
        eprintln!("[predict]   phase B {} reads, {} clusters ...", end, clusterer.n_clusters());
        start = end;
    }
    let n_clusters = clusterer.n_clusters();
    // Free the seed tables before building the member lists.
    drop(clusterer);
    drop(informative);

    // Replay the ranks in file order to map each read to its cluster id.
    let mut ranks = LengthRanks::new(store);
    let mut cluster_members: Vec<Vec<u32>> = vec![Vec::new(); n_clusters];
    for read_idx in 0..n_reads {
        let c = slot[ranks.rank(store.read_len(read_idx))];
        cluster_members[c as usize].push(read_idx as u32);
    }
    cluster_members
}

/// Phase C output. Clusters finish in any order but are written in `targets` order, so the
/// output files do not depend on thread timing: a finished cluster waits here until every
/// cluster before it has been written.
struct OrderedWriter {
    out: OutputFiles,
    /// Finished clusters by position in `targets`, with their isoform counts.
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

/// Phase C for one cluster, on a worker thread: unpack its reads, detect its isoforms, and
/// format its output records. Returns the records and the isoform count.
fn detect_and_format_isoforms(
    cid: u32,
    member_indices: &[u32],
    store: &ReadStore,
    options: &IsoformOptions,
    ns_unpack: &AtomicU64,
    ns_format: &AtomicU64,
) -> (ClusterRecords, usize) {
    let lap_start = Instant::now();
    let reads: Vec<Read> = member_indices
        .iter()
        .map(|&read_idx| Read {
            name: store.read_name(read_idx as usize).to_string(),
            seq: store.read_seq(read_idx as usize),
        })
        .collect();
    ns_unpack.fetch_add(lap_start.elapsed().as_nanos() as u64, Relaxed);

    let isoforms = isoform::detect_isoforms(&reads, options);

    let lap_start = Instant::now();
    let records = isoform::format_cluster_records(cid, &reads, &isoforms);
    ns_format.fetch_add(lap_start.elapsed().as_nanos() as u64, Relaxed);
    (records, isoforms.len())
}
