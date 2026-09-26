//! The `predict` subcommand: reads -> clusters -> isoforms in one process.
//!
//! The input file is read exactly once. Every read is kept 2-bit packed in memory
//! ([`SeqStore`]), so clustering hands sequences straight to isoform resolution. Three phases:
//!
//!   A. Read and pack every read, counting how many reads contain each seed.
//!   B. Cluster the stored reads greedily, longest first ([`GreedyClusterer`]).
//!   C. Resolve each cluster into isoforms ([`isoform::resolve_cluster`]), in parallel across
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
use crate::hashmap::{fmap, FastMap};
use crate::io::SeqReader;
use crate::isoform::{self, parse_flag, Cfg, OutputFiles, Read, Rendered, RESOLVE_HELP};
use crate::seqstore::SeqStore;

pub struct Config {
    pub reads: String,
    pub out_dir: String,
    pub threads: usize,
    pub batch: usize,
    pub max_reads: usize,
    // Clustering.
    pub k: usize,
    pub s: usize,
    pub t: usize,
    pub levels: usize,
    pub min_occ: u32,
    pub max_occ: u32,
    pub min_shared: u32,
    pub min_shared_frac: f64,
    pub weights: Vec<u32>,
    // Which clusters to resolve, by read count.
    pub min_size: usize,
    pub max_size: usize,
    // Write each read's phase B cluster to `gene_clusters.tsv`.
    pub dump_clusters: bool,
    // Isoform resolution.
    pub cfg: Cfg,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            reads: String::new(),
            out_dir: "lcpsynciso_out".into(),
            threads: 0,
            batch: 500_000,
            max_reads: usize::MAX,
            k: 15,
            s: 9,
            t: usize::MAX,
            levels: 3,
            min_occ: 2,
            max_occ: 100_000,
            min_shared: 3,
            min_shared_frac: 0.10,
            weights: vec![1, 2, 5],
            min_size: 1,
            max_size: usize::MAX,
            dump_clusters: false,
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
            "-o" | "--out" => a.out_dir = next(),
            "--threads" => a.threads = next().parse().unwrap(),
            "--batch" => a.batch = next().parse().unwrap(),
            "--max-reads" => a.max_reads = next().parse().unwrap(),
            "--k" => a.k = next().parse().unwrap(),
            "--s" => a.s = next().parse().unwrap(),
            "--t" => a.t = next().parse().unwrap(),
            "--levels" => a.levels = next().parse().unwrap(),
            "--min-occ" => a.min_occ = next().parse().unwrap(),
            "--max-occ" => a.max_occ = next().parse().unwrap(),
            "--min-shared" => a.min_shared = next().parse().unwrap(),
            "--min-shared-frac" => a.min_shared_frac = next().parse().unwrap(),
            "--level-weights" => {
                a.weights = next().split(',').map(|x| x.trim().parse().unwrap()).collect()
            }
            "--min-size" => a.min_size = next().parse().unwrap(),
            "--max-size" => a.max_size = next().parse().unwrap(),
            "--iso-k" => a.cfg.k = next().parse().unwrap(),
            "--iso-w" => a.cfg.w = next().parse().unwrap(),
            "--no-consensus" => a.cfg.consensus = false,
            "--dump-clusters" => a.dump_clusters = true,
            "-h" | "--help" => {
                usage();
                std::process::exit(0);
            }
            other => {
                if !parse_flag(&mut a.cfg, other, &mut next) {
                    if other.starts_with('-') {
                        eprintln!("unknown flag {other}");
                        std::process::exit(2);
                    }
                    a.reads = other.into();
                }
            }
        }
        i += 1;
    }
    if a.reads.is_empty() {
        usage();
        std::process::exit(2);
    }
    if a.t == usize::MAX {
        a.t = (a.k - a.s) / 2;
    }
    a
}

pub fn usage() {
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
    --no-consensus          write each isoform's longest read instead of a consensus

{RESOLVE_HELP}

OUTPUT (in DIR):
    isoform_assignments.tsv   read_name <tab> cluster_id <tab> isoform_id
    isoform_summary.tsv       cluster_id <tab> isoform_id <tab> n_reads <tab> length
    isoforms.fasta            refined consensus per isoform
    gene_clusters.tsv         read_name <tab> cluster_id <tab> cluster_size (--dump-clusters)"
    );
}

pub fn run(a: Config) {
    if a.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(a.threads).build_global().ok();
    }
    let params = SeedParams { k: a.k, s: a.s, t: a.t, levels: a.levels };
    let t_start = Instant::now();

    // ── Phase A: read, pack, count seeds ──
    eprintln!("[predict] phase A: reading {} once (pack + seed frequencies) ...", a.reads);
    let (store, counts, [d_read, d_seed, d_pack]) = load_reads(&a, params);
    let n = store.len();
    let total_distinct = counts.len();
    let informative = informative_seeds(&counts, a.min_occ, a.max_occ);
    drop(counts);
    let t_a = t_start.elapsed();
    eprintln!(
        "[predict] {} reads, distinct seeds={}, informative={}",
        n,
        total_distinct,
        informative.len()
    );
    eprintln!(
        "[predict] === phase A took {:.1}s  [read/decompress={:.1}s  seed-freq={:.1}s  2bit-pack={:.1}s] ===",
        t_a.as_secs_f64(),
        d_read.as_secs_f64(),
        d_seed.as_secs_f64(),
        d_pack.as_secs_f64()
    );

    // ── Phase B: cluster ──
    let t_b0 = Instant::now();
    eprintln!(
        "[predict] phase B: greedy clustering (min-shared={}, min-shared-frac={}) ...",
        a.min_shared, a.min_shared_frac
    );
    let members = cluster_reads(&a, params, &store, informative);
    let n_clusters = members.len();
    let t_b = t_b0.elapsed();
    eprintln!("[predict] === phase B (clustering) took {:.1}s ===", t_b.as_secs_f64());
    if a.dump_clusters {
        write_clusters(&a.out_dir, &members, &store);
    }

    // ── Phase C: resolve isoforms ──
    let t_c0 = Instant::now();
    let mut targets: Vec<u32> = (0..n_clusters as u32)
        .filter(|&c| {
            let size = members[c as usize].len();
            size >= a.min_size && size <= a.max_size
        })
        .collect();
    targets.sort_unstable_by_key(|&c| Reverse(members[c as usize].len()));
    eprintln!(
        "[predict] phase C: resolving isoforms for {} clusters (parallel) ...",
        targets.len()
    );

    let writer = Mutex::new(OrderedWriter::new(OutputFiles::create(&a.out_dir)));
    let t_unpack = AtomicU64::new(0);
    let t_render = AtomicU64::new(0);
    let finished = AtomicUsize::new(0);
    // Largest cluster first, each to whichever thread is free: `par_bridge` hands out one
    // cluster at a time, so the largest clusters spread over all threads. (An indexed
    // `par_iter` splits the list into contiguous ranges, and one thread would get the first
    // range: every one of the largest clusters.)
    targets.iter().enumerate().par_bridge().for_each(|(k, &cid)| {
        let (records, n_iso) =
            resolve_one(cid, &members[cid as usize], &store, &a.cfg, &t_unpack, &t_render);
        writer.lock().unwrap().push(k, records, n_iso);
        let done = finished.fetch_add(1, Relaxed) + 1;
        if done % 10_000 == 0 {
            eprintln!("[predict]   phase C {} of {} clusters ...", done, targets.len());
        }
    });
    let (total_iso, n_resolved) = writer.into_inner().unwrap().finish();
    let t_c = t_c0.elapsed();

    // ── Report ──
    let (ns_maps, ns_split, ns_cons) = isoform::stats::timings_ns();
    let secs = |ns: u64| ns as f64 / 1e9;
    eprintln!(
        "[predict] === phase C took {:.1}s wall  [aggregate CPU across threads: unpack={:.0}s  minimizer-maps={:.0}s  split={:.0}s  consensus={:.0}s  render={:.0}s] ===",
        t_c.as_secs_f64(),
        secs(t_unpack.load(Relaxed)),
        secs(ns_maps),
        secs(ns_split),
        secs(ns_cons),
        secs(t_render.load(Relaxed)),
    );
    isoform::stats::print_report();
    eprintln!(
        "[predict] done: {} reads -> {} clusters -> {} isoforms ({} clusters resolved). out: {}/",
        n, n_clusters, total_iso, n_resolved, a.out_dir
    );
    eprintln!(
        "[predict] TIMING  A(read+freq)={:.1}s  B(cluster)={:.1}s  C(isoforms)={:.1}s  total={:.1}s",
        t_a.as_secs_f64(),
        t_b.as_secs_f64(),
        t_c.as_secs_f64(),
        t_start.elapsed().as_secs_f64()
    );
}

/// Write every read's phase B cluster, including reads that phase C later leaves out of
/// `isoform_assignments.tsv`.
fn write_clusters(dir: &str, members: &[Vec<u32>], store: &SeqStore) {
    fs::create_dir_all(dir).expect("mkdir out");
    let mut w = BufWriter::new(File::create(format!("{dir}/gene_clusters.tsv")).unwrap());
    writeln!(w, "read_name\tcluster_id\tcluster_size").unwrap();
    for (cid, reads) in members.iter().enumerate() {
        for &ri in reads {
            writeln!(w, "{}\t{cid}\t{}", store.name(ri as usize), reads.len()).unwrap();
        }
    }
    w.flush().unwrap();
}

/// Phase A: pack every read into a [`SeqStore`] while counting seeds, a batch at a time.
/// Also returns the time spent reading, counting, and packing.
fn load_reads(a: &Config, params: SeedParams) -> (SeqStore, SeedCounts, [Duration; 3]) {
    let mut store = SeqStore::default();
    let mut counts: SeedCounts = fmap();
    let (mut d_read, mut d_seed, mut d_pack) = (Duration::ZERO, Duration::ZERO, Duration::ZERO);
    let mut reader = SeqReader::open(&a.reads);
    let mut batch: Vec<(String, Vec<u8>)> = Vec::with_capacity(a.batch);
    loop {
        batch.clear();
        let remaining = a.max_reads.saturating_sub(store.len());
        if remaining == 0 {
            break;
        }
        let want = a.batch.min(remaining);
        let t = Instant::now();
        while batch.len() < want {
            match reader.next() {
                Some(record) => batch.push(record),
                None => break,
            }
        }
        d_read += t.elapsed();
        if batch.is_empty() {
            break;
        }

        let t = Instant::now();
        let part = count_seeds(batch.par_iter().map(|(_, seq)| seq.as_slice()), params);
        merge_seed_counts(&mut counts, part);
        d_seed += t.elapsed();

        let t = Instant::now();
        for (name, seq) in batch.drain(..) {
            store.push(name, &seq);
        }
        d_pack += t.elapsed();
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
    fn new(store: &SeqStore) -> Self {
        let mut next: BTreeMap<usize, u32> = BTreeMap::new();
        for i in 0..store.len() {
            *next.entry(store.read_len(i)).or_insert(0) += 1;
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
    a: &Config,
    params: SeedParams,
    store: &SeqStore,
    informative: FastMap<u8>,
) -> Vec<Vec<u32>> {
    let n = store.len();
    let mut slot: Vec<u32> = vec![0; n];
    let mut ranks = LengthRanks::new(store);
    for i in 0..n {
        slot[ranks.rank(store.read_len(i))] = i as u32;
    }

    let mut clusterer = GreedyClusterer::new(a.min_shared, a.min_shared_frac);
    let mut start = 0;
    while start < n {
        let end = (start + a.batch).min(n);
        let seeds: Vec<Vec<(u64, u32)>> = slot[start..end]
            .par_iter()
            .map(|&i| weighted_seeds(&store.get(i as usize), params, &informative, &a.weights))
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
    let mut members: Vec<Vec<u32>> = vec![Vec::new(); n_clusters];
    for i in 0..n {
        let c = slot[ranks.rank(store.read_len(i))];
        members[c as usize].push(i as u32);
    }
    members
}

/// Phase C output. Clusters finish in any order but are written in `targets` order, so the
/// output files do not depend on thread timing: a finished cluster waits here until every
/// cluster before it has been written.
struct OrderedWriter {
    out: OutputFiles,
    /// Finished clusters by position in `targets`, with their isoform counts.
    waiting: BTreeMap<usize, (Rendered, usize)>,
    /// Position of the next cluster to write.
    next: usize,
    total_iso: usize,
}

impl OrderedWriter {
    fn new(out: OutputFiles) -> Self {
        OrderedWriter { out, waiting: BTreeMap::new(), next: 0, total_iso: 0 }
    }

    /// Take the finished cluster at position `k`, then write every cluster now next in line.
    fn push(&mut self, k: usize, records: Rendered, n_iso: usize) {
        self.waiting.insert(k, (records, n_iso));
        while let Some((records, n_iso)) = self.waiting.remove(&self.next) {
            self.out.append(&records);
            self.total_iso += n_iso;
            self.next += 1;
        }
    }

    /// Flush the files. Returns the isoform count and the number of clusters written.
    fn finish(mut self) -> (usize, usize) {
        assert!(self.waiting.is_empty(), "a cluster before position {} never finished", self.next);
        self.out.flush();
        (self.total_iso, self.next)
    }
}

/// Phase C for one cluster, on a worker thread: unpack its reads, resolve its isoforms, and
/// render its output records. Returns the records and the isoform count.
fn resolve_one(
    cid: u32,
    member_ids: &[u32],
    store: &SeqStore,
    cfg: &Cfg,
    t_unpack: &AtomicU64,
    t_render: &AtomicU64,
) -> (Rendered, usize) {
    let t = Instant::now();
    let reads: Vec<Read> = member_ids
        .iter()
        .map(|&i| Read { name: store.name(i as usize).to_string(), seq: store.get(i as usize) })
        .collect();
    t_unpack.fetch_add(t.elapsed().as_nanos() as u64, Relaxed);

    let isos = isoform::resolve_cluster(&reads, cfg);

    let t = Instant::now();
    let records = isoform::render(cid, &reads, &isos);
    t_render.fetch_add(t.elapsed().as_nanos() as u64, Relaxed);
    (records, isos.len())
}
