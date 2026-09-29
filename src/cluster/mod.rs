//! Greedy single-pass clustering of reads by shared seeds: the `cluster` subcommand, and
//! the phases A and B that `predict` shares with it.
//!
//! Seeds are synpact's LCP block hashes from levels L1..L`levels`. Only mid-frequency seeds
//! are informative: found in at least `min_occ` reads (shared) and at most `max_occ`
//! (not a repeat). Reads are then clustered one at a time, longest first (ties in file
//! order), so each cluster is seeded by its most complete read:
//!
//!   for each read:
//!     each informative seed already claimed by a cluster votes for that cluster,
//!       weighted by its level (higher levels are more specific)
//!     if the best cluster's votes reach min_shared, and also min_shared_frac of
//!       the read's total seed weight                 -> join it
//!     otherwise                                      -> start a new cluster
//!     the read's unclaimed seeds are claimed by its cluster
//!
//! A read of the same gene shares most of its seeds with that gene's cluster; a read of
//! another gene may share a few by chance, such as a repeat in its UTR. The fraction keeps
//! those few from merging unrelated genes, which a fixed `min_shared` cannot do for long
//! reads without also splitting reads that have few seeds. Longest-first order needs it,
//! because the longest reads start clusters first.
//!
//! The file is read once. Phase A (`load_reads`) packs every read into a [`ReadStore`]
//! while counting how many reads contain each seed; phase B (`cluster_reads`) clusters the
//! stored reads. The `cluster` subcommand stops there; `predict` goes on to detect isoforms.
//! Seeds are computed in parallel per batch; assignment is sequential because it depends on
//! read order.

use std::collections::BTreeMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::read_store::ReadStore;
use crate::seedmap::{new_seed_map, SeedMap};
use crate::seqio::SeqReader;
use crate::seeds::read_seed_levels;

// ── Building blocks shared with `predict` ────────────────────────────────────

/// Seed hash -> (number of reads containing it, level).
pub(crate) type SeedCounts = SeedMap<(u32, u8)>;

/// Seed extraction settings.
#[derive(Clone, Copy)]
pub(crate) struct SeedParams {
    pub k: usize,
    pub s: usize,
    /// Open-syncmer offset of the minimal s-mer.
    pub t: usize,
    pub levels: usize,
}

/// Count, for each seed, how many of `seqs` contain it.
pub(crate) fn count_seeds<'a>(
    seqs: impl ParallelIterator<Item = &'a [u8]>,
    params: SeedParams,
) -> SeedCounts {
    seqs.fold(new_seed_map, |mut m, seq| {
        for (h, level) in read_seed_levels(seq, params.k, params.s, params.t, params.levels) {
            m.entry(h).or_insert((0, level)).0 += 1;
        }
        m
    })
    .reduce(new_seed_map, |mut acc, b| {
        if acc.is_empty() {
            return b;
        }
        merge_seed_counts(&mut acc, b);
        acc
    })
}

/// Add `part`'s counts into `total`.
pub(crate) fn merge_seed_counts(total: &mut SeedCounts, part: SeedCounts) {
    if total.is_empty() {
        *total = part;
        return;
    }
    for (h, (c, level)) in part {
        total.entry(h).or_insert((0, level)).0 += c;
    }
}

/// Seeds found in between `min_occ` and `max_occ` reads, as hash -> level.
pub(crate) fn informative_seeds(counts: &SeedCounts, min_occ: u32, max_occ: u32) -> SeedMap<u8> {
    let mut informative: SeedMap<u8> = new_seed_map();
    for (h, (c, level)) in counts.iter() {
        if *c >= min_occ && *c <= max_occ {
            informative.insert(*h, *level);
        }
    }
    informative
}

/// A read's informative seeds with their vote weights. `level_weights[i]` is the weight of
/// level `i + 1`; levels beyond it weigh their level number.
pub(crate) fn weighted_seeds(
    seq: &[u8],
    params: SeedParams,
    informative: &SeedMap<u8>,
    level_weights: &[u32],
) -> Vec<(u64, u32)> {
    read_seed_levels(seq, params.k, params.s, params.t, params.levels)
        .into_iter()
        .filter_map(|(h, _)| {
            informative.get(&h).map(|&level| {
                (h, level_weights.get(level as usize - 1).copied().unwrap_or(level as u32))
            })
        })
        .collect()
}

/// Assigns reads to clusters one at a time, in the order given (see the module docs).
pub(crate) struct GreedyClusterer {
    min_shared: u32,
    /// Fraction of the read's total seed weight the votes must also reach.
    min_shared_frac: f64,
    /// Seed hash -> the cluster that claimed it first.
    owner: SeedMap<u32>,
    /// Scratch for one read: cluster id -> total vote weight.
    votes: SeedMap<u32>,
    /// Reads per cluster.
    cluster_sizes: Vec<u32>,
}

impl GreedyClusterer {
    pub(crate) fn new(min_shared: u32, min_shared_frac: f64) -> Self {
        GreedyClusterer {
            min_shared,
            min_shared_frac,
            owner: new_seed_map(),
            votes: new_seed_map(),
            cluster_sizes: Vec::new(),
        }
    }

    /// Assign a read, given its weighted informative seeds, and return its cluster id.
    pub(crate) fn assign(&mut self, seeds: &[(u64, u32)]) -> u32 {
        self.votes.clear();
        let total: u32 = seeds.iter().map(|&(_, weight)| weight).sum();
        let need = self.min_shared.max((self.min_shared_frac * total as f64).ceil() as u32);
        let mut best = u32::MAX;
        let mut best_votes = 0u32;
        for (h, weight) in seeds {
            if let Some(&cid) = self.owner.get(h) {
                let e = self.votes.entry(cid as u64).or_insert(0);
                *e += *weight;
                if *e > best_votes {
                    best_votes = *e;
                    best = cid;
                }
            }
        }
        let cid = if best_votes >= need && best != u32::MAX {
            best
        } else {
            let new = self.cluster_sizes.len() as u32;
            self.cluster_sizes.push(0);
            new
        };
        self.cluster_sizes[cid as usize] += 1;
        for (h, _) in seeds {
            self.owner.entry(*h).or_insert(cid);
        }
        cid
    }

    pub(crate) fn n_clusters(&self) -> usize {
        self.cluster_sizes.len()
    }
}

// ── Phases A and B, shared by `cluster` and `predict` ────────────────────────

/// Phase A: pack every read into a [`ReadStore`] while counting seeds, a batch at a time.
/// `command` names the subcommand in progress messages. Also returns the time spent reading,
/// counting, and packing.
pub(crate) fn load_reads(
    reads_path: &str,
    batch_size: usize,
    max_reads: usize,
    seed_params: SeedParams,
    command: &str,
) -> (ReadStore, SeedCounts, [Duration; 3]) {
    let mut store = ReadStore::default();
    let mut counts: SeedCounts = new_seed_map();
    let (mut d_read, mut d_seed, mut d_pack) = (Duration::ZERO, Duration::ZERO, Duration::ZERO);
    let mut reader = SeqReader::open(reads_path);
    let mut batch: Vec<(String, Vec<u8>)> = Vec::with_capacity(batch_size);
    loop {
        batch.clear();
        let remaining = max_reads.saturating_sub(store.len());
        if remaining == 0 {
            break;
        }
        let want = batch_size.min(remaining);
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
        eprintln!("[{command}]   phase A {} reads ...", store.len());
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

/// Phase B: cluster the stored reads with `clusterer`, longest first, and return each
/// cluster's member read indices in file order. Seeds are computed in parallel per batch;
/// assignment is sequential. `command` names the subcommand in progress messages.
///
/// The reads stay where they are in the store. One `u32` per read, `slot`, first holds the
/// read at each rank, then, once that read is assigned, its cluster id. So sorting costs
/// no memory beyond the per-read cluster id that clustering needs anyway.
pub(crate) fn cluster_reads(
    store: &ReadStore,
    seed_params: SeedParams,
    informative: SeedMap<u8>,
    level_weights: &[u32],
    mut clusterer: GreedyClusterer,
    batch_size: usize,
    command: &str,
) -> Vec<Vec<u32>> {
    let n_reads = store.len();
    let mut slot: Vec<u32> = vec![0; n_reads];
    let mut ranks = LengthRanks::new(store);
    for read_idx in 0..n_reads {
        slot[ranks.rank(store.read_len(read_idx))] = read_idx as u32;
    }

    let mut start = 0;
    while start < n_reads {
        let end = (start + batch_size).min(n_reads);
        let seeds: Vec<Vec<(u64, u32)>> = slot[start..end]
            .par_iter()
            .map(|&read_idx| {
                weighted_seeds(
                    &store.read_seq(read_idx as usize),
                    seed_params,
                    &informative,
                    level_weights,
                )
            })
            .collect();
        for (offset, read_seeds) in seeds.iter().enumerate() {
            slot[start + offset] = clusterer.assign(read_seeds);
        }
        eprintln!("[{command}]   phase B {} reads, {} clusters ...", end, clusterer.n_clusters());
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

// ── The `cluster` subcommand ─────────────────────────────────────────────────

pub struct Config {
    pub reads_path: String,
    pub seed_k: usize,
    pub seed_s: usize,
    pub seed_t: usize,
    pub levels: usize,
    pub min_occ: u32,
    pub max_occ: u32,
    /// Total vote weight needed to join a cluster.
    pub min_shared: u32,
    /// ... and the fraction of the read's total seed weight it must also reach.
    pub min_shared_frac: f64,
    /// Vote weight per level (index 0 = L1).
    pub level_weights: Vec<u32>,
    pub batch_size: usize,
    pub max_reads: usize,
    pub threads: usize,
    pub out_dir: String,
    pub emit_assignments: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            reads_path: String::new(),
            seed_k: 15,
            seed_s: 9,
            seed_t: usize::MAX,
            levels: 3,
            min_occ: 2,
            max_occ: 100_000,
            min_shared: 3,
            min_shared_frac: 0.10,
            level_weights: vec![1, 2, 5],
            batch_size: 500_000,
            max_reads: usize::MAX,
            threads: 0,
            out_dir: "lcpsynciso_cluster_out".into(),
            emit_assignments: false,
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
            "--batch" => config.batch_size = next().parse().unwrap(),
            "--max-reads" => config.max_reads = next().parse().unwrap(),
            "--threads" => config.threads = next().parse().unwrap(),
            "--emit-assignments" => config.emit_assignments = true,
            "-o" | "--out" => config.out_dir = next(),
            "-h" | "--help" => {
                print_usage();
                std::process::exit(0);
            }
            other => {
                if other.starts_with('-') {
                    eprintln!("unknown flag {other}");
                    std::process::exit(2);
                }
                config.reads_path = other.into();
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
        "lcpsynciso cluster — greedy clustering of reads by shared, level-weighted seeds
(longest read first; the file is read once and the reads are held in RAM, 2-bit packed)

USAGE:
    lcpsynciso cluster <reads.fq[.gz]|.fa[.gz]> [options]

OPTIONS:
    --k N            k-mer length                              [15]
    --s N            s-mer length                              [9]
    --levels N       use seeds from levels L1..N               [3]
    --min-occ N      keep seed if it occurs in >= N reads      [2]
    --max-occ N      ... and in <= N reads                     [100000]
    --level-weights a,b,c   per-level vote weights (L1,L2,..)  [1,2,5]
    --min-shared N   join a cluster if weighted shared >= N    [3]
    --min-shared-frac F  ... and >= F of the read's seed weight [0.1]
    --batch N        reads per batch                           [500000]
    --max-reads N    stop after N reads                        [all]
    --threads N      worker threads (0 = all)                  [0]
    --emit-assignments  also write read_name -> cluster TSV
    -o, --out DIR    output directory                          [lcpsynciso_cluster_out]"
    );
}

pub fn run(config: Config) {
    if config.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(config.threads).build_global().ok();
    }
    let params =
        SeedParams { k: config.seed_k, s: config.seed_s, t: config.seed_t, levels: config.levels };

    // ── Phase A: read, pack, count seeds ──
    eprintln!(
        "[cluster] phase A: reading {} once (pack + seed frequencies; k={} s={} L1..{}) ...",
        config.reads_path, config.seed_k, config.seed_s, config.levels
    );
    let (store, counts, _) =
        load_reads(&config.reads_path, config.batch_size, config.max_reads, params, "cluster");
    let n_reads = store.len() as u64;
    let n_distinct_seeds = counts.len();
    let informative = informative_seeds(&counts, config.min_occ, config.max_occ);
    drop(counts);
    let n_informative_seeds = informative.len();
    eprintln!(
        "[cluster] distinct seeds={n_distinct_seeds}  informative (in [{},{}])={n_informative_seeds}",
        config.min_occ, config.max_occ
    );

    // ── Phase B: greedy clustering, longest read first ──
    eprintln!(
        "[cluster] phase B: greedy clustering (weights={:?}, min-shared={}, min-shared-frac={}) ...",
        config.level_weights, config.min_shared, config.min_shared_frac
    );
    let clusterer = GreedyClusterer::new(config.min_shared, config.min_shared_frac);
    let cluster_members = cluster_reads(
        &store,
        params,
        informative,
        &config.level_weights,
        clusterer,
        config.batch_size,
        "cluster",
    );

    if config.emit_assignments {
        write_assignments(&config.out_dir, &store, &cluster_members);
    }
    let cluster_sizes: Vec<u32> = cluster_members.iter().map(|m| m.len() as u32).collect();
    write_outputs(&config, n_reads, &cluster_sizes, n_distinct_seeds, n_informative_seeds);
}

/// Write `assignments.tsv`: every read's cluster, in file order.
fn write_assignments(dir: &str, store: &ReadStore, cluster_members: &[Vec<u32>]) {
    let mut cluster_of: Vec<u32> = vec![0; store.len()];
    for (cid, member_indices) in cluster_members.iter().enumerate() {
        for &read_idx in member_indices {
            cluster_of[read_idx as usize] = cid as u32;
        }
    }
    fs::create_dir_all(dir).expect("mkdir out");
    let mut w = BufWriter::new(
        fs::File::create(format!("{dir}/assignments.tsv")).expect("create assignments.tsv"),
    );
    writeln!(w, "read_name\tcluster_id").unwrap();
    for (read_idx, cid) in cluster_of.iter().enumerate() {
        writeln!(w, "{}\t{cid}", store.read_name(read_idx)).unwrap();
    }
    w.flush().unwrap();
}

/// Write summary.tsv and cluster_size_hist.tsv, and print a summary to stdout.
fn write_outputs(
    config: &Config,
    n_reads: u64,
    cluster_sizes: &[u32],
    n_distinct_seeds: usize,
    n_informative_seeds: usize,
) {
    fs::create_dir_all(&config.out_dir).expect("mkdir out");
    let n_clusters = cluster_sizes.len();
    let mut singletons = 0u64;
    let mut max_sz = 0u32;
    let mut reads_in_clusters: u64 = 0;
    let mut hist: BTreeMap<u32, u64> = BTreeMap::new();
    for &sz in cluster_sizes {
        if sz == 1 {
            singletons += 1;
        }
        if sz > max_sz {
            max_sz = sz;
        }
        reads_in_clusters += sz as u64;
        *hist.entry(sz).or_insert(0) += 1;
    }
    let mut body = String::from("cluster_size\tnum_clusters\n");
    for (sz, n) in &hist {
        body.push_str(&format!("{sz}\t{n}\n"));
    }
    fs::write(format!("{}/cluster_size_hist.tsv", config.out_dir), body).unwrap();

    let mut sizes: Vec<u32> = cluster_sizes.to_vec();
    sizes.sort_unstable_by(|x, y| y.cmp(x));
    let top: Vec<u32> = sizes.iter().take(20).copied().collect();
    let reads_in_non_singleton_clusters = reads_in_clusters - singletons;
    let percent_reads_in_non_singleton_clusters =
        100.0 * reads_in_non_singleton_clusters as f64 / n_reads.max(1) as f64;
    let mean_non_singleton_cluster_size: f64 = {
        let non_singleton_sizes: Vec<u32> =
            cluster_sizes.iter().copied().filter(|&size| size > 1).collect();
        if non_singleton_sizes.is_empty() {
            0.0
        } else {
            let total_reads = non_singleton_sizes.iter().map(|&size| size as u64).sum::<u64>();
            total_reads as f64 / non_singleton_sizes.len() as f64
        }
    };

    let mut summary = String::new();
    summary.push_str(&format!("reads\t{n_reads}\n"));
    summary.push_str(&format!(
        "k\t{}\ns\t{}\nlevels\t{}\n",
        config.seed_k, config.seed_s, config.levels
    ));
    summary.push_str(&format!(
        "min_occ\t{}\nmax_occ\t{}\nmin_shared\t{}\nmin_shared_frac\t{}\n",
        config.min_occ, config.max_occ, config.min_shared, config.min_shared_frac
    ));
    summary.push_str(&format!(
        "level_weights\t{}\n",
        config.level_weights.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
    ));
    summary.push_str(&format!("distinct_seeds\t{n_distinct_seeds}\n"));
    summary.push_str(&format!("informative_seeds\t{n_informative_seeds}\n"));
    summary.push_str(&format!("clusters\t{n_clusters}\n"));
    summary.push_str(&format!("singleton_clusters\t{singletons}\n"));
    summary.push_str(&format!("multi_read_clusters\t{}\n", n_clusters as u64 - singletons));
    summary.push_str(&format!("largest_cluster\t{max_sz}\n"));
    summary.push_str(&format!("mean_size_multi_clusters\t{mean_non_singleton_cluster_size:.2}\n"));
    summary.push_str(&format!(
        "pct_reads_in_multi_clusters\t{percent_reads_in_non_singleton_clusters:.2}\n"
    ));
    summary.push_str(&format!(
        "top20_sizes\t{}\n",
        top.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
    ));
    fs::write(format!("{}/summary.tsv", config.out_dir), &summary).unwrap();

    println!();
    println!(
        "  reads={n_reads}  k={} s={} L1..{}  band=[{},{}]  weights={:?}  min_shared={}  min_shared_frac={}",
        config.seed_k,
        config.seed_s,
        config.levels,
        config.min_occ,
        config.max_occ,
        config.level_weights,
        config.min_shared,
        config.min_shared_frac
    );
    println!("  distinct seeds={n_distinct_seeds}  informative={n_informative_seeds}");
    println!(
        "  clusters={n_clusters}  (singletons={singletons}, multi={})",
        n_clusters as u64 - singletons
    );
    println!("  reads in multi-read clusters: {reads_in_non_singleton_clusters} ({percent_reads_in_non_singleton_clusters:.1}%)");
    println!("  largest cluster={max_sz}  mean(multi)={mean_non_singleton_cluster_size:.1}");
    println!(
        "  top sizes: {}",
        top.iter().take(12).map(|x| x.to_string()).collect::<Vec<_>>().join(", ")
    );
    println!();
    println!("  wrote summary.tsv + cluster_size_hist.tsv to {}/", config.out_dir);
}

#[cfg(test)]
mod tests {
    use rayon::prelude::*;

    use super::{
        cluster_reads, count_seeds, informative_seeds, GreedyClusterer, LengthRanks, SeedParams,
    };
    use crate::read_store::ReadStore;

    /// A fixed pseudo-random sequence, so seeds are unique and well spread.
    fn random_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    /// Twenty unit-weight seeds, the first `shared` of which the first read also has.
    fn read(shared: u64) -> Vec<(u64, u32)> {
        (0..20).map(|i| (if i < shared { i } else { 1000 + i }, 1)).collect()
    }

    #[test]
    fn min_shared_frac_needs_a_share_of_the_reads_own_seeds() {
        // (min_shared_frac, seeds shared with cluster 0, joins cluster 0)
        for (frac, shared, joins) in [(0.0, 3, true), (0.25, 3, false), (0.25, 5, true)] {
            let mut c = GreedyClusterer::new(3, frac);
            assert_eq!(c.assign(&read(20)), 0);
            assert_eq!(c.assign(&read(shared)) == 0, joins, "frac={frac} shared={shared}");
        }
    }

    #[test]
    fn length_ranks_put_the_longest_first_and_keep_file_order_among_equals() {
        let mut store = ReadStore::default();
        for (i, len) in [5usize, 9, 5, 9, 7].into_iter().enumerate() {
            store.push(format!("r{i}"), &vec![b'A'; len]);
        }
        let mut ranks = LengthRanks::new(&store);
        let ranked: Vec<usize> = (0..store.len()).map(|i| ranks.rank(store.read_len(i))).collect();
        // Lengths 9, 9, 7, 5, 5 take ranks 0..5; equal lengths keep file order.
        assert_eq!(ranked, vec![3, 0, 4, 1, 2]);
    }

    #[test]
    fn the_longest_read_seeds_the_first_cluster_wherever_it_sits_in_the_file() {
        // Two genes, three overlapping copies each. The longest read of all, the full-length
        // copy of gene 2, is last in the file, so cluster 0 belongs to gene 2 only if the
        // reads are clustered longest first; in file order it would belong to gene 1.
        let gene1 = random_seq(1200, 1);
        let gene2 = random_seq(1800, 2);
        let seqs = vec![
            gene1.clone(),
            gene1[100..].to_vec(),
            gene1[200..].to_vec(),
            gene2[300..].to_vec(),
            gene2[600..].to_vec(),
            gene2.clone(),
        ];
        let mut store = ReadStore::default();
        for (i, seq) in seqs.iter().enumerate() {
            store.push(format!("r{i}"), seq);
        }
        let params = SeedParams { k: 15, s: 9, t: 3, levels: 3 };
        let counts = count_seeds(seqs.par_iter().map(|seq| seq.as_slice()), params);
        let informative = informative_seeds(&counts, 2, 100_000);

        // A batch size of 4 splits the six reads over two batches.
        let clusters = cluster_reads(
            &store,
            params,
            informative,
            &[1, 2, 5],
            GreedyClusterer::new(3, 0.1),
            4,
            "test",
        );
        // Members are listed in file order.
        assert_eq!(clusters, vec![vec![3, 4, 5], vec![0, 1, 2]]);
    }
}
