//! Greedy single-pass clustering of reads by shared seeds: the `cluster` subcommand, and
//! the building blocks `predict` uses for its phases A and B.
//!
//! Seeds are synpact's LCP block hashes from levels L1..L`levels`. Only mid-frequency seeds
//! are informative: found in at least `min_occ` reads (shared) and at most `max_occ`
//! (not a repeat). Reads are then clustered one at a time, in file order here and longest
//! first in `predict`:
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
//! reads without also splitting reads that have few seeds.
//!
//! The subcommand reads the file twice: pass 1 counts seed frequencies, pass 2 assigns.
//! Seeds are computed in parallel per batch; assignment is sequential because it depends on
//! read order.

use std::collections::BTreeMap;
use std::fs;
use std::io::{BufWriter, Write};

use rayon::prelude::*;

use crate::hashmap::{fmap, FastMap};
use crate::io::SeqReader;
use crate::seeds::read_seed_levels;

// ── Building blocks shared with `predict` ────────────────────────────────────

/// Seed hash -> (number of reads containing it, level).
pub(crate) type SeedCounts = FastMap<(u32, u8)>;

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
    p: SeedParams,
) -> SeedCounts {
    seqs.fold(fmap, |mut m, seq| {
        for (h, lv) in read_seed_levels(seq, p.k, p.s, p.t, p.levels) {
            m.entry(h).or_insert((0, lv)).0 += 1;
        }
        m
    })
    .reduce(fmap, |mut acc, b| {
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
    for (h, (c, lv)) in part {
        total.entry(h).or_insert((0, lv)).0 += c;
    }
}

/// Seeds found in between `min_occ` and `max_occ` reads, as hash -> level.
pub(crate) fn informative_seeds(counts: &SeedCounts, min_occ: u32, max_occ: u32) -> FastMap<u8> {
    let mut informative: FastMap<u8> = fmap();
    for (h, (c, lv)) in counts.iter() {
        if *c >= min_occ && *c <= max_occ {
            informative.insert(*h, *lv);
        }
    }
    informative
}

/// A read's informative seeds with their vote weights. `weights[i]` is the weight of level
/// `i + 1`; levels beyond it weigh their level number.
pub(crate) fn weighted_seeds(
    seq: &[u8],
    p: SeedParams,
    informative: &FastMap<u8>,
    weights: &[u32],
) -> Vec<(u64, u32)> {
    read_seed_levels(seq, p.k, p.s, p.t, p.levels)
        .into_iter()
        .filter_map(|(h, _)| {
            informative
                .get(&h)
                .map(|&lv| (h, weights.get(lv as usize - 1).copied().unwrap_or(lv as u32)))
        })
        .collect()
}

/// Assigns reads to clusters one at a time, in the order given (see the module docs).
pub(crate) struct GreedyClusterer {
    min_shared: u32,
    /// Fraction of the read's total seed weight the votes must also reach.
    min_shared_frac: f64,
    /// Seed hash -> the cluster that claimed it first.
    owner: FastMap<u32>,
    /// Scratch for one read: cluster id -> total vote weight.
    votes: FastMap<u32>,
    /// Reads per cluster.
    sizes: Vec<u32>,
}

impl GreedyClusterer {
    pub(crate) fn new(min_shared: u32, min_shared_frac: f64) -> Self {
        GreedyClusterer {
            min_shared,
            min_shared_frac,
            owner: fmap(),
            votes: fmap(),
            sizes: Vec::new(),
        }
    }

    /// Assign a read, given its weighted informative seeds, and return its cluster id.
    pub(crate) fn assign(&mut self, seeds: &[(u64, u32)]) -> u32 {
        self.votes.clear();
        let total: u32 = seeds.iter().map(|&(_, w)| w).sum();
        let need = self.min_shared.max((self.min_shared_frac * total as f64).ceil() as u32);
        let mut best = u32::MAX;
        let mut best_w = 0u32;
        for (h, w) in seeds {
            if let Some(&cid) = self.owner.get(h) {
                let e = self.votes.entry(cid as u64).or_insert(0);
                *e += *w;
                if *e > best_w {
                    best_w = *e;
                    best = cid;
                }
            }
        }
        let cid = if best_w >= need && best != u32::MAX {
            best
        } else {
            let new = self.sizes.len() as u32;
            self.sizes.push(0);
            new
        };
        self.sizes[cid as usize] += 1;
        for (h, _) in seeds {
            self.owner.entry(*h).or_insert(cid);
        }
        cid
    }

    pub(crate) fn n_clusters(&self) -> usize {
        self.sizes.len()
    }

    pub(crate) fn sizes(&self) -> &[u32] {
        &self.sizes
    }
}

// ── The `cluster` subcommand ─────────────────────────────────────────────────

pub struct Config {
    pub reads: String,
    pub k: usize,
    pub s: usize,
    pub t: usize,
    pub levels: usize,
    pub min_occ: u32,
    pub max_occ: u32,
    /// Total vote weight needed to join a cluster.
    pub min_shared: u32,
    /// ... and the fraction of the read's total seed weight it must also reach.
    pub min_shared_frac: f64,
    /// Vote weight per level (index 0 = L1).
    pub weights: Vec<u32>,
    pub batch: usize,
    pub max_reads: usize,
    pub threads: usize,
    pub out_dir: String,
    pub emit_assignments: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            reads: String::new(),
            k: 15,
            s: 9,
            t: usize::MAX,
            levels: 3,
            min_occ: 2,
            max_occ: 100_000,
            min_shared: 8,
            // Off: in file order with min_shared 8 it only trades merges for splits.
            min_shared_frac: 0.0,
            weights: vec![1, 2, 5],
            batch: 500_000,
            max_reads: usize::MAX,
            threads: 0,
            out_dir: "lcpsynciso_cluster_out".into(),
            emit_assignments: false,
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
            "--batch" => a.batch = next().parse().unwrap(),
            "--max-reads" => a.max_reads = next().parse().unwrap(),
            "--threads" => a.threads = next().parse().unwrap(),
            "--emit-assignments" => a.emit_assignments = true,
            "-o" | "--out" => a.out_dir = next(),
            "-h" | "--help" => {
                usage();
                std::process::exit(0);
            }
            other => {
                if other.starts_with('-') {
                    eprintln!("unknown flag {other}");
                    std::process::exit(2);
                }
                a.reads = other.into();
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
        "lcpsynciso cluster — greedy clustering of reads by shared, level-weighted seeds

USAGE:
    lcpsynciso cluster <reads.fq[.gz]|.fa[.gz]> [options]

OPTIONS:
    --k N            k-mer length                              [15]
    --s N            s-mer length                              [9]
    --levels N       use seeds from levels L1..N               [3]
    --min-occ N      keep seed if it occurs in >= N reads      [2]
    --max-occ N      ... and in <= N reads                     [100000]
    --level-weights a,b,c   per-level vote weights (L1,L2,..)  [1,2,5]
    --min-shared N   join a cluster if weighted shared >= N    [8]
    --min-shared-frac F  ... and >= F of the read's seed weight [0 = off]
    --batch N        reads per batch                           [500000]
    --max-reads N    stop after N reads                        [all]
    --threads N      worker threads (0 = all)                  [0]
    --emit-assignments  also write read_name -> cluster TSV
    -o, --out DIR    output directory                          [lcpsynciso_cluster_out]"
    );
}

pub fn run(a: Config) {
    if a.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(a.threads).build_global().ok();
    }
    let params = SeedParams { k: a.k, s: a.s, t: a.t, levels: a.levels };

    // ── Pass 1: seed frequencies ──
    eprintln!(
        "[cluster] pass 1: counting seed frequencies (k={} s={} L1..{}) ...",
        a.k, a.s, a.levels
    );
    let mut counts: SeedCounts = fmap();
    let mut n_reads: u64 = 0;
    {
        let mut reader = SeqReader::open(&a.reads);
        let mut batch: Vec<Vec<u8>> = Vec::with_capacity(a.batch);
        let mut done = false;
        while !done {
            batch.clear();
            while batch.len() < a.batch {
                match reader.next() {
                    Some((_, seq)) => {
                        batch.push(seq);
                        n_reads += 1;
                        if n_reads as usize >= a.max_reads {
                            done = true;
                            break;
                        }
                    }
                    None => {
                        done = true;
                        break;
                    }
                }
            }
            if batch.is_empty() {
                break;
            }
            let part = count_seeds(batch.par_iter().map(|seq| seq.as_slice()), params);
            merge_seed_counts(&mut counts, part);
            eprintln!("[cluster]   pass1 {} reads ...", n_reads);
        }
    }
    let total_distinct = counts.len();
    let informative = informative_seeds(&counts, a.min_occ, a.max_occ);
    drop(counts);
    let informative_n = informative.len();
    eprintln!(
        "[cluster] distinct seeds={total_distinct}  informative (in [{},{}])={informative_n}",
        a.min_occ, a.max_occ
    );

    // ── Pass 2: greedy assignment ──
    eprintln!(
        "[cluster] pass 2: greedy clustering (weights={:?}, min-shared={}, min-shared-frac={}) ...",
        a.weights, a.min_shared, a.min_shared_frac
    );
    let mut clusterer = GreedyClusterer::new(a.min_shared, a.min_shared_frac);

    fs::create_dir_all(&a.out_dir).expect("mkdir out");
    let mut assignments = if a.emit_assignments {
        let mut w = BufWriter::new(
            fs::File::create(format!("{}/assignments.tsv", a.out_dir))
                .expect("create assignments.tsv"),
        );
        writeln!(w, "read_name\tcluster_id").unwrap();
        Some(w)
    } else {
        None
    };

    let mut reader = SeqReader::open(&a.reads);
    let mut batch: Vec<(String, Vec<u8>)> = Vec::with_capacity(a.batch);
    let mut seen: u64 = 0;
    let mut done = false;
    while !done {
        batch.clear();
        while batch.len() < a.batch {
            match reader.next() {
                Some(record) => {
                    batch.push(record);
                    seen += 1;
                    if seen as usize >= a.max_reads {
                        done = true;
                        break;
                    }
                }
                None => {
                    done = true;
                    break;
                }
            }
        }
        if batch.is_empty() {
            break;
        }

        let seeds: Vec<Vec<(u64, u32)>> = batch
            .par_iter()
            .map(|(_, seq)| weighted_seeds(seq, params, &informative, &a.weights))
            .collect();
        for (i, read_seeds) in seeds.iter().enumerate() {
            let cid = clusterer.assign(read_seeds);
            if let Some(w) = assignments.as_mut() {
                writeln!(w, "{}\t{}", batch[i].0, cid).unwrap();
            }
        }
        eprintln!("[cluster]   pass2 {} reads, {} clusters ...", seen, clusterer.n_clusters());
    }
    if let Some(mut w) = assignments {
        w.flush().unwrap();
    }

    write_outputs(&a, n_reads, clusterer.sizes(), total_distinct, informative_n);
}

/// Write summary.tsv and cluster_size_hist.tsv, and print a summary to stdout.
fn write_outputs(
    a: &Config,
    n_reads: u64,
    cluster_size: &[u32],
    total_distinct: usize,
    informative_n: usize,
) {
    fs::create_dir_all(&a.out_dir).expect("mkdir out");
    let n_clusters = cluster_size.len();
    let mut singletons = 0u64;
    let mut max_sz = 0u32;
    let mut sum: u64 = 0;
    let mut hist: BTreeMap<u32, u64> = BTreeMap::new();
    for &sz in cluster_size {
        if sz == 1 {
            singletons += 1;
        }
        if sz > max_sz {
            max_sz = sz;
        }
        sum += sz as u64;
        *hist.entry(sz).or_insert(0) += 1;
    }
    let mut body = String::from("cluster_size\tnum_clusters\n");
    for (sz, n) in &hist {
        body.push_str(&format!("{sz}\t{n}\n"));
    }
    fs::write(format!("{}/cluster_size_hist.tsv", a.out_dir), body).unwrap();

    let mut sizes: Vec<u32> = cluster_size.to_vec();
    sizes.sort_unstable_by(|x, y| y.cmp(x));
    let top: Vec<u32> = sizes.iter().take(20).copied().collect();
    let reads_in_clustered = sum - singletons;
    let pct_clustered = 100.0 * reads_in_clustered as f64 / n_reads.max(1) as f64;
    let mean_nonsingle: f64 = {
        let multi: Vec<u32> = cluster_size.iter().copied().filter(|&s| s > 1).collect();
        if multi.is_empty() {
            0.0
        } else {
            multi.iter().map(|&x| x as u64).sum::<u64>() as f64 / multi.len() as f64
        }
    };

    let mut summary = String::new();
    summary.push_str(&format!("reads\t{n_reads}\n"));
    summary.push_str(&format!("k\t{}\ns\t{}\nlevels\t{}\n", a.k, a.s, a.levels));
    summary.push_str(&format!(
        "min_occ\t{}\nmax_occ\t{}\nmin_shared\t{}\nmin_shared_frac\t{}\n",
        a.min_occ, a.max_occ, a.min_shared, a.min_shared_frac
    ));
    summary.push_str(&format!(
        "level_weights\t{}\n",
        a.weights.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
    ));
    summary.push_str(&format!("distinct_seeds\t{total_distinct}\n"));
    summary.push_str(&format!("informative_seeds\t{informative_n}\n"));
    summary.push_str(&format!("clusters\t{n_clusters}\n"));
    summary.push_str(&format!("singleton_clusters\t{singletons}\n"));
    summary.push_str(&format!("multi_read_clusters\t{}\n", n_clusters as u64 - singletons));
    summary.push_str(&format!("largest_cluster\t{max_sz}\n"));
    summary.push_str(&format!("mean_size_multi_clusters\t{mean_nonsingle:.2}\n"));
    summary.push_str(&format!("pct_reads_in_multi_clusters\t{pct_clustered:.2}\n"));
    summary.push_str(&format!(
        "top20_sizes\t{}\n",
        top.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
    ));
    fs::write(format!("{}/summary.tsv", a.out_dir), &summary).unwrap();

    println!();
    println!(
        "  reads={n_reads}  k={} s={} L1..{}  band=[{},{}]  weights={:?}  min_shared={}  min_shared_frac={}",
        a.k, a.s, a.levels, a.min_occ, a.max_occ, a.weights, a.min_shared, a.min_shared_frac
    );
    println!("  distinct seeds={total_distinct}  informative={informative_n}");
    println!("  clusters={n_clusters}  (singletons={singletons}, multi={})", n_clusters as u64 - singletons);
    println!("  reads in multi-read clusters: {reads_in_clustered} ({pct_clustered:.1}%)");
    println!("  largest cluster={max_sz}  mean(multi)={mean_nonsingle:.1}");
    println!(
        "  top sizes: {}",
        top.iter().take(12).map(|x| x.to_string()).collect::<Vec<_>>().join(", ")
    );
    println!();
    println!("  wrote summary.tsv + cluster_size_hist.tsv to {}/", a.out_dir);
}

#[cfg(test)]
mod tests {
    use super::GreedyClusterer;

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
}
