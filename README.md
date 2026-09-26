# lcpsynciso

Cluster long transcript reads (PacBio HiFi or ONT) by shared LCP-syncmer seeds, then
resolve the isoforms within each cluster and build a consensus sequence for each. One
binary, pure Rust, no external tools.

| Command | What it does |
|---|---|
| `predict` | Clusters the reads and resolves isoforms in a single run. The input file is read once and held in memory, 2-bit packed. |
| `cluster` | Clustering only. |
| `find-isoforms` | Isoforms for an existing cluster assignment. The reads file does not need to be sorted by cluster. Single-threaded. |

## Build

Requires Rust 1.80 or newer (the minimum for the locked `rayon`) and an x86_64 or aarch64
target; no other architecture is configured.

```sh
cargo build --release      # -> target/release/lcpsynciso
```

The aligner dependency, block-aligner, is built with AVX2 on x86_64 and NEON on aarch64.
It does not check the CPU at run time, so on x86_64 the binary needs an AVX2-capable CPU.

## Usage

```sh
# Everything in one run
lcpsynciso predict reads.fastq.gz --threads 16 -o out/

# Or in two stages
lcpsynciso cluster reads.fastq.gz --emit-assignments -o cluster_out/
lcpsynciso find-isoforms --clusters cluster_out/assignments.tsv \
    --reads reads.fastq.gz -o isoform_out/
```

Input is FASTA or FASTQ, plain or gzip/bgzf; compression is recognised by a `.gz` or
`.bgz` file extension. Reads are compared on the forward strand only, so they should
already be oriented (for example PacBio FLNC reads).

`--threads` (default 0, meaning all cores) applies to `predict` and `cluster`. Run
`lcpsynciso <command> --help` for every option and its default. Note that the clustering
thresholds differ: `--min-shared` is 3 for `predict` and 8 for `cluster`, and
`--min-shared-frac` is 0.1 for `predict` and off for `cluster`.

### Outputs

`predict` and `find-isoforms` write, in the `-o` directory:

| File | Contents |
|---|---|
| `isoform_assignments.tsv` | read_name, cluster_id, isoform_id (with a header row) |
| `isoform_summary.tsv` | cluster_id, isoform_id, n_reads, length (with a header row) |
| `isoforms.fasta` | consensus sequence per isoform; the longest member read with `--no-consensus` |

Not every read is assigned. A read is left out of `isoform_assignments.tsv` if its group
has fewer than `--min-iso` reads, if it fits no isoform's ends, or if its cluster is
outside `--min-size`/`--max-size`. On the HiFi SIRV test data at `--min-iso 5`, 2.3% of
reads were left out.

With `--dump-clusters`, `predict` also writes `gene_clusters.tsv` (read_name, cluster_id,
cluster_size, with a header row): the phase B cluster of every read, including the reads
left out of `isoform_assignments.tsv`.

`cluster` writes `summary.tsv`, `cluster_size_hist.tsv`, and, with `--emit-assignments`,
`assignments.tsv` (read_name, cluster_id).

## How it works

**Clustering.** Seeds are LCP-syncmer block hashes from several levels. Only seeds found
in between `--min-occ` and `--max-occ` reads are used. Reads are clustered greedily: in
`predict`, longest read first (ties in file order), so each cluster is seeded by its most
complete read; in `cluster`, which does not hold the reads, in file order. Each seed
already claimed by a cluster votes for it, weighted by the seed's level. A read joins the
winning cluster if its votes reach `--min-shared` and also `--min-shared-frac` of the
read's own total seed weight, and otherwise starts a new cluster. The fraction stops a read
from joining an unrelated gene's cluster through a few seeds shared by chance, such as a
repeat in a UTR; longest-first order needs it, because the longest reads start clusters
first.

**Isoforms**, per cluster:

1. **Group by structure.** The longest unassigned read is the backbone, and every other
   unassigned read is tested against it. Minimizer anchors show where the two agree.
   Gaps between anchors longer than `--max-gap` are aligned, and so are the read ends
   beyond the first and last anchor. A read is rejected if the sequences diverge, a gap
   holds one long indel (such as an exon or retained intron), or too little of the read
   is matched. Rejected reads seed later groups.
2. **Split on a recurrent indel** (`--min-variant-frac`, off by default): isoforms that
   differ by a few bases, such as a shifted splice site, share the same indel position.
3. **Split by read ends.** Reads that share a start peak and an end peak form an isoform
   if there are at least `--min-iso` of them. Other reads join the best-supported isoform
   that contains them, or are dropped. Groups smaller than `--min-iso` are dropped.
   `--start-split` sets how starts count. In cDNA most reads are 5'-truncated, so start
   peaks mostly mark truncation, and splitting by every one (`on`, the default) makes the
   truncated majority an isoform of its own and drops the few full-length reads. With `off`
   only 3' end peaks split, and a read reaching further 5' than an isoform joins it. `auto`
   is `off` plus: a downstream start peak that looks like a real transcription start keeps
   its reads as an isoform of its own, beside the full-length one. In capped cDNA (template
   switching) that means most of its reads carry an untemplated 5' G; without that signal
   (direct RNA, spike-ins), a sharp peak.
4. **Consensus.** The longest read is cut into windows at minimizers carried by a majority
   of the reads covering that position (at least 2); each window takes the most frequent
   read substring. Where fewer than 3 of the isoform's own reads cover a position, the other
   reads of its structure group (same exons, other ends) vote too. No multiple alignment is
   needed.
5. **Merge duplicates:** isoforms whose consensuses show they are the same transcript. When
   the merged isoform is the more complete one (it reaches further 5', or on to the 3' end
   past reads that stop early), the kept isoform's consensus is rebuilt from all its reads.

## Notes

- Defaults were tuned on PacBio HiFi SIRV reads and have not been calibrated for ONT.
- Output can vary slightly between runs, because Rust's standard `HashMap` is seeded
  randomly per process. Isoform ids are not stable across runs; compare results by
  sequence.
- `src/seeds/` is duplicated from [synpact](https://github.com/mahmudsami/synpact) so that
  seed hashes match exactly. Keep it in sync with synpact rather than editing it here.

## Layout

```
src/
  main.rs         subcommand dispatch
  lib.rs          library root
  pipeline/       predict
  cluster/        cluster, plus the seed counting and greedy clustering predict uses
  isoform/        isoform resolution and find-isoforms
    mod.rs          overview of the five steps
    fit.rs          step 1: does a read match the backbone
    variant.rs      step 2: recurrent-indel split
    ends.rs         step 3: split by read ends
    consensus.rs    step 4: consensus
    collapse.rs     step 5: merge duplicates
    anchors.rs      minimizers and anchor chaining
    ba.rs           block-aligner calls
    options.rs      settings, defaults, shared flags and help
    output.rs       output files
    stats.rs        diagnostic counters
    cli.rs          find-isoforms
  seeds/          LCP-syncmer seeds (from synpact)
  io/             FASTA/FASTQ reader
  seqstore.rs     2-bit packed read store
  hashmap.rs      map keyed by seed hashes
```
