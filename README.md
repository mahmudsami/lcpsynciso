# lcpsynciso

Cluster long transcript reads (PacBio HiFi, ONT) by shared LCP-syncmer seeds, then
resolve the isoforms within each cluster and build a consensus sequence for each. One
binary, pure Rust, no external tools.

| Command | What it does |
|---|---|
| `predict` | Reads to clusters to isoforms in one run. The input is read once and kept 2-bit packed in memory. |
| `cluster` | Clustering only. |
| `find-isoforms` | Isoforms for an existing cluster assignment. The reads file does not need to be sorted by cluster. |

## Build

Requires Rust 1.79 or newer.

```sh
cargo build --release      # -> target/release/lcpsynciso
```

The aligner dependency (block-aligner) uses AVX2 on x86_64 and NEON on aarch64.

## Usage

```sh
# Everything in one run
lcpsynciso predict reads.fastq.gz --threads 16 -o out/

# Or in two stages
lcpsynciso cluster reads.fastq.gz --emit-assignments -o cluster_out/
lcpsynciso find-isoforms --clusters cluster_out/assignments.tsv \
    --reads reads.fastq.gz -o isoform_out/
```

Input is FASTA or FASTQ, plain, gzip or bgzf. Reads are compared on the forward strand
only, so they should already be oriented (for example PacBio FLNC reads).

Run `lcpsynciso <command> --help` for every option and its default. Note that the default
`--min-shared` differs: 3 for `predict`, 8 for `cluster`.

### Outputs

`predict` and `find-isoforms` write, in the `-o` directory:

| File | Columns |
|---|---|
| `isoform_assignments.tsv` | read_name, cluster_id, isoform_id |
| `isoform_summary.tsv` | cluster_id, isoform_id, n_reads, length |
| `isoforms.fasta` | consensus sequence per isoform |

`cluster` writes `summary.tsv`, `cluster_size_hist.tsv`, and, with `--emit-assignments`,
`assignments.tsv` (read_name, cluster_id).

## How it works

**Clustering.** Seeds are LCP-syncmer block hashes from several levels. Only seeds found
in a moderate number of reads are used. Reads are clustered greedily in file order: each
read joins the cluster its seeds vote for most, weighted by seed level, or starts a new
cluster.

**Isoforms**, per cluster:

1. **Group by structure.** The longest unassigned read is the backbone, and every other
   read is tested against it. Minimizer anchors show where the two agree; the gaps
   between anchors and the read ends beyond them are then aligned. A read is rejected if
   the sequences diverge, a gap holds one long indel (an exon or retained intron), or
   too little of the read is matched. Rejected reads seed later groups.
2. **Split on a recurrent indel** (`--min-variant-frac`, off by default): isoforms that
   differ by a few bases, such as a shifted splice site, share the same indel position.
3. **Split by read ends:** reads that share a start peak and an end peak form an isoform;
   truncated reads join the isoform that contains them.
4. **Consensus:** windows between minimizers shared by most reads, each filled with the
   most frequent read substring. No multiple alignment is needed.
5. **Merge duplicates:** isoforms whose consensuses show they are the same transcript.

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
