# lcpsynciso

Cluster long transcript reads (PacBio HiFi or ONT) by shared LCP-syncmer seeds, then
detect the isoforms within each cluster and build a consensus sequence for each. One
binary, pure Rust, no external tools.

| Command | What it does |
|---|---|
| `predict` | Clusters the reads and detects isoforms in a single run. The input file is read once and held in memory, 2-bit packed. |
| `cluster` | Clustering only, with the same clustering settings as `predict`. Also reads the input once and holds it in memory, 2-bit packed. |
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
`lcpsynciso <command> --help` for every option and its default. `predict` and `cluster`
cluster identically by default (`--min-shared` 3, `--min-shared-frac` 0.1).

### Outputs

`predict` and `find-isoforms` write, in the `-o` directory:

| File | Contents |
|---|---|
| `isoform_assignments.tsv` | read_name, cluster_id, isoform_id (with a header row) |
| `isoform_summary.tsv` | cluster_id, isoform_id, n_reads, length (with a header row) |
| `isoforms.fasta` | consensus sequence per isoform |
| `start_variants.fasta` | consensus per start variant: an isoform's other start in the same promoter (header: `variant_of=`, `start_offset=`) |
| `start_variant_assignments.tsv` | read_name, cluster_id, variant_id, isoform_id (with a header row) |

A start variant's reads count for its isoform in the first three files, so they describe the
isoforms alone. `--start-variants keep` reports start variants as isoforms instead.

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
in between `--min-occ` and `--max-occ` reads are used; for large inputs the upper cap
rises to `--max-occ-frac` of the reads (default 0.01) when that is more than `--max-occ`,
e.g. 138,395 for 13.8M reads. Reads are clustered greedily,
longest read first (ties in file order), so each cluster is seeded by its most complete
read; `predict` and `cluster` share this step. Each seed already claimed by a cluster
votes for it, weighted by the seed's level. A read joins the winning cluster if its votes
reach `--min-shared` and also `--min-shared-frac` of the read's own total seed weight, and
otherwise starts a new cluster. The fraction stops a read from joining an unrelated gene's
cluster through a few seeds shared by chance, such as a repeat in a UTR; longest-first
order needs it, because the longest reads start clusters first.

**Isoforms**, per cluster:

1. **Group by structure.** The longest unassigned read is the backbone, and every other
   unassigned read is tested against it. Minimizer anchors show where the two agree.
   Gaps between anchors longer than `--max-gap` are aligned, and so are the read ends
   beyond the first and last anchor. A read is rejected if the sequences diverge, a gap
   holds one long indel (such as an exon or retained intron), or too little of the read
   is matched. Rejected reads seed later groups.
2. **Split by read ends.** Reads that share a start peak and an end peak form an isoform
   if there are at least `--min-iso` of them. Other reads join the best-supported isoform
   that contains them, or are dropped. Groups smaller than `--min-iso` are dropped.
   `--start-split` sets how starts count. In cDNA most reads are 5'-truncated, so start
   peaks mostly mark truncation, and splitting by every one (`on`) makes the truncated
   majority an isoform of its own and drops the few full-length reads. With `off` only 3'
   end peaks split, and a read reaching further 5' than an isoform joins it. `auto` (the
   default) is `off` plus: a downstream start peak that looks like a real transcription
   start keeps its reads as an isoform of its own, beside the full-length one. In capped cDNA (template
   switching) that means most of its reads carry an untemplated 5' G; without that signal
   (direct RNA, spike-ins), a high stop fraction: of the molecules reaching the peak (reverse
   transcription runs 3' to 5'), the share that stop there. A truncation pile stops a few
   percent of them; at a real start the RNA itself ends. `predict` decides both from the
   whole library before detecting isoforms: the cap signal is on when at least 20% of the
   reads that can show the G carry it (and capped peaks are sharper than the others), and
   each share a start needs is the cut between the two modes of that share over all start
   peaks (Otsu). `--cap-signal`, `--cap-threshold` and `--stop-threshold` set them.
3. **Consensus.** The longest read is cut into windows at minimizers carried by a majority
   of the reads covering that position (at least 2); each window takes the most frequent
   read substring. Where fewer than 3 of the isoform's own reads cover a position, the other
   reads of its structure group (same exons, other ends) vote too. No multiple alignment is
   needed. The 5' end is the longest read's; in a library with a cap signal, its leading Gs
   are trimmed while more reads start at the G than carry one more G before it, since
   template switching adds one untemplated G at most capped ends (`--cap-g`). The 3' end is
   the longest read's. `--consensus-tail vote` takes instead the most frequent read tail from
   the last window boundary most reads hold: on UHRR HiFi it ends exactly at the annotated end
   more often (67% vs 56% of full-length matches), but recovers about 1% fewer transcripts
   full-length and 4% fewer at 99.9% coverage.
4. **Merge duplicates:** isoforms whose consensuses show they are the same transcript. When
   the merged isoform is the more complete one (it reaches further 5', or on to the 3' end
   past reads that stop early), the kept isoform's consensus is rebuilt from all its reads.
   An isoform at a validated start is never merged into the isoform it lies inside
   (`--alt-starts keep`, the default; `merge` lets it merge like any other). Isoforms that
   differ only at the 3' end merge when the smaller has at most `--collapse-ratio` of the
   larger's reads; with `--end-merge polya` they merge instead only if the longer isoform
   is A-rich just past the earlier end (internal priming), whatever the read counts. Two
   isoforms kept apart that share the structure and 3' end peak, with starts at most
   `--start-variant-window` bp apart (50; one promoter), are one isoform and a start variant
   of it: the one with fewer reads is reported in `start_variants.fasta`.

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
  predict/        predict
  cluster/        cluster, plus the reading, seed counting and greedy clustering predict shares
  isoform/        isoform detection and find-isoforms
    mod.rs          overview of the four steps
    fit.rs          step 1: does a read match the backbone
    ends.rs         step 2: split by read ends
    consensus.rs    step 3: consensus
    collapse.rs     step 4: merge duplicates
    anchors.rs      minimizers and anchor chaining
    align.rs        block-aligner calls
    options.rs      settings, defaults, shared flags and help
    output.rs       output files
    stats.rs        diagnostic counters
    find_isoforms.rs find-isoforms
  seeds/          LCP-syncmer seeds (from synpact)
  seqio/          FASTA/FASTQ reader
  read_store.rs   2-bit packed read store
  seedmap.rs      map keyed by seed hashes
```
