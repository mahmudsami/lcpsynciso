//! lcpsynciso: cluster long transcript reads by shared LCP-syncmer seeds, then resolve the
//! isoforms within each cluster.
//!
//! Subcommands (dispatched in `main.rs`):
//!   * `predict`       ([`pipeline`]): reads -> clusters -> isoforms, reading the input once.
//!   * `cluster`       ([`cluster`]): clustering only.
//!   * `find-isoforms` ([`isoform`]): isoforms for an existing cluster assignment.
//!
//! Building blocks: [`seeds`] (LCP-syncmer block hashes, duplicated from synpact), [`io`]
//! (FASTA/FASTQ reading), [`seqstore`] (2-bit packed in-memory reads), and [`hashmap`] (a map
//! that uses seed hashes as their own hash).

pub mod cluster;
pub mod hashmap;
pub mod io;
pub mod isoform;
pub mod pipeline;
pub mod seeds;
pub mod seqstore;
