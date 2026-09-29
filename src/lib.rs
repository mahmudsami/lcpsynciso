//! lcpsynciso: cluster long transcript reads by shared LCP-syncmer seeds, then detect the
//! isoforms within each cluster.
//!
//! Subcommands (dispatched in `main.rs`):
//!   * `predict`       ([`predict`]): reads -> clusters -> isoforms, reading the input once.
//!   * `cluster`       ([`cluster`]): clustering only.
//!   * `find-isoforms` ([`isoform`]): isoforms for an existing cluster assignment.
//!
//! Building blocks: [`seeds`] (LCP-syncmer block hashes, duplicated from synpact), [`seqio`]
//! (FASTA/FASTQ reading), [`read_store`] (2-bit packed in-memory reads), and [`seedmap`] (a map
//! that uses seed hashes as their own hash).

pub mod cluster;
pub mod isoform;
pub mod predict;
pub mod read_store;
pub mod seedmap;
pub mod seeds;
pub mod seqio;
