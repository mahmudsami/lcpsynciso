//! The output files of `predict` and `find-isoforms`:
//!
//!   isoform_assignments.tsv          read_name, cluster_id, isoform_id
//!   isoform_summary.tsv              cluster_id, isoform_id, n_reads, length
//!   isoforms.fasta                   one consensus per isoform
//!   start_variants.fasta             one consensus per start variant (`--start-variants report`)
//!   start_variant_assignments.tsv    read_name, cluster_id, variant_id, isoform_id
//!
//! Isoform ids are `<cluster_id>.<index within the cluster>`, start variant ids
//! `<cluster_id>.v<its isoform's index>_<k>`. A start variant's reads count for its isoform in
//! the first three files, so they describe the isoforms alone; the last two give the variants
//! and their reads.
//! FASTA headers: `reads=` and `len=`, `alt_start=1` at a validated downstream start, and for an
//! isoform `start_variants=` (how many it has), for a variant `variant_of=` and `start_offset=`
//! (its start relative to the isoform's, bp, + = downstream).

use std::fs::{self, File};
use std::io::{BufWriter, Write};

use super::{Isoform, Read};

/// One cluster's records, formatted in memory so they can be built on a worker thread.
pub(crate) struct ClusterRecords {
    assignments: Vec<u8>,
    summary: Vec<u8>,
    fasta: Vec<u8>,
    variant_fasta: Vec<u8>,
    variant_assignments: Vec<u8>,
}

/// Format the records for `isoforms`, the isoforms of cluster `cid` followed by their start
/// variants (as `merge_duplicates` returns them), whose members index `reads`.
pub(crate) fn format_cluster_records(
    cid: u32,
    reads: &[Read],
    isoforms: &[Isoform],
) -> ClusterRecords {
    let mut records = ClusterRecords {
        assignments: Vec::new(),
        summary: Vec::new(),
        fasta: Vec::new(),
        variant_fasta: Vec::new(),
        variant_assignments: Vec::new(),
    };
    // Each isoform's start variants, by index in `isoforms`.
    let mut variants: Vec<Vec<usize>> = vec![Vec::new(); isoforms.len()];
    for (index, isoform) in isoforms.iter().enumerate() {
        if let Some((of, _)) = isoform.variant_of {
            variants[of].push(index);
        }
    }
    let alt = |isoform: &Isoform| if isoform.alt_start { " alt_start=1" } else { "" };
    for (isoform_index, isoform) in isoforms.iter().enumerate() {
        if isoform.variant_of.is_some() {
            continue;
        }
        let isoform_id = format!("{cid}.{isoform_index}");
        let n_reads = isoform.members.len()
            + variants[isoform_index].iter().map(|&v| isoforms[v].members.len()).sum::<usize>();
        writeln!(records.summary, "{cid}\t{isoform_id}\t{n_reads}\t{}", isoform.consensus.len())
            .unwrap();
        write!(records.fasta, ">{isoform_id} reads={n_reads} len={}", isoform.consensus.len())
            .unwrap();
        if !variants[isoform_index].is_empty() {
            write!(records.fasta, " start_variants={}", variants[isoform_index].len()).unwrap();
        }
        writeln!(records.fasta, "{}", alt(isoform)).unwrap();
        records.fasta.extend_from_slice(&isoform.consensus);
        records.fasta.push(b'\n');
        for &read_idx in &isoform.members {
            writeln!(records.assignments, "{}\t{cid}\t{isoform_id}", reads[read_idx].name).unwrap();
        }
        for (k, &v) in variants[isoform_index].iter().enumerate() {
            let variant = &isoforms[v];
            let variant_id = format!("{cid}.v{isoform_index}_{k}");
            let (_, offset) = variant.variant_of.unwrap();
            writeln!(
                records.variant_fasta,
                ">{variant_id} variant_of={isoform_id} start_offset={offset:+} reads={} len={}{}",
                variant.members.len(),
                variant.consensus.len(),
                alt(variant)
            )
            .unwrap();
            records.variant_fasta.extend_from_slice(&variant.consensus);
            records.variant_fasta.push(b'\n');
            for &read_idx in &variant.members {
                let name = &reads[read_idx].name;
                writeln!(records.assignments, "{name}\t{cid}\t{isoform_id}").unwrap();
                writeln!(records.variant_assignments, "{name}\t{cid}\t{variant_id}\t{isoform_id}")
                    .unwrap();
            }
        }
    }
    records
}

/// Isoforms in a cluster's list that are not start variants.
pub(crate) fn count_isoforms(isoforms: &[Isoform]) -> usize {
    isoforms.iter().filter(|isoform| isoform.variant_of.is_none()).count()
}

pub(crate) struct OutputFiles {
    assignments: BufWriter<File>,
    summary: BufWriter<File>,
    fasta: BufWriter<File>,
    variant_fasta: BufWriter<File>,
    variant_assignments: BufWriter<File>,
}

impl OutputFiles {
    /// Create `dir` and the files, writing the TSV headers.
    pub(crate) fn create(dir: &str) -> Self {
        fs::create_dir_all(dir).expect("mkdir out");
        let open = |name: &str| BufWriter::new(File::create(format!("{dir}/{name}")).unwrap());
        let mut files = OutputFiles {
            assignments: open("isoform_assignments.tsv"),
            summary: open("isoform_summary.tsv"),
            fasta: open("isoforms.fasta"),
            variant_fasta: open("start_variants.fasta"),
            variant_assignments: open("start_variant_assignments.tsv"),
        };
        writeln!(files.assignments, "read_name\tcluster_id\tisoform_id").unwrap();
        writeln!(files.summary, "cluster_id\tisoform_id\tn_reads\tlength").unwrap();
        writeln!(files.variant_assignments, "read_name\tcluster_id\tvariant_id\tisoform_id").unwrap();
        files
    }

    pub(crate) fn append(&mut self, records: &ClusterRecords) {
        self.assignments.write_all(&records.assignments).unwrap();
        self.summary.write_all(&records.summary).unwrap();
        self.fasta.write_all(&records.fasta).unwrap();
        self.variant_fasta.write_all(&records.variant_fasta).unwrap();
        self.variant_assignments.write_all(&records.variant_assignments).unwrap();
    }

    pub(crate) fn flush(&mut self) {
        self.assignments.flush().unwrap();
        self.summary.flush().unwrap();
        self.fasta.flush().unwrap();
        self.variant_fasta.flush().unwrap();
        self.variant_assignments.flush().unwrap();
    }
}
