//! The three output files of `predict` and `find-isoforms`:
//!
//!   isoform_assignments.tsv   read_name, cluster_id, isoform_id
//!   isoform_summary.tsv       cluster_id, isoform_id, n_reads, length
//!   isoforms.fasta            one consensus per isoform
//!
//! Isoform ids are `<cluster_id>.<index within the cluster>`.

use std::fs::{self, File};
use std::io::{BufWriter, Write};

use super::{Isoform, Read};

/// One cluster's records, formatted in memory so they can be built on a worker thread.
pub(crate) struct ClusterRecords {
    assignments: Vec<u8>,
    summary: Vec<u8>,
    fasta: Vec<u8>,
}

/// Format the records for `isoforms`, the isoforms of cluster `cid`, whose members index `reads`.
pub(crate) fn format_cluster_records(
    cid: u32,
    reads: &[Read],
    isoforms: &[Isoform],
) -> ClusterRecords {
    let mut records =
        ClusterRecords { assignments: Vec::new(), summary: Vec::new(), fasta: Vec::new() };
    for (isoform_index, isoform) in isoforms.iter().enumerate() {
        let isoform_id = format!("{cid}.{isoform_index}");
        writeln!(
            records.summary,
            "{cid}\t{isoform_id}\t{}\t{}",
            isoform.members.len(),
            isoform.consensus.len()
        )
        .unwrap();
        writeln!(
            records.fasta,
            ">{isoform_id} reads={} len={}",
            isoform.members.len(),
            isoform.consensus.len()
        )
        .unwrap();
        records.fasta.extend_from_slice(&isoform.consensus);
        records.fasta.push(b'\n');
        for &read_idx in &isoform.members {
            writeln!(records.assignments, "{}\t{cid}\t{isoform_id}", reads[read_idx].name).unwrap();
        }
    }
    records
}

pub(crate) struct OutputFiles {
    assignments: BufWriter<File>,
    summary: BufWriter<File>,
    fasta: BufWriter<File>,
}

impl OutputFiles {
    /// Create `dir` and the three files, writing the TSV headers.
    pub(crate) fn create(dir: &str) -> Self {
        fs::create_dir_all(dir).expect("mkdir out");
        let open = |name: &str| BufWriter::new(File::create(format!("{dir}/{name}")).unwrap());
        let mut files = OutputFiles {
            assignments: open("isoform_assignments.tsv"),
            summary: open("isoform_summary.tsv"),
            fasta: open("isoforms.fasta"),
        };
        writeln!(files.assignments, "read_name\tcluster_id\tisoform_id").unwrap();
        writeln!(files.summary, "cluster_id\tisoform_id\tn_reads\tlength").unwrap();
        files
    }

    pub(crate) fn append(&mut self, records: &ClusterRecords) {
        self.assignments.write_all(&records.assignments).unwrap();
        self.summary.write_all(&records.summary).unwrap();
        self.fasta.write_all(&records.fasta).unwrap();
    }

    pub(crate) fn flush(&mut self) {
        self.assignments.flush().unwrap();
        self.summary.flush().unwrap();
        self.fasta.flush().unwrap();
    }
}
