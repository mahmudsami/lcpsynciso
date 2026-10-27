//! The output files of `predict` and `find-isoforms`:
//!
//!   isoform_assignments.tsv   read_name, cluster_id, isoform_id
//!   isoform_summary.tsv       cluster_id, isoform_id, n_reads, length
//!   isoforms.fasta            one consensus per isoform
//!
//! and the same three for the isoforms whose 3' end step 5 could not validate (no
//! polyadenylation signal and no isoform to merge into), with `.unvalidated_end` before the
//! extension. Isoform ids are `<cluster_id>.<index within the cluster>`, numbered across both.

use std::fs::{self, File};
use std::io::{BufWriter, Write};

use super::{Isoform, Read};

/// The three files' records, in memory.
#[derive(Default)]
struct Records {
    assignments: Vec<u8>,
    summary: Vec<u8>,
    fasta: Vec<u8>,
}

/// One cluster's records, formatted in memory so they can be built on a worker thread.
pub(crate) struct ClusterRecords {
    main: Records,
    unvalidated: Records,
}

/// Format the records for `isoforms`, the isoforms of cluster `cid`, whose members index `reads`.
pub(crate) fn format_cluster_records(
    cid: u32,
    reads: &[Read],
    isoforms: &[Isoform],
) -> ClusterRecords {
    let mut records = ClusterRecords { main: Records::default(), unvalidated: Records::default() };
    for (isoform_index, isoform) in isoforms.iter().enumerate() {
        let isoform_id = format!("{cid}.{isoform_index}");
        let out = if isoform.end_validated { &mut records.main } else { &mut records.unvalidated };
        writeln!(
            out.summary,
            "{cid}\t{isoform_id}\t{}\t{}",
            isoform.members.len(),
            isoform.consensus.len()
        )
        .unwrap();
        writeln!(
            out.fasta,
            ">{isoform_id} reads={} len={}",
            isoform.members.len(),
            isoform.consensus.len()
        )
        .unwrap();
        out.fasta.extend_from_slice(&isoform.consensus);
        out.fasta.push(b'\n');
        for &read_idx in &isoform.members {
            writeln!(out.assignments, "{}\t{cid}\t{isoform_id}", reads[read_idx].name).unwrap();
        }
    }
    records
}

/// The three files, opened with their TSV headers.
struct Files {
    assignments: BufWriter<File>,
    summary: BufWriter<File>,
    fasta: BufWriter<File>,
}

impl Files {
    fn create(dir: &str, infix: &str) -> Self {
        let open = |name: &str, ext: &str| {
            BufWriter::new(File::create(format!("{dir}/{name}{infix}.{ext}")).unwrap())
        };
        let mut files = Files {
            assignments: open("isoform_assignments", "tsv"),
            summary: open("isoform_summary", "tsv"),
            fasta: open("isoforms", "fasta"),
        };
        writeln!(files.assignments, "read_name\tcluster_id\tisoform_id").unwrap();
        writeln!(files.summary, "cluster_id\tisoform_id\tn_reads\tlength").unwrap();
        files
    }

    fn append(&mut self, records: &Records) {
        self.assignments.write_all(&records.assignments).unwrap();
        self.summary.write_all(&records.summary).unwrap();
        self.fasta.write_all(&records.fasta).unwrap();
    }

    fn flush(&mut self) {
        self.assignments.flush().unwrap();
        self.summary.flush().unwrap();
        self.fasta.flush().unwrap();
    }
}

pub(crate) struct OutputFiles {
    main: Files,
    unvalidated: Files,
}

impl OutputFiles {
    /// Create `dir` and the six files, writing the TSV headers.
    pub(crate) fn create(dir: &str) -> Self {
        fs::create_dir_all(dir).expect("mkdir out");
        OutputFiles { main: Files::create(dir, ""), unvalidated: Files::create(dir, ".unvalidated_end") }
    }

    pub(crate) fn append(&mut self, records: &ClusterRecords) {
        self.main.append(&records.main);
        self.unvalidated.append(&records.unvalidated);
    }

    pub(crate) fn flush(&mut self) {
        self.main.flush();
        self.unvalidated.flush();
    }
}
