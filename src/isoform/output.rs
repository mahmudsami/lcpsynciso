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

/// One cluster's records, rendered in memory so they can be built on a worker thread.
pub(crate) struct Rendered {
    assignments: Vec<u8>,
    summary: Vec<u8>,
    fasta: Vec<u8>,
}

/// Render the records for `isos`, the isoforms of cluster `cid`, whose members index `reads`.
pub(crate) fn render(cid: u32, reads: &[Read], isos: &[Isoform]) -> Rendered {
    let mut r = Rendered { assignments: Vec::new(), summary: Vec::new(), fasta: Vec::new() };
    for (k, iso) in isos.iter().enumerate() {
        let iso_id = format!("{cid}.{k}");
        writeln!(r.summary, "{cid}\t{iso_id}\t{}\t{}", iso.members.len(), iso.consensus.len())
            .unwrap();
        writeln!(r.fasta, ">{iso_id} reads={} len={}", iso.members.len(), iso.consensus.len())
            .unwrap();
        r.fasta.extend_from_slice(&iso.consensus);
        r.fasta.push(b'\n');
        for &ri in &iso.members {
            writeln!(r.assignments, "{}\t{cid}\t{iso_id}", reads[ri].name).unwrap();
        }
    }
    r
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

    pub(crate) fn append(&mut self, r: &Rendered) {
        self.assignments.write_all(&r.assignments).unwrap();
        self.summary.write_all(&r.summary).unwrap();
        self.fasta.write_all(&r.fasta).unwrap();
    }

    pub(crate) fn flush(&mut self) {
        self.assignments.flush().unwrap();
        self.summary.flush().unwrap();
        self.fasta.flush().unwrap();
    }
}
