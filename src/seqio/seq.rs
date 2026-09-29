//! FASTQ/FASTA reader for plain, gzip and bgzf input (bgzf is gzip-compatible, so
//! `MultiGzDecoder` reads it).

use std::fs::File;
use std::io::{BufRead, BufReader};

use flate2::read::MultiGzDecoder;

/// Reverse-complement of a DNA sequence (non-ACGT bytes pass through unchanged).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            x => x,
        })
        .collect()
}

/// Yields `(read_name, sequence)` for every record in `path`.
///
/// Format (FASTQ vs FASTA) is auto-detected from the first record; plain,
/// `.gz`, and `.bgz` inputs are all supported. FASTA sequences may span
/// multiple lines. Sequences are upper-cased.
pub struct SeqReader {
    inner: Box<dyn BufRead>,
    is_fasta: bool,
    pending: Option<String>, // one-line lookahead (needed for multi-line FASTA)
    buf: String,
}

impl SeqReader {
    pub fn open(path: &str) -> Self {
        let file = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
        let is_gzip = path.ends_with(".gz") || path.ends_with(".bgz");
        let mut inner: Box<dyn BufRead> = if is_gzip {
            Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(1 << 20, file))
        };
        // Peek first non-empty line to decide the format.
        let mut first = String::new();
        loop {
            first.clear();
            if inner.read_line(&mut first).unwrap_or(0) == 0 {
                break;
            }
            if !first.trim().is_empty() {
                break;
            }
        }
        let is_fasta = first.starts_with('>');
        let pending = if first.is_empty() { None } else { Some(first) };
        SeqReader { inner, is_fasta, pending, buf: String::new() }
    }

    /// Read the next line, consuming the lookahead if present.
    fn next_line(&mut self) -> Option<String> {
        if let Some(p) = self.pending.take() {
            return Some(p);
        }
        self.buf.clear();
        if self.inner.read_line(&mut self.buf).unwrap_or(0) == 0 {
            return None;
        }
        Some(std::mem::take(&mut self.buf))
    }
}

impl Iterator for SeqReader {
    type Item = (String, Vec<u8>);
    fn next(&mut self) -> Option<(String, Vec<u8>)> {
        if self.is_fasta {
            // header (skip blanks)
            let header = loop {
                let line = self.next_line()?;
                if line.trim().is_empty() {
                    continue;
                }
                break line;
            };
            if !header.starts_with('>') {
                return None;
            }
            let name = header[1..].trim_end().split_whitespace().next().unwrap_or("").to_string();
            // sequence lines until next '>' or EOF
            let mut seq: Vec<u8> = Vec::new();
            loop {
                let line = match self.next_line() {
                    Some(l) => l,
                    None => break,
                };
                if line.starts_with('>') {
                    self.pending = Some(line);
                    break;
                }
                seq.extend(line.trim_end().as_bytes().iter().map(|b| b.to_ascii_uppercase()));
            }
            Some((name, seq))
        } else {
            // FASTQ: @header / seq / + / qual
            let header = self.next_line()?;
            if !header.starts_with('@') {
                return None;
            }
            let name = header[1..].trim_end().split_whitespace().next().unwrap_or("").to_string();
            let seq_line = self.next_line()?;
            let seq: Vec<u8> =
                seq_line.trim_end().as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect();
            let _plus = self.next_line()?;
            let _qual = self.next_line()?;
            Some((name, seq))
        }
    }
}
