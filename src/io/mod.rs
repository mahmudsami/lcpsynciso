//! Sequence I/O shared by every subcommand.

mod seq;

pub use seq::{revcomp, SeqReader};
