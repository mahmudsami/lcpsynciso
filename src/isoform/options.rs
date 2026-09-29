//! Isoform-resolution settings: [`Cfg`], its defaults, and the command-line flags and help
//! text that `predict` and `find-isoforms` share.

/// How read starts split a structure group (`--start-split`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StartSplit {
    /// Every start peak splits, like every end peak.
    On,
    /// Only 3' ends split; a read reaching further 5' than an isoform joins it.
    Off,
    /// As `Off`, except that a downstream start peak that looks like a real transcription
    /// start, not 5' truncation, becomes an isoform of its own (see `ends.rs`).
    Auto,
}

/// Settings for [`super::resolve_cluster`], grouped by the step that reads them.
#[derive(Clone)]
pub struct Cfg {
    // Minimizer anchors.
    /// Minimizer k-mer length.
    pub k: usize,
    /// Minimizer window.
    pub w: usize,

    // Step 1: does a read fit the backbone? (`fit.rs`)
    /// Minimum shared anchors, and minimum chain length, needed to compare two reads.
    pub min_anchors: usize,
    /// Minimum fraction of the read that must be matched; unmatched flanks count against it.
    pub min_cov: f64,
    /// A gap between consecutive anchors longer than this (bp) is aligned; shorter gaps are
    /// accepted as sequencing error.
    pub max_gap: u32,
    /// Minimum identity Eq/(Eq+X) of an aligned gap, indel columns excluded.
    pub min_gap_ident: f64,
    /// An indel of at least this many bp in an aligned gap is extra sequence, not error.
    pub max_indel_run: u32,
    /// Unmatched flank (bp) that means a different structure: at one end if the sequences
    /// diverge there, at both ends if the read merely overhangs the backbone.
    pub max_flank: u32,

    // Step 2: split by read ends. (`ends.rs`)
    /// Whether read starts split a group as well as read ends. In cDNA most reads are
    /// 5'-truncated, so most start peaks mark truncation, not transcription start sites.
    pub start_split: StartSplit,
    /// Window (bp) used to find end peaks.
    pub peak_width: u32,
    /// How far (bp) a read end may sit from its peak, and the slack when folding a partial
    /// read into an isoform.
    pub boundary_tol: u32,

    // Step 3: isoforms.
    /// Minimum reads for an isoform, and for an end peak.
    pub min_iso: usize,

    // Step 4: merge duplicates. (`collapse.rs`)
    /// Largest gap, indel or length difference (bp) between two isoforms that still merge.
    pub collapse_gap: u32,
    /// Merge by containment only if the smaller isoform has <= this fraction of the larger's
    /// reads.
    pub collapse_ratio: f64,
    /// Consensus identity at which two isoforms merge regardless of `collapse_ratio`.
    pub collapse_ident: f64,
}

impl Default for Cfg {
    fn default() -> Self {
        Cfg {
            k: 15,
            w: 10,
            min_anchors: 3,
            min_cov: 0.80,
            max_gap: 40,
            min_gap_ident: 0.80,
            max_indel_run: 4,
            max_flank: 25,
            start_split: StartSplit::On,
            peak_width: 10,
            boundary_tol: 150,
            min_iso: 3,
            collapse_gap: 3,
            collapse_ratio: 0.5,
            collapse_ident: 0.98,
        }
    }
}

/// Apply `flag`, if it is one of the shared isoform flags, taking its value from `next`.
/// Returns false for any other flag. The minimizer k/w flags are named
/// differently per subcommand, so each parses those itself.
pub(crate) fn parse_flag(cfg: &mut Cfg, flag: &str, mut next: impl FnMut() -> String) -> bool {
    match flag {
        "--min-anchors" => cfg.min_anchors = next().parse().unwrap(),
        "--min-cov" => cfg.min_cov = next().parse().unwrap(),
        "--max-gap" => cfg.max_gap = next().parse().unwrap(),
        "--min-gap-ident" => cfg.min_gap_ident = next().parse().unwrap(),
        "--max-indel-run" => cfg.max_indel_run = next().parse().unwrap(),
        "--max-flank" => cfg.max_flank = next().parse().unwrap(),
        "--start-split" => {
            cfg.start_split = match next().as_str() {
                "on" => StartSplit::On,
                "off" => StartSplit::Off,
                "auto" => StartSplit::Auto,
                other => {
                    eprintln!("--start-split takes on, off or auto, not {other}");
                    std::process::exit(2);
                }
            }
        }
        "--split-starts" => cfg.start_split = StartSplit::On,
        "--no-split-starts" => cfg.start_split = StartSplit::Off,
        "--peak-width" => cfg.peak_width = next().parse().unwrap(),
        "--boundary-tol" => cfg.boundary_tol = next().parse().unwrap(),
        "--min-iso" => cfg.min_iso = next().parse().unwrap(),
        "--collapse-gap" => cfg.collapse_gap = next().parse().unwrap(),
        "--collapse-ratio" => cfg.collapse_ratio = next().parse().unwrap(),
        "--collapse-ident" => cfg.collapse_ident = next().parse().unwrap(),
        _ => return false,
    }
    true
}

/// Help for the flags handled by [`parse_flag`], printed by both subcommands.
pub(crate) const RESOLVE_HELP: &str = "\
STRUCTURE (each read is tested against the longest unassigned read of its cluster):
    --min-iso N             min reads for a reported isoform             [3]
    --min-anchors N         min collinear anchors to trust an alignment  [3]
    --min-cov F             read fraction its alignment must cover       [0.80]
    --max-gap N             gap between anchors (bp) above which it is   [40]
                            aligned; shorter ones count as sequencing error
    --min-gap-ident F       min identity Eq/(Eq+X) of an aligned gap,    [0.80]
                            indel columns excluded
    --max-indel-run N       one indel of >= N bp in an aligned gap is    [4]
                            extra sequence, not error
    --max-flank N           unmatched read end (bp) meaning a different  [25]
                            structure: one end if the sequences diverge there,
                            both ends if the read only overhangs the longest read

READ ENDS (a structure group splits into isoforms by where its reads start and end):
    --start-split on|off|auto
                            how read starts split a group: on = every    [on]
                            start peak; off = only 3' ends split, and a read
                            reaching further 5' than an isoform joins it;
                            auto = off, plus a downstream start that looks like
                            a real transcription start (capped reads carry an
                            untemplated 5' G; otherwise a sharp peak) becomes
                            its own isoform. --split-starts / --no-split-starts
                            are on / off
    --peak-width N          window (bp) for finding end peaks; real ends [10]
                            can sit ~20 bp apart
    --boundary-tol N        max bp between a read end and its peak, and  [150]
                            the slack when folding partial reads into an isoform

MERGING (near-identical isoforms, after consensus):
    --collapse-gap N        max gap or length difference (bp) to merge   [3]
    --collapse-ratio F      merge only if the smaller has <= F x the     [0.5]
                            reads of the larger
    --collapse-ident F      consensus identity at which isoforms merge   [0.98]
                            regardless of --collapse-ratio";
