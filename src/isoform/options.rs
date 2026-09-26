//! Isoform-resolution settings: [`Cfg`], its defaults, and the command-line flags and help
//! text that `predict` and `find-isoforms` share.

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
    /// Minimum identity Eq/(Eq+X) of an aligned gap, indel columns excluded. 0 disables.
    pub min_gap_ident: f64,
    /// An indel of at least this many bp in an aligned gap is extra sequence, not error.
    /// 0 disables.
    pub max_indel_run: u32,
    /// Unmatched flank (bp) that means a different structure: at one end if the sequences
    /// diverge there, at both ends if the read merely overhangs the backbone.
    pub max_flank: u32,
    /// Treat a homopolymer (polyA) overhang as library artifact: it is not counted as an
    /// unmatched flank or against `min_cov`, and does not extend the read's endpoint.
    pub polya_clamp: bool,

    // Step 2: split on a recurrent indel. (`variant.rs`)
    /// Fraction of a group that must share one indel at one position to split it. 0 disables.
    pub min_variant_frac: f64,

    // Step 3: split by read ends. (`ends.rs`)
    /// Split by where reads start as well as where they end. False splits by 3' ends only,
    /// and a read that reaches further 5' than an isoform still joins it: in cDNA most reads
    /// are 5'-truncated, so start peaks mark truncation, not transcription start sites.
    pub split_starts: bool,
    /// Group read ends by the peaks they form; false uses a fixed `boundary_tol` grid.
    pub end_modes: bool,
    /// Window (bp) used to find end peaks.
    pub peak_width: u32,
    /// How far (bp) a read end may sit from its peak, and the slack when folding a partial
    /// read into an isoform. The grid cell size when `end_modes` is off.
    pub boundary_tol: u32,

    // Step 4: isoforms.
    /// Minimum reads for an isoform, and for an end peak or either side of a variant split.
    pub min_iso: usize,
    /// Build a consensus per isoform; false emits the longest member read instead.
    pub consensus: bool,

    // Step 5: merge duplicates. (`collapse.rs`)
    /// Merge near-identical isoforms.
    pub collapse: bool,
    /// Largest gap, indel or length difference (bp) between two isoforms that still merge.
    pub collapse_gap: u32,
    /// Merge by containment only if the smaller isoform has <= this fraction of the larger's
    /// reads.
    pub collapse_ratio: f64,
    /// Consensus identity at which two isoforms merge regardless of `collapse_ratio`.
    /// 0 disables.
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
            polya_clamp: true,
            min_variant_frac: 0.0,
            split_starts: true,
            end_modes: true,
            peak_width: 10,
            boundary_tol: 150,
            min_iso: 3,
            consensus: true,
            collapse: true,
            collapse_gap: 3,
            collapse_ratio: 0.5,
            collapse_ident: 0.98,
        }
    }
}

/// Apply `flag`, if it is one of the shared isoform flags, taking its value from `next`.
/// Returns false for any other flag. The minimizer k/w and consensus flags are named
/// differently per subcommand, so each parses those itself.
pub(crate) fn parse_flag(cfg: &mut Cfg, flag: &str, mut next: impl FnMut() -> String) -> bool {
    match flag {
        "--min-anchors" => cfg.min_anchors = next().parse().unwrap(),
        "--min-cov" => cfg.min_cov = next().parse().unwrap(),
        "--max-gap" => cfg.max_gap = next().parse().unwrap(),
        "--min-gap-ident" => cfg.min_gap_ident = next().parse().unwrap(),
        "--max-indel-run" => cfg.max_indel_run = next().parse().unwrap(),
        "--max-flank" => cfg.max_flank = next().parse().unwrap(),
        "--polya-clamp" => cfg.polya_clamp = true,
        "--no-polya-clamp" => cfg.polya_clamp = false,
        "--min-variant-frac" => cfg.min_variant_frac = next().parse().unwrap(),
        "--split-starts" => cfg.split_starts = true,
        "--no-split-starts" => cfg.split_starts = false,
        "--end-modes" => cfg.end_modes = true,
        "--no-end-modes" => cfg.end_modes = false,
        "--peak-width" => cfg.peak_width = next().parse().unwrap(),
        "--boundary-tol" => cfg.boundary_tol = next().parse().unwrap(),
        "--min-iso" => cfg.min_iso = next().parse().unwrap(),
        "--no-collapse" => cfg.collapse = false,
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
                            indel columns excluded; 0 = off
    --max-indel-run N       one indel of >= N bp in an aligned gap is    [4]
                            extra sequence, not error; 0 = off
    --max-flank N           unmatched read end (bp) meaning a different  [25]
                            structure: one end if the sequences diverge there,
                            both ends if the read only overhangs the longest read
    --polya-clamp / --no-polya-clamp
                            ignore homopolymer (polyA) overhang in the   [on]
                            flank and coverage tests
    --min-variant-frac F    split a group when >= F of its reads share   [0]
                            the SAME >=2 bp indel at the SAME position; 0 = off

READ ENDS (a structure group splits into isoforms by where its reads start and end):
    --split-starts / --no-split-starts
                            split by read starts as well as ends; with   [on]
                            --no-split-starts only 3' ends split, and a read
                            reaching further 5' than an isoform joins it
    --end-modes / --no-end-modes
                            group ends by the peaks they form, or on a   [peaks]
                            fixed --boundary-tol grid
    --peak-width N          window (bp) for finding end peaks; real ends [10]
                            can sit ~20 bp apart
    --boundary-tol N        max bp between a read end and its peak, and  [150]
                            the slack when folding partial reads into an isoform;
                            grid cell size with --no-end-modes

MERGING (near-identical isoforms, after consensus):
    --no-collapse           keep near-identical isoforms separate
    --collapse-gap N        max gap or length difference (bp) to merge   [3]
    --collapse-ratio F      merge only if the smaller has <= F x the     [0.5]
                            reads of the larger
    --collapse-ident F      consensus identity at which isoforms merge   [0.98]
                            regardless of --collapse-ratio; 0 = off";
