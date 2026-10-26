//! Isoform-detection settings: [`IsoformOptions`], its defaults, and the command-line flags and
//! help text that `predict` and `find-isoforms` share.

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

/// Whether reads in this library carry the untemplated 5' G of a capped end (`--cap-signal`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CapSignal {
    /// Decided from the library: `predict` measures it over all structure groups first;
    /// `find-isoforms`, which streams clusters, decides per structure group.
    Auto,
    On,
    Off,
}

/// How start validation works in one library: whether the library shows the cap signal at
/// all, the share of a start candidate's labelled reads that must carry the untemplated G (with
/// a cap signal), and how sharp a candidate must be (without one).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StartLibrary {
    pub cap_signal: bool,
    pub threshold: f64,
    pub sharpness: f64,
}

/// Whether merging may fold an isoform at a validated start into the one it lies inside
/// (`--alt-starts`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AltStarts {
    /// By containment, like any other isoform, when it has at most `collapse_ratio` of the
    /// other's reads.
    Merge,
    /// Never by containment: it starts at a transcription start of its own.
    Keep,
}

/// How merging decides between isoforms that differ only at the 3' end (`--end-merge`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EndMerge {
    /// By read counts: the shorter-ended one merges when it has at most `collapse_ratio` of
    /// the other's reads.
    Ratio,
    /// By internal priming: the earlier 3' end merges, whatever the read counts, if the
    /// longer isoform is A-rich just past it; otherwise both are real ends and stay apart.
    Polya,
}

/// Settings for [`super::detect_isoforms`], grouped by the step that reads them.
#[derive(Clone)]
pub struct IsoformOptions {
    // Minimizer anchors.
    /// Minimizer k-mer length.
    pub minimizer_k: usize,
    /// Minimizer window.
    pub minimizer_w: usize,

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
    /// Whether the library shows the cap signal (`--cap-signal`).
    pub cap_signal: CapSignal,
    /// Share of a candidate's labelled reads that must be capped (`--cap-threshold`); None
    /// sets it from the library.
    pub cap_threshold: Option<f64>,
    /// How sharp a start candidate must be without a cap signal (`--sharp-threshold`); None
    /// sets it from the library.
    pub sharp_threshold: Option<f64>,
    /// The library's cap signal and thresholds, once known; None decides per structure group.
    pub start_library: Option<StartLibrary>,

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
    /// Whether an isoform at a validated start may merge into one it lies inside.
    pub alt_starts: AltStarts,
    /// How isoforms that differ only at the 3' end are merged.
    pub end_merge: EndMerge,
    /// With `EndMerge::Polya`: bases past the earlier 3' end that are tested for A content.
    pub priming_window: usize,
    /// With `EndMerge::Polya`: A fraction of that window that means internal priming.
    pub priming_a_frac: f64,
}

impl Default for IsoformOptions {
    fn default() -> Self {
        IsoformOptions {
            minimizer_k: 15,
            minimizer_w: 10,
            min_anchors: 3,
            min_cov: 0.80,
            max_gap: 40,
            min_gap_ident: 0.80,
            max_indel_run: 4,
            max_flank: 25,
            start_split: StartSplit::Auto,
            peak_width: 10,
            boundary_tol: 150,
            cap_signal: CapSignal::Auto,
            cap_threshold: None,
            sharp_threshold: None,
            start_library: None,
            min_iso: 3,
            collapse_gap: 3,
            collapse_ratio: 0.5,
            collapse_ident: 0.98,
            alt_starts: AltStarts::Keep,
            end_merge: EndMerge::Ratio,
            priming_window: 20,
            priming_a_frac: 0.6,
        }
    }
}

/// Apply `flag`, if it is one of the shared isoform flags, taking its value from `next`.
/// Returns false for any other flag. The minimizer k/w flags are named
/// differently per subcommand, so each parses those itself.
pub(crate) fn apply_isoform_flag(
    options: &mut IsoformOptions,
    flag: &str,
    mut next: impl FnMut() -> String,
) -> bool {
    match flag {
        "--min-anchors" => options.min_anchors = next().parse().unwrap(),
        "--min-cov" => options.min_cov = next().parse().unwrap(),
        "--max-gap" => options.max_gap = next().parse().unwrap(),
        "--min-gap-ident" => options.min_gap_ident = next().parse().unwrap(),
        "--max-indel-run" => options.max_indel_run = next().parse().unwrap(),
        "--max-flank" => options.max_flank = next().parse().unwrap(),
        "--start-split" => {
            options.start_split = match next().as_str() {
                "on" => StartSplit::On,
                "off" => StartSplit::Off,
                "auto" => StartSplit::Auto,
                other => {
                    eprintln!("--start-split takes on, off or auto, not {other}");
                    std::process::exit(2);
                }
            }
        }
        "--split-starts" => options.start_split = StartSplit::On,
        "--no-split-starts" => options.start_split = StartSplit::Off,
        "--peak-width" => options.peak_width = next().parse().unwrap(),
        "--boundary-tol" => options.boundary_tol = next().parse().unwrap(),
        "--cap-signal" => {
            options.cap_signal = match next().as_str() {
                "auto" => CapSignal::Auto,
                "on" => CapSignal::On,
                "off" => CapSignal::Off,
                other => {
                    eprintln!("--cap-signal takes auto, on or off, not {other}");
                    std::process::exit(2);
                }
            }
        }
        "--cap-threshold" => {
            options.cap_threshold = match next().as_str() {
                "auto" => None,
                value => Some(value.parse().unwrap()),
            }
        }
        "--sharp-threshold" => {
            options.sharp_threshold = match next().as_str() {
                "auto" => None,
                value => Some(value.parse().unwrap()),
            }
        }
        "--min-iso" => options.min_iso = next().parse().unwrap(),
        "--collapse-gap" => options.collapse_gap = next().parse().unwrap(),
        "--collapse-ratio" => options.collapse_ratio = next().parse().unwrap(),
        "--collapse-ident" => options.collapse_ident = next().parse().unwrap(),
        "--alt-starts" => {
            options.alt_starts = match next().as_str() {
                "merge" => AltStarts::Merge,
                "keep" => AltStarts::Keep,
                other => {
                    eprintln!("--alt-starts takes merge or keep, not {other}");
                    std::process::exit(2);
                }
            }
        }
        "--end-merge" => {
            options.end_merge = match next().as_str() {
                "ratio" => EndMerge::Ratio,
                "polya" => EndMerge::Polya,
                other => {
                    eprintln!("--end-merge takes ratio or polya, not {other}");
                    std::process::exit(2);
                }
            }
        }
        "--priming-window" => options.priming_window = next().parse().unwrap(),
        "--priming-a-frac" => options.priming_a_frac = next().parse().unwrap(),
        _ => return false,
    }
    true
}

/// Help for the flags handled by [`apply_isoform_flag`], printed by both subcommands.
pub(crate) const ISOFORM_HELP: &str = "\
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
                            how read starts split a group: on = every    [auto]
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
    --cap-signal auto|on|off
                            whether reads carry the untemplated 5' G of  [auto]
                            a capped end (template switching); auto: on if
                            >= 20% of the reads that can show it do (predict
                            measures this over the whole library)
    --cap-threshold F|auto  share of a start peak's reads that must      [auto]
                            carry the G for a real start; auto: the cut between
                            the two modes of this share over the library's start
                            peaks (Otsu)
    --sharp-threshold F|auto
                            without a cap signal, the share of the reads [auto]
                            starting within 300 bp of a start peak that must
                            start within 25 bp of it; auto: the Otsu cut over the
                            library's start peaks (predict)

MERGING (near-identical isoforms, after consensus):
    --collapse-gap N        max gap or length difference (bp) to merge   [3]
    --collapse-ratio F      merge only if the smaller has <= F x the     [0.5]
                            reads of the larger
    --collapse-ident F      consensus identity at which isoforms merge   [0.98]
                            regardless of --collapse-ratio
    --alt-starts keep|merge whether an isoform at a validated start may  [keep]
                            merge into one it lies inside (keep: never;
                            merge: by containment, like any other isoform)
    --end-merge ratio|polya isoforms differing only at the 3' end: ratio [ratio]
                            = merge by --collapse-ratio; polya = merge the
                            earlier end, whatever the read counts, only if the
                            longer isoform is A-rich just past it (internal
                            priming), else keep both
    --priming-window N      bases past the earlier end tested (polya)    [20]
    --priming-a-frac F      A fraction there meaning internal priming    [0.6]";
