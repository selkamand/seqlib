use crate::{base::Base, sequence::Seq};

/// A reference context window fetched around a mutation.
///
/// `ContextWindow` stores a stretch of reference sequence plus the metadata needed
/// to interpret it relative to a mutation.
///
/// This type is intended to support downstream context-derived classification
/// (e.g. trinucleotide/pentanucleotide context for SBS-style tallies, or longer
/// flanks for indel classification) **without re-fetching** from a FASTA.
///
/// # Coordinates and anchoring
/// - `start` is the reference coordinate corresponding to the **first base** of `seq`.
/// - `anchor` is the reference coordinate that the window was fetched *around* (e.g.
///   the mutation position). This anchor is expected to lie within the window.
///
/// With these fields, callers can map from reference coordinate → sequence index:
/// `index = anchor - start` (after bounds checks).
///
/// # Orientation
/// `orientation` describes the direction the stored sequence is written in.
/// - [`Orientation::Forward`]: `seq` is in the same orientation as the reference
///   sequence used for fetching.
/// - [`Orientation::Reverse`]: `seq` is stored reverse-oriented relative to the
///   reference sequence.
///
/// **Important:** `ContextWindow` does not perform strand normalization for mutational
/// signatures (e.g. pyrimidine-centering for SBS96). Downstream code should derive
/// normalized k-mers from the stored window as required.
///
/// # Invariants (expected by downstream classification)
/// - `start` and `anchor` are in the same coordinate system as the mutation being annotated.
/// - The window should be large enough for intended k-mer extraction (e.g. at least
///   3 bases for trinucleotide context).
/// - The anchor should be within the window bounds.
///
/// # Display
/// The [`std::fmt::Display`] implementation renders a short provenance string useful
/// for logs and debugging. It is not a stable serialization format.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContextWindow<B: Base> {
    /// The reference bases for this window.
    seq: Seq<B>,

    /// Reference coordinate corresponding to `seq[0]`.
    start: u64,

    /// Reference coordinate that this window is centered around (typically the
    /// mutation position). Expected to lie within the window.
    anchor: u64,

    /// Orientation of the stored sequence relative to the reference used to fetch it.
    orientation: Orientation,
}

impl<B: Base> std::fmt::Display for ContextWindow<B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let end = if self.seq.is_empty() {
            self.start
        } else {
            self.start + self.seq.len() as u64 - 1
        };

        write!(
            f,
            "Sequence: {} | from {}-{} (both 1-based) | Orientation: {}",
            self.seq, self.start, end, self.orientation
        )
    }
}

impl<B: Base> ContextWindow<B> {
    /// Return the number of bases in the stored context window.
    pub fn len(&self) -> usize {
        self.seq.len()
    }
    /// Return true if the stored context window contains no bases.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// Orientation of a stored context window.
///
/// This describes how the `ContextWindow.seq` is written relative to the reference
/// sequence used to fetch it.
///
/// Most pipelines should store contexts as [`Forward`], and perform any strand
/// normalization (e.g. reverse-complementing k-mers for pyrimidine-centering)
/// during classification rather than during fetching.
///
/// [`Forward`]: Orientation::Forward
/// [`Reverse`]: Orientation::Reverse
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Orientation {
    /// As provided by the reference sequence.
    Forward,
    /// Reverse-oriented relative to the reference sequence.
    Reverse,
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Orientation::Forward => "Forward",
            Orientation::Reverse => "Reverse",
        };
        write!(f, "{s}")
    }
}
