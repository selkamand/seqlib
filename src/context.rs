use crate::{base::Base, coord::Pos, sequence::Seq};

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
/// - `anchor` is a **sequence-local 1-based position within `seq`** identifying the base (or
///   anchoring point) that the window was fetched around. For SNVs this is typically the mutated
///   reference base; for indels it may represent the left-anchored VCF position.
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
///
/// # Display
/// The [`std::fmt::Display`] implementation renders a short provenance string useful
/// for logs and debugging. It is not a stable serialization format.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContextWindow<B: Base> {
    /// The reference bases for this window.
    seq: Seq<B>,

    /// Reference coordinate corresponding to `seq[0]`.
    start: Pos,

    /// Sequence-local 1-based position within `seq` that the window is anchored around.
    ///
    /// This is **not** an external coordinate. It is a position within this window.
    /// For example, `anchor of Pos(1) refers to `seq[0]`.
    ///
    anchor: Pos,

    /// Orientation of the stored sequence relative to the reference used to fetch it.
    orientation: Orientation,
}

impl<B: Base> std::fmt::Display for ContextWindow<B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Sequence: {} | from {}-{} (both 1-based) in external sequence | Orientation: {}",
            self.seq,
            self.start,
            self.end(),
            self.orientation
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

    /// Coordinate (1-based, inclusive) in external sequence that corresponds to the end position. Defaults to start
    /// position if sequence is empty
    pub fn end(&self) -> Pos {
        match self.is_empty() {
            true => self.start,
            false => self.start.saturating_add(self.len() - 1usize),
        }
    }

    /// Get anchor position as a 0-based usize, or None if anchor position is greater than sequence
    /// length
    pub fn anchor_index0(&self) -> Option<usize> {
        // anchor is Pos (1-based inside window)
        let a = self.anchor.get();
        if a == 0 {
            return None;
        } // should never happen
        let idx0 = a - 1usize;
        (idx0 < self.seq.len()).then_some(idx0)
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
