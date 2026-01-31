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

    pub fn is_forward(&self) -> bool {
        matches!(self.orientation, Orientation::Forward)
    }

    pub fn is_reverse(&self) -> bool {
        matches!(self.orientation, Orientation::Reverse)
    }

    // < Getters (Read only)>
    /// Sequence-local anchor position (1-based).
    pub fn anchor(&self) -> Pos {
        self.anchor
    }

    /// Does the anchor occur within the confines of the sequence
    pub fn sequence_contains_anchor(&self) -> bool {
        self.seq().sequence_contains_position(self.anchor())
    }

    /// Orientation of the stored sequence relative to the reference.
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// Borrow the underlying sequence.
    pub fn seq(&self) -> &Seq<B> {
        &self.seq
    }

    /// Reference coordinate (1-based) corresponding to the first base of the window.
    pub fn start(&self) -> Pos {
        self.start
    }

    /// Reference coordinate (1-based, inclusive) of the last base in the window.
    pub fn end(&self) -> Pos {
        if self.seq.is_empty() {
            self.start
        } else {
            self.start.saturating_add(self.seq.len() - 1)
        }
    }

    /// Return the anchor index as a 0-based index into `seq`.
    ///
    /// Returns `None` if the anchor lies outside the stored sequence.
    pub fn anchor_index0(&self) -> Option<usize> {
        let idx0 = self.anchor.get().checked_sub(1)?;
        (idx0 < self.seq.len()).then_some(idx0)
    }

    /// External reference coordinate (1-based) corresponding to the anchor.
    pub fn anchor_ref_pos(&self) -> Option<Pos> {
        let idx0 = self.anchor_index0()?;
        Some(self.start.saturating_add(idx0))
    }

    /// Return the base at the anchor position.
    pub fn anchor_base(&self) -> Option<&B> {
        let idx = self.anchor_index0()?;
        self.seq.get(idx)
    }

    /// Borrow a reference k-mer of length `k` centered on the anchor position.
    ///
    /// The anchor base lies exactly at the center of the returned k-mer.
    /// `k` must be **odd and non-zero**.
    ///
    /// ## Returns
    /// - `Some(&[B])` if a centered k-mer of length `k` can be represented
    ///   by the stored context window
    /// - `None` otherwise
    ///
    /// ## `None` is returned when:
    /// - `k` is zero or even
    /// - the anchor lies outside the stored sequence
    /// - there is insufficient flanking sequence on either side of the anchor
    ///
    /// ## Notes
    /// - No mutation-type checks are performed (e.g. SNV vs indel)
    /// - No strand normalization or reverse-complementing is applied
    /// - The returned slice is borrowed; no allocation or copying occurs
    pub fn kmer_centered_on_anchor(&self, k: usize) -> Option<&[B]> {
        if k == 0 || k % 2 == 0 {
            return None;
        }
        let center = self.anchor_index0()?;
        let half = k / 2;
        let start = center.checked_sub(half)?;
        let end = center + half + 1;
        if end > self.seq.len() {
            return None;
        }
        // subseq_slice returns Result<&[B]>; convert to Option
        self.seq.subseq_slice(start, end).ok()
    }

    /// Construct a new `ContextWindow`.
    ///
    /// This constructor performs **no validation** of coordinate or anchoring
    /// semantics. In particular, it does not check that:
    ///
    /// - `anchor` lies within the stored sequence
    /// - the window is correctly anchored to a mutation position
    ///
    /// Methods that depend on these invariants return `Option` to signal when they
    /// cannot be satisfied.
    ///
    /// Callers that require stronger guarantees should validate these conditions
    /// externally or provide a checked constructor.
    pub fn new(seq: Seq<B>, start: Pos, anchor: Pos, orientation: Orientation) -> Self {
        Self {
            seq,
            start,
            anchor,
            orientation,
        }
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
