//! Data structures for describing sequences that originate from some other source (usually a
//! reference genome / transcriptome or another sequence)

use crate::base::Base;
use crate::coords::{Region, Strand};
use crate::error::SequenceError;
use crate::sequences::Seq;

/// A sequence that originates from some other source like a reference genome / transcriptome or
/// another sequence.
///
/// The region uses reference-forward coordinates. The bases are written 5′→3′
/// in the orientation indicated by `strand`; `None` means an unstranded molecule.
///
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourcedSeq<B: Base> {
    /// The sequence sourced from some other source
    seq: Seq<B>,

    /// Region of the original sequence that contained `seq`.
    /// Always reference-forward irrespective of strand.
    region: Region,

    /// Orientation of the stored bases, or `None` for an unstranded molecule.
    strand: Option<Strand>,
}

impl<B: Base> std::fmt::Display for SourcedSeq<B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //TODO: sequences can be quite long so we actually want a pretty print (where for long
        //sequences we just print the first 20 or so bases then add elipses). I think we should add
        //this as method to Seq
        write!(
            f,
            "{} | from: {} ({})",
            self.seq,
            self.region,
            self.strand
                .map(|x| x.to_string())
                .unwrap_or("no strand".to_string())
        )
    }
}

impl<B: Base> SourcedSeq<B> {
    /// Create a new [`SourcedSeq`] object representing a sequence originating from some other
    /// source (e.g. a reference genome / transcriptome / other sequence).
    ///
    pub fn new<S, R>(seq: S, region: R, strand: Option<Strand>) -> Self
    where
        S: Into<Seq<B>>,
        R: Into<Region>,
    {
        Self {
            seq: seq.into(),
            region: region.into(),
            strand,
        }
    }

    /// Returns the sourced [`Seq`].
    pub fn seq(&self) -> &Seq<B> {
        &self.seq
    }

    /// Returns the source [`Region`] in reference-forward coordinates.
    /// For a negative strand sequence, its stored bases run in the opposite
    /// direction from the coordinates in this region.
    pub fn region(&self) -> &Region {
        &self.region
    }

    /// Returns the orientation of the stored bases on a double-stranded molecule.
    /// Returns `None` for an unstranded molecule.
    pub fn strand(&self) -> Option<&Strand> {
        self.strand.as_ref()
    }

    /// Return a new SourcedSeq representing the reverse complement.
    ///
    /// Reverse complements the [`SourcedSeq::seq()`] and flips the strand.
    /// [`SourcedSeq::region()`] is intentionally left unchanged since its always reference-forward.
    ///
    /// Returns an error for an unstranded source, which has no opposite strand
    /// to describe with the same reference coordinates.
    pub fn reverse_complement(&self) -> Result<Self, SequenceError> {
        let strand = self
            .strand
            .ok_or(SequenceError::CannotReverseComplementUnstrandedSourcedSeq)?;

        Ok(Self {
            seq: self.seq().reverse_complement(),
            region: self.region().to_owned(),
            strand: Some(strand.flip()),
        })
    }

    /// Reverse Complements the SourcedSeq in place
    ///
    /// Reverse complements the [`SourcedSeq::seq()`] and flips the strand.
    /// [`SourcedSeq::region()`] is intentionally left unchanged. Returns an error
    /// without changing `self` when the source is unstranded.
    pub fn reverse_complement_in_place(&mut self) -> Result<(), SequenceError> {
        let strand = self
            .strand
            .ok_or(SequenceError::CannotReverseComplementUnstrandedSourcedSeq)?;

        // Reverse complement in place
        self.seq.reverse_complement_in_place();

        // Flip strand
        self.strand = Some(strand.flip());

        Ok(())
    }

    /// Format the sourced sequence with ansi colouring of the sequence
    pub fn format_with_colour(&self) -> String {
        format!(
            "{} | from: {} ({})",
            self.seq.format_with_colour(),
            self.region,
            self.strand
                .map(|x| x.to_string())
                .unwrap_or("no strand".to_string())
        )
    }
}
