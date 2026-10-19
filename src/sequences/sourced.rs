//! Data structures for describing sequences that originate from some other source (usually a
//! reference genome / transcriptome or another sequence)

use crate::base::Base;
use crate::coords::{Region, Strand};
use crate::sequences::Seq;

/// A sequence that originates from some other source like a reference genome / transcriptome or
/// another sequence.
///
/// Describes the sequence itself ([`Seq`]), and the [`Region`] and [`Strand`] of the original sequence which contained the sequence.
///
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourcedSeq<B: Base> {
    /// The sequence sourced from some other source
    seq: Seq<B>,

    /// Region of the original sequence that contained `seq`
    region: Region,

    /// Strand of the original sequence which contained `seq`
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

    /// Returns the source [`Region`] this sequence came from.
    /// Always given in the forward orientation, so to actually get the [`seq()`] from
    /// the original sequence, check [`strand()`] to figure out whether the original sequence was
    /// reverse complemented
    pub fn region(&self) -> &Region {
        &self.region
    }

    /// Returns the [`Strand`] this sequence originated from if sequence originates from a double stranded
    /// molecule.
    pub fn strand(&self) -> Option<&Strand> {
        self.strand.as_ref()
    }

    /// Return a new SourcedSeq representing the reverse complement.
    ///
    /// Reverse complements the [`SourcedSeq::seq()`] and flips the strand.
    /// [`SourcedSeq::region()`] is intentionally left unchanged.
    pub fn reverse_complement(&self) -> Self {
        Self {
            seq: self.seq().reverse_complement(),
            region: self.region().to_owned(),
            strand: self.strand().map(|s| s.flip()),
        }
    }

    /// Reverse Complements the SourcedSeq in place
    ///
    /// Reverse complements the [`SourcedSeq::seq()`] and flips the strand.
    /// [`SourcedSeq::region()`] is intentionally left unchanged.
    pub fn reverse_complement_in_place(&mut self) {
        // Reverse complement in place
        self.seq.reverse_complement_in_place();

        // Flip strand
        if let Some(unwrapped_strand) = self.strand {
            self.strand = Some(unwrapped_strand.flip())
        };
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
