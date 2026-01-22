use crate::base::{Alphabet, Base, ChemClass, DnaBase, RnaBase};
use crate::errors::SeqError;
use core::fmt;

/// A biological sequence with an associated alphabet (DNA/RNA).
///
/// `Seq` represents a validated nucleotide sequence (DNA or RNA).
/// The alphabet determines which characters are considered valid and
/// influences downstream operations (e.g. reverse complementing).
///
/// This struct is not appropriate for larger-than-memory sequences.
/// It also completely ignores softmasks (for now)
#[derive(Debug, Clone)]
pub struct Seq<B: Base> {
    seq: Vec<B>, // A vector of objects with the Base Trait
}

impl<B: Base> fmt::Display for Seq<B> {
    /// Formats the sequence as its string representation.
    ///
    /// This prints only the underlying sequence characters, without
    /// additional metadata.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use core::fmt::Write;
        for b in &self.seq {
            f.write_char(b.to_char())?;
        }
        Ok(())
    }
}

/// A DNA sequence (`Seq<DnaBase>`).
///
/// `DnaSeq` is the primary and recommended type for working with DNA sequences.
/// It represents a validated sequence of DNA bases, including IUPAC ambiguity codes,
/// backed by a compact in-memory representation.
///
/// Most users of this crate should prefer `DnaSeq` over the generic `Seq<B>` type
/// unless they are defining new alphabets or writing generic sequence algorithms.
///
/// # Examples
///
/// ```rust
/// use seqlib::sequences::DnaSeq;
///
/// let seq = DnaSeq::new("ACGTN").unwrap();
/// println!("{}", seq);
/// ```
///
/// Internally, this is just a type alias:
/// ```text
/// type DnaSeq = Seq<DnaBase>
/// ```
pub type DnaSeq = Seq<DnaBase>;

/// An RNA sequence (`Seq<RnaBase>`).
///
/// `RnaSeq` is the primary and recommended type for working with RNA sequences.
/// It represents a validated sequence of RNA bases, including IUPAC ambiguity codes,
/// using `U` instead of `T`.
///
/// As with [`DnaSeq`], most users should prefer this alias rather than constructing
/// a generic `Seq<B>` directly.
///
/// # Examples
///
/// ```rust
/// use seqlib::sequences::RnaSeq;
///
/// let seq = RnaSeq::new("ACGUN").unwrap();
/// println!("{}", seq);
/// ```
///
/// Internally, this is just a type alias:
/// ```text
/// type RnaSeq = Seq<RnaBase>
/// ```
pub type RnaSeq = Seq<RnaBase>;

// Other functions we can run on Seq
impl<B: Base> Seq<B> {
    /// Returns the sequence as a `String` using uppercase IUPAC symbols.
    ///
    /// This is a convenience method for turning an in-memory `Seq<B>` back into a normal
    /// string representation (e.g. for printing, logging, or writing FASTA).
    ///
    /// Note: this always uses the uppercase representation, even if the original input
    /// contained lowercase characters.
    pub fn to_string_upper(&self) -> String {
        self.seq.iter().map(|b| b.to_char()).collect()
    }

    /// Returns `true` if the sequence contains zero bases.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Returns the number of bases in the sequence.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Returns the middle base of the sequence.
    ///
    /// If the sequence length is odd, this returns a reference to the base
    /// at the center of the sequence. If the length is zero or even,
    /// `None` is returned.
    ///
    /// # Examples
    ///
    /// ```text
    /// Length 5: index 2 is returned
    /// Length 4: returns None
    /// ```
    ///
    /// # Returns
    ///
    /// - `Some(&B)` if the sequence length is odd
    /// - `None` if the sequence is empty or even-length/// Returns the middle base of the sequence.
    pub fn middlebase(&self) -> Option<&B> {
        let len = self.len();

        if len == 0 || len.is_multiple_of(2) {
            return None;
        }

        let idx_middle = len / 2;
        self.seq.get(idx_middle)
    }

    /// Returns the complement of the sequence.
    ///
    /// This creates a new `Seq<B>` where each base is replaced with its complement
    /// (e.g. DNA: A↔T, C↔G, including IUPAC ambiguity complements).
    ///
    /// The original sequence is not modified.
    pub fn complement(&self) -> Seq<B> {
        let newseq: Vec<B> = self.seq.iter().map(|c| c.complement()).collect();
        Seq { seq: newseq }
    }

    /// Returns the reverse-complement of the sequence.
    ///
    /// This creates a new `Seq<B>` where the sequence order is reversed and each
    /// base is complemented (DNA: A↔T, C↔G; RNA: A↔U, C↔G), including IUPAC ambiguity codes.
    ///
    /// The original sequence is not modified.
    pub fn reverse_complement(&self) -> Seq<B> {
        let newseq: Vec<B> = self.seq.iter().map(|c| c.complement()).rev().collect();
        Seq { seq: newseq }
    }

    /// Returns the alphabet for this sequence (DNA or RNA).
    ///
    /// This is derived from the base type parameter `B` (e.g. `DnaBase` → `Alphabet::DNA`)
    /// and is therefore always consistent with the underlying sequence representation.
    pub fn alphabet(&self) -> Alphabet {
        B::ALPHABET
    }

    /// Returns `true` if any base in the sequence is ambiguous.
    ///
    /// Ambiguous bases include IUPAC codes such as `N`, `R`, `Y`, etc.
    /// This is often useful for guarding operations that require unambiguous input
    /// (e.g. translation).
    pub fn any_ambiguous(&self) -> bool {
        self.seq.iter().any(|b| b.is_ambiguous())
    }

    /// Returns `true` if the middle base of the sequence is a pyrimidine.
    ///
    /// This is a convenience predicate commonly used for motif / context logic.
    ///
    /// Returns `true` if:
    /// - the sequence has an odd length, **and**
    /// - the middle base can be classified as a pyrimidine (C/T for DNA, C/U for RNA),
    ///   including “unambiguous ambiguity” codes such as `Y`.
    ///
    /// Returns `false` if:
    /// - the sequence has no middle base (empty or even length), or
    /// - the middle base is a purine (A/G), or
    /// - the middle base cannot be classified (e.g. `N`).
    ///
    pub fn pyrimidine_centered(&self) -> bool {
        let middlebase = self.middlebase();

        match middlebase {
            Some(base) => base
                .chemical_class()
                .is_some_and(|class| class == ChemClass::Pyrimidine),
            None => false,
        }
    }

    /// Returns a human-readable multi-line description of the sequence.
    ///
    /// This is intended for debugging, logging, and CLI output. It includes basic
    /// derived properties such as length, alphabet, ambiguity, and middle base.
    ///
    /// Note: this returns an owned `String` (it allocates).
    pub fn describe(&self) -> String {
        let heading = "Sequence summary";
        let len = self.len();
        let alphabet = B::ALPHABET;
        let middlebase = match self.middlebase() {
            Some(base) => base.to_string(),
            None => "No middle base".to_string(),
        };
        let pyrimidine_centered = self.pyrimidine_centered();

        let any_ambiguous = self.any_ambiguous();

        format!(
            "----------------\n\
             {heading}\n\
             ----------------\n\
             Alphabet : {alphabet}\n\
             Sequence : {self}\n\
             Length   : {len}\n\
             Any Ambiguous: {any_ambiguous}\n\
             Middle Base: {middlebase}\n\
             Pyrimidine Centered: {pyrimidine_centered}\n\
             "
        )
    }

    /// Parses and validates a sequence from a string slice.
    ///
    /// This is the main construction “gatekeeper”: it converts each ASCII character
    /// into a base of type `B` using [`Base::try_from_ascii`].
    ///
    /// # Errors
    ///
    /// Returns a `SeqError` if any character is not valid for the alphabet implied by `B`.
    ///
    /// # Notes
    ///
    /// - Parsing is byte-based (`&str` is interpreted as ASCII nucleotide symbols).
    /// - Lowercase letters are accepted if `try_from_ascii` is case-insensitive.
    pub fn new(sequence: &str) -> Result<Self, SeqError> {
        let mut seq = Vec::with_capacity(sequence.len());
        for &byte in sequence.as_bytes() {
            seq.push(B::try_from_ascii(byte)?);
        }

        Ok(Self { seq })
    }
}
