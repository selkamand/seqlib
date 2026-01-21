use crate::base::{self, Alphabet, Base};
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
    alphabet: Alphabet,
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

// Other functions we can run on Seq
impl<B: Base> Seq<B> {
    pub fn to_string_upper(&self) -> String {
        self.seq.iter().map(|b| b.to_char()).collect()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }
    pub fn new(sequence: &str, alphabet: Alphabet) -> Result<Seq<B>, SeqError> {
        let mut seq = Vec::with_capacity(sequence.len());

        for &byte in sequence.as_bytes() {
            // Parse one base at a time (case-insensitive handled inside try_from_ascii)
            let base = B::try_from_ascii(byte)?;
            seq.push(base);
        }

        Ok(Seq { seq, alphabet })
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
    /// - `None` if the sequence is empty or even-length
    pub fn middlebase(&self) -> Option<&B> {
        let len = self.len();

        if len == 0 || len.is_multiple_of(2) {
            return None;
        }

        let idx_middle = len / 2;
        self.seq.get(idx_middle)
    }
}
