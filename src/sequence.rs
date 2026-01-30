use crate::base::{Alphabet, Base, ChemClass, DnaBase, RnaBase};
use crate::coord::{Pos, Region};
use crate::error::{Error, Result};
use core::fmt;

/// A biological sequence with an associated alphabet (DNA/RNA).
///
/// `Seq` represents a validated nucleotide sequence (DNA or RNA).
/// The alphabet determines which characters are considered valid and
/// influences downstream operations (e.g. reverse complementing).
///
/// This struct is not appropriate for larger-than-memory sequences.
/// It also completely ignores softmasks (for now)
#[derive(Debug, Clone, PartialEq, Eq)]
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
/// use seqlib::sequence::DnaSeq;
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
/// use seqlib::sequence::RnaSeq;
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
    /// Returns the middle base of the sequence.
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

    /// Returns the complement of the sequence.
    ///
    /// Modifies an existing `Seq<B>` where each base is replaced with its complement
    /// (e.g. DNA: A↔T, C↔G, including IUPAC ambiguity complements).
    ///
    ///
    pub fn complement_in_place(&mut self) {
        for base in &mut self.seq {
            *base = base.complement();
        }
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

    /// Reverse Complements the sequence in place
    pub fn reverse_complement_in_place(&mut self) {
        self.complement_in_place();
        self.rev_in_place();
    }

    /// Reverse the sequence (in place)
    pub fn rev(&self) -> Seq<B> {
        let newseq: Vec<B> = self.seq.iter().rev().copied().collect();
        Seq { seq: newseq }
    }

    /// Reverse sequence in place
    pub fn rev_in_place(&mut self) {
        self.seq.reverse();
    }

    /// Returns the alphabet for this sequence (DNA or RNA).
    ///
    /// This is derived from the base type parameter `B` (e.g. `DnaBase` → `Alphabet::DNA`)
    /// and is therefore always consistent with the underlying sequence representation.
    pub fn alphabet(&self) -> Alphabet {
        B::ALPHABET
    }

    /// Returns a new sequence with the bases in reverse order.
    ///
    /// This operation does **not** modify the original sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::sequence::DnaSeq;
    ///
    /// let seq = DnaSeq::new("ACGT").unwrap();
    /// assert_eq!(seq.reverse().to_string(), "TGCA");
    /// ```
    pub fn reverse(&self) -> Seq<B> {
        let newseq: Vec<B> = self.seq.iter().copied().rev().collect();
        Seq { seq: newseq }
    }

    /// Returns `true` if any base in the sequence is ambiguous.
    ///
    /// Ambiguous bases include IUPAC codes such as `N`, `R`, `Y`, etc.
    /// This is often useful for guarding operations that require unambiguous input
    /// (e.g. translation).
    pub fn any_ambiguous(&self) -> bool {
        self.seq.iter().any(|b| b.is_ambiguous())
    }

    /// Returns `true` if all bases in the sequence are unambiguous`
    pub fn all_unambiguous(&self) -> bool {
        !self.any_ambiguous()
    }

    /// Returns `true` if the sequence contains zero bases.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Returns `true` if the sequence is **certainly palindromic**.
    ///
    /// A sequence is considered *palindromic* if it is identical to its
    /// **reverse complement** at the level of *concrete nucleotide bases*.
    /// This method uses a **strict, conservative definition** and only returns
    /// `true` when palindromicity can be established with **100% certainty**.
    ///
    /// ## Certainty guarantees
    ///
    /// This method returns `true` **only if all of the following hold**:
    ///
    /// 1. **The sequence contains no ambiguous bases**
    ///    - Any IUPAC ambiguity code (e.g. `N`, `R`, `Y`, `S`, etc.) makes it
    ///      impossible to determine palindromicity with certainty, because such
    ///      symbols represent multiple possible concrete bases.
    ///    - Sequences containing *any* ambiguous base always return `false`.
    ///
    /// 2. **The sequence length is non-zero and even**
    ///    - Empty sequences are not considered palindromic.
    ///    - For DNA/RNA, no unambiguous base is self-complementary, so
    ///      odd-length sequences cannot form concrete palindromes.
    ///
    /// 3. **Each base matches the complement of its mirrored base**
    ///    - For every position `i` in the first half of the sequence,
    ///      `seq[i] == complement(seq[n - 1 - i])` must hold.
    ///
    /// ## What this method does *not* do
    ///
    /// - It does **not** treat symbolically palindromic sequences as palindromes
    ///   if ambiguity is present (e.g. `NNNNNN` or `SAAS`).
    /// - It does **not** consider “possibly palindromic” sequences; the result
    ///   reflects *certainty*, not potential.
    ///
    /// ## Intended use
    ///
    /// This definition is designed for **counting and filtering palindromic
    /// sequences without overcounting**, making it suitable for statistical
    /// analyses, motif discovery, and other contexts where false positives
    /// caused by ambiguity must be avoided.
    ///
    /// ## Examples
    ///
    /// ```rust
    /// use seqlib::sequence::DnaSeq;
    ///
    /// // A concrete, unambiguous palindrome
    /// assert!(DnaSeq::new("GAATTC").unwrap().is_palindromic());
    ///
    /// // Symbolically palindromic but ambiguous → false
    /// assert!(!DnaSeq::new("NNNNNN").unwrap().is_palindromic());
    ///
    /// // Odd length → false
    /// assert!(!DnaSeq::new("AAA").unwrap().is_palindromic());
    /// ```
    pub fn is_palindromic(&self) -> bool {
        // Any ambiguous characters make it impossible to identify palindromes with certainty.
        // We set all these sequences to `false``
        if self.any_ambiguous() {
            return false;
        }

        let n = self.len();

        // Empty sequences are not considered palindromes
        if n == 0 {
            return false;
        };

        // Only even numbered sequences can be palindromes
        if n % 2 != 0 {
            return false;
        }

        // Actually check palindrome status
        for i in 0..(n / 2) {
            if self.seq[i] != self.seq[n - 1 - i].complement() {
                return false;
            }
        }
        true
    }

    pub fn max_pos(&self) -> Pos {
        Pos::new(self.len()).unwrap_or_default()
    }

    /// Does region span a range that exists in this sequence. If sequence is empty - by definition region can Note be
    /// valid
    pub fn is_region_valid(&self, region: &Region) -> bool {
        match self.is_empty() {
            true => false,
            false => region.end() < self.max_pos(),
        }
    }

    /// Does the sequence contain a particular position, or does it fall outside of the sequence
    /// length
    pub fn sequence_contains_position(&self, pos: Pos) -> bool {
        pos <= self.max_pos()
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
            Some(base) => base.chemical_class().eq(&ChemClass::Pyrimidine),
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
        let palindrome = self.is_palindromic();
        let complement = self.complement();
        let reverse_complement = self.reverse_complement();

        format!(
            "----------------\n\
             {heading}\n\
             ----------------\n\
             Alphabet : {alphabet}\n\
             Sequence : {self}\n\
             Reverse Complement: {reverse_complement}\n\
             Complement: {complement}\n\
             Length   : {len}\n\
             Any Ambiguous: {any_ambiguous}\n\
             Middle Base: {middlebase}\n\
             Pyrimidine Centered: {pyrimidine_centered}\n\
             Palindrome: {palindrome}\n
             "
        )
    }

    /// Extract a subsequence as a new, independent `Seq`.
    ///
    /// This method uses **0-based indexing** and a **half-open interval**: `[start, end)`.
    /// That means:
    /// - `start` is included
    /// - `end` is excluded
    ///
    /// This is the **safe, ergonomic default** for most users: the returned subsequence
    /// is an owned copy, so it can be stored, returned from functions, or modified
    /// without affecting the original sequence.
    ///
    /// # Arguments
    ///
    /// * `start` - Start index (inclusive, 0-based)
    /// * `end` - End index (exclusive, 0-based)
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `start > end`
    /// - `end > self.len()`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use seqlib::sequence::DnaSeq;
    ///
    /// let seq = DnaSeq::new("ACGTAC").unwrap();
    /// let sub = seq.subseq(1, 4).unwrap(); // bases 1,2,3
    /// assert_eq!(sub.to_string(), "CGT");
    /// ```
    pub fn subseq(&self, start: usize, end: usize) -> Result<Seq<B>> {
        if start > end || end > self.len() {
            return Err(Error::InvalidSlice {
                start,
                end,
                len: self.len(),
            });
        }

        Ok(Seq {
            seq: self.seq[start..end].to_vec(),
        })
    }

    pub fn subseq_in_place(&mut self, start: usize, end: usize) -> Result<()> {
        if start > end || end > self.len() {
            return Err(Error::InvalidSlice {
                start,
                end,
                len: self.len(),
            });
        }
        // // Remove everything after `end`
        // self.seq.drain(end..);
        // // Remove everything before `start`
        // self.seq.drain(..start);
        self.seq.copy_within(start..end, 0);
        self.seq.truncate(end - start);
        Ok(())
    }

    /// Extract a subsequence as a **borrowed view** (`&[B]`) with **no copying**.
    ///
    /// This method uses **0-based indexing** and a **half-open interval**: `[start, end)`.
    ///
    /// Compared to [`Seq::subseq`], this version:
    /// - does **not** allocate
    /// - does **not** copy bases
    /// - is ideal for quickly computing statistics on many subsequences
    ///
    /// Because it returns a borrowed slice, the returned value is only valid while
    /// the original `Seq` is still alive. Also, Rust will prevent you from calling
    /// in-place mutation methods (like `reverse_in_place`) while this slice is in use.
    ///
    /// # Arguments
    ///
    /// * `start` - Start index (inclusive, 0-based)
    /// * `end` - End index (exclusive, 0-based)
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `start > end`
    /// - `end > self.len()`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use seqlib::sequence::DnaSeq;
    ///
    /// let seq = DnaSeq::new("ACGTAC").unwrap();
    /// let sub = seq.subseq_slice(1, 4).unwrap();
    /// assert_eq!(sub.len(), 3);
    /// ```
    pub fn subseq_slice(&self, start: usize, end: usize) -> Result<&[B]> {
        if start > end || end > self.len() {
            return Err(Error::InvalidSlice {
                start,
                end,
                len: self.len(),
            });
        }

        Ok(&self.seq[start..end])
    }

    // Conversions to other data types

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

    /// Returns the sequence as a read-only slice of bases.
    ///
    /// This provides **borrowed access** to the underlying contiguous storage
    /// without allocating or copying.
    ///
    /// ## What this allows
    ///
    /// The returned slice can be used to:
    /// - index individual bases
    /// - iterate efficiently over the sequence
    /// - take subslices (e.g. for k-mer extraction)
    /// - use standard slice methods such as `windows`, `chunks`, and `split_at`
    ///
    /// ## What this does *not* allow
    ///
    /// - mutation of the sequence
    /// - resizing or reallocation
    /// - violating any invariants of `Seq`
    ///
    /// ## Lifetime and safety
    ///
    /// The returned slice is valid only for the lifetime of `&self`.
    /// Rust’s borrow checker guarantees that the sequence cannot be mutated
    /// while the slice is in use.
    ///
    /// ## Performance
    ///
    /// This method is **zero-cost**:
    /// - no allocation
    /// - no copying
    /// - compiles down to returning a pointer and a length
    ///
    /// ## Examples
    ///
    /// ```rust
    /// use seqlib::sequence::DnaSeq;
    /// use seqlib::base::Base;
    /// let seq = DnaSeq::new("ACGT").unwrap();
    /// let bases = seq.as_slice();
    ///
    /// assert_eq!(bases.len(), 4);
    /// assert_eq!(bases[0].to_char(), 'A');
    /// ```
    pub fn as_slice(&self) -> &[B] {
        &self.seq
    }

    // <- Formatters ->
    ///  the sequence as a string, optionally highlighting a 0-based position.
    ///
    /// If `position` is `Some(i)` and `i < self.len()`, the base at `i` is wrapped
    /// in square brackets like `[...]`. If `position` is out of range, no base is
    /// highlighted.
    pub fn format_with_highlight_index(&self, position: Option<usize>) -> String {
        let mut out = String::new();

        for (i, b) in self.as_slice().iter().enumerate() {
            if position == Some(i) {
                out.push('[');
                out.push_str(&b.to_string());
                out.push(']');
            } else {
                out.push_str(&b.to_string());
            }
        }

        out
    }

    /// Highlight a base using a 1-based sequence-local [`Pos`].
    ///
    /// If the position is out of bounds, no base is highlighted.
    pub fn format_with_highlight_pos(&self, pos: Option<Pos>) -> String {
        let idx = pos.map(|position| position.as_0based_index());
        self.format_with_highlight_index(idx)
    }

    /// Highlight a series of bases using a region. If region end falls outside of sequence length
    /// it will be annotated with ]>EndPosition
    pub fn format_with_highlight_region(&self, region: Option<Region>) -> String {
        if let Some(reg) = region {
            let (start, end) = reg.as_0based_indices();
            let mut s = self.to_string();

            if self.sequence_contains_position(reg.start()) {
                s.insert(start, '[');
            }

            if self.sequence_contains_position(reg.end()) {
                s.insert(end + 1, ']');
            } else if !self.is_empty() {
                s.push_str(&format!("{}{}", "]>", reg.end()));
            };
            s
        } else {
            self.to_string()
        }
    }
    // <- Constructors ->

    /// Parses and validates a sequence from a string slice.
    ///
    /// This is the main construction “gatekeeper”: it converts each ASCII character
    /// into a base of type `B` using [`Base::try_from_ascii`].
    ///
    /// # Errors
    ///
    /// Returns a `Error` if any character is not valid for the alphabet implied by `B`.
    ///
    /// # Notes
    ///
    /// - Parsing is byte-based (`&str` is interpreted as ASCII nucleotide symbols).
    /// - Lowercase letters are accepted if `try_from_ascii` is case-insensitive.
    pub fn new(sequence: &str) -> Result<Self> {
        let mut seq = Vec::with_capacity(sequence.len());
        for &byte in sequence.as_bytes() {
            seq.push(B::try_from_ascii(byte)?);
        }

        Ok(Self { seq })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::DnaBase;

    // --- Helpers ---

    fn dna(s: &str) -> DnaSeq {
        DnaSeq::new(s).unwrap()
    }

    fn rna(s: &str) -> RnaSeq {
        RnaSeq::new(s).unwrap()
    }

    // --- Construction / strict alphabets ---

    #[test]
    fn new_rejects_invalid_characters_dna() {
        assert!(DnaSeq::new("ACGTX").is_err());
    }

    #[test]
    fn new_rejects_u_in_dna_strict() {
        // Strict DNA: U is not allowed
        assert!(DnaSeq::new("ACGU").is_err());
    }

    #[test]
    fn new_rejects_t_in_rna_strict() {
        // Strict RNA: T is not allowed
        assert!(RnaSeq::new("ACGT").is_err());
    }

    #[test]
    fn new_accepts_lowercase() {
        // try_from_ascii is case-insensitive in your Base impls
        let s = dna("acgtn");
        assert_eq!(s.to_string_upper(), "ACGTN");
    }

    // --- Basic properties ---

    #[test]
    fn len_and_is_empty_work() {
        let s = dna("");
        assert_eq!(s.len(), 0);
        assert!(s.is_empty());

        let s2 = dna("A");
        assert_eq!(s2.len(), 1);
        assert!(!s2.is_empty());
    }

    #[test]
    fn alphabet_is_correct() {
        assert_eq!(dna("AC").alphabet(), Alphabet::DNA);
        assert_eq!(rna("AC").alphabet(), Alphabet::RNA);
    }

    // --- middlebase / pyrimidine_centered ---

    #[test]
    fn middlebase_none_for_empty_or_even() {
        assert!(dna("").middlebase().is_none());
        assert!(dna("AC").middlebase().is_none());
        assert!(dna("ACGT").middlebase().is_none());
    }

    #[test]
    fn middlebase_some_for_odd() {
        let s = dna("AGACT"); // len 5, middle index 2 => A
        assert_eq!(*s.middlebase().unwrap(), DnaBase::A);
    }

    #[test]
    fn pyrimidine_centered_true_only_when_middle_is_pyrimidine() {
        // middle is C (pyrimidine)
        assert!(dna("AACAA").pyrimidine_centered()); // middle = C

        // middle is A (purine)
        assert!(!dna("AAGAA").pyrimidine_centered()); // middle = G? Actually "AAGAA" middle is G (purine)
        assert!(!dna("AAAAA").pyrimidine_centered()); // middle = A

        // even length => false
        assert!(!dna("AACC").pyrimidine_centered());
    }

    // --- Complement / reverse / reverse-complement (copying variants) ---

    #[test]
    fn complement_produces_expected_dna() {
        let s = dna("AGACT");
        assert_eq!(s.complement().to_string_upper(), "TCTGA");
    }

    #[test]
    fn rev_and_reverse_match_and_do_not_modify_original() {
        let s = dna("ACGT");
        assert_eq!(s.rev().to_string_upper(), "TGCA");
        assert_eq!(s.reverse().to_string_upper(), "TGCA");
        assert_eq!(s.to_string_upper(), "ACGT"); // original unchanged
    }

    #[test]
    fn reverse_complement_produces_expected_dna() {
        let s = dna("ACGT");
        assert_eq!(s.reverse_complement().to_string_upper(), "ACGT"); // ACGT is its own revcomp
    }

    // --- In-place complement / reverse / reverse-complement ---

    #[test]
    fn complement_in_place_mutates_sequence() {
        let mut s = dna("ACGT");
        s.complement_in_place();
        assert_eq!(s.to_string_upper(), "TGCA");
    }

    #[test]
    fn rev_in_place_mutates_sequence() {
        let mut s = dna("ACGT");
        s.rev_in_place();
        assert_eq!(s.to_string_upper(), "TGCA");
    }

    #[test]
    fn reverse_complement_in_place_matches_copying_version() {
        let s = dna("AGACT");
        let mut t = s.clone();
        t.reverse_complement_in_place();
        assert_eq!(
            t.to_string_upper(),
            s.reverse_complement().to_string_upper()
        );
    }

    // --- Ambiguity predicates ---

    #[test]
    fn ambiguity_checks_work() {
        assert!(dna("ACGT").all_unambiguous());
        assert!(!dna("ACNT").all_unambiguous());

        assert!(!dna("ACGT").any_ambiguous());
        assert!(dna("ACNT").any_ambiguous());
    }

    // --- Palindrome (revcomp symmetry) ---

    #[test]
    fn is_palindromic_true_for_simple_palindrome() {
        // GAATTC is a classic restriction site palindrome (EcoRI)
        assert!(dna("GAATTC").is_palindromic());
    }

    #[test]
    fn is_palindromic_false_for_non_palindrome() {
        assert!(!dna("AGACT").is_palindromic());
    }

    #[test]
    fn is_palindromic_handles_edge_cases() {
        // Empty: empty sequences are not considered palindromic
        assert!(!dna("").is_palindromic());

        // Length-1 DNA cannot be palindromic unless the base equals its own complement.
        // For DNA A<->T and C<->G, so no unambiguous base is self-complementary.
        assert!(!dna("A").is_palindromic());
        assert!(!dna("C").is_palindromic());

        // Ambiguous IUPAC symbols like  S (C/G) can never be classified as symbolic with
        // certainty since each A might be a different base!
        assert!(!dna("SAAS").is_palindromic());

        // Odd length sequences are never palindromes since the middle base will always break the
        // palindrome (it cannot be identical when reverse complemented)
        assert!(!dna("AAA").is_palindromic());
    }

    // --- Subsequence methods ---

    #[test]
    fn subseq_returns_expected_owned_copy() {
        let s = dna("ACGTAC");
        let sub = s.subseq(1, 4).unwrap();
        assert_eq!(sub.to_string_upper(), "CGT");
        assert_eq!(s.to_string_upper(), "ACGTAC"); // original unchanged
    }

    #[test]
    fn subseq_slice_returns_expected_view() {
        let s = dna("ACGTAC");
        let sub = s.subseq_slice(1, 4).unwrap();
        let as_string: String = sub.iter().map(|b| b.to_char()).collect();
        assert_eq!(as_string, "CGT");
    }

    #[test]
    fn subseq_in_place_mutates_without_allocating_new_seq() {
        let mut s = dna("ACGTAC");
        s.subseq_in_place(1, 4).unwrap();
        assert_eq!(s.to_string_upper(), "CGT");
    }

    #[test]
    fn subseq_errors_on_invalid_ranges() {
        let s = dna("ACGT");
        assert!(s.subseq(3, 2).is_err());
        assert!(s.subseq(0, 10).is_err());
        assert!(s.subseq_slice(3, 2).is_err());
        assert!(s.subseq_slice(0, 10).is_err());

        let mut t = dna("ACGT");
        assert!(t.subseq_in_place(3, 2).is_err());
        assert!(t.subseq_in_place(0, 10).is_err());
    }

    // --- as_slice ---

    #[test]
    fn as_slice_exposes_bases_read_only() {
        let s = dna("ACGT");
        let slice = s.as_slice();
        assert_eq!(slice.len(), 4);
        assert_eq!(slice[0], DnaBase::A);
        assert_eq!(slice[3], DnaBase::T);
    }
}
