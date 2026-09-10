use crate::base::{
    Alphabet, Base, ChemClass, ConcreteBase, DegenerateBase, DnaBase, IupacDnaBase, IupacRnaBase,
    RnaBase,
};
use crate::coords::{BaseInterval, BasePos};
use crate::error::SequenceError;
use crate::render::SeqStyler;
use core::fmt;

pub(crate) type Result<T> = std::result::Result<T, SequenceError>;

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

/// A DNA sequence (`Seq<IupacDnaBase>`).
///
/// `IupacDnaSeq` is a flexible type for working with DNA sequences.
/// It represents a validated sequence of DNA bases, including IUPAC ambiguity codes,
/// backed by a compact in-memory representation.
///
/// Most users of this crate should prefer `IupacDnaSeq` over the generic `Seq<B>` type
/// unless they are defining new alphabets or writing generic sequence algorithms.
///
/// If you know your sequence will not include any IUPAC ambiguity codes, consider using [`DnaSeq`]
/// instead.
///
/// # Examples
///
/// ```rust
/// use seqlib::sequences::IupacDnaSeq;
///
/// let seq = IupacDnaSeq::new("ACGTN").unwrap();
/// println!("{}", seq);
/// ```
///
/// Internally, this is just a type alias:
/// ```text
/// type DnaSeq = Seq<DnaBase>
/// ```
pub type IupacDnaSeq = Seq<IupacDnaBase>;

/// A DNA sequence
///
/// `DnaSeq` is a struct type for working with DNA sequences.
/// It represents a validated sequence of DNA bases (ACTG) and does NOT allow IUPAC ambiguity codes
///
/// Most users of this crate should prefer `DnaSeq` over the generic `Seq<B>` type
/// unless they are defining new alphabets or writing generic sequence algorithms.
///
/// # Examples
///
/// ```rust
/// use seqlib::sequences::DnaSeq;
///
/// let seq = DnaSeq::new("ACGT").unwrap();
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
/// `IupacRnaSeq` is a flexible type for working with RNA sequences.
/// It represents a validated sequence of RNA bases, including IUPAC ambiguity codes,
/// using `U` instead of `T`.
///
/// As with [`IupacDnaSeq`], most users should prefer this alias rather than constructing
/// a generic `Seq<B>` directly.
///
/// # Examples
///
/// ```rust
/// use seqlib::sequences::IupacRnaSeq;
///
/// let seq = IupacRnaSeq::new("ACGUN").unwrap();
/// println!("{}", seq);
/// ```
///
/// Internally, this is just a type alias:
/// ```text
/// type IupacRnaSeq = Seq<IupacRnaBase>
/// ```
pub type IupacRnaSeq = Seq<IupacRnaBase>;

/// A RNA sequence
///
/// `RnaSeq` is a struct type for working with RNA sequences.
/// It represents a validated sequence of RNA bases (ACUG) and does NOT allow IUPAC ambiguity codes
///
/// Most users of this crate should prefer `RnaSeq` over the generic `Seq<B>` type
/// unless they are defining new alphabets or writing generic sequence algorithms.
///
/// # Examples
///
/// ```rust
/// use seqlib::sequences::RnaSeq;
///
/// let seq = RnaSeq::new("ACGU").unwrap();
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
    /// Safe indexing; returns None if out of bounds.
    pub fn get(&self, idx: usize) -> Option<&B> {
        self.as_slice().get(idx)
    }

    /// Safe mutable indexing, if Seq is mutable internally.
    pub fn get_mut(&mut self, idx: usize) -> Option<&mut B> {
        self.as_mut_slice().get_mut(idx)
    }

    pub fn as_mut_slice(&mut self) -> &mut [B] {
        &mut self.seq
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
    /// use seqlib::sequences::DnaSeq;
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

    /// Returns a shared reference to the underlying vector of bases.
    ///
    /// This exposes the concrete backing container (`Vec<B>`).
    pub fn as_vec(&self) -> &Vec<B> {
        &self.seq
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
    /// use seqlib::sequences::DnaSeq;
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

    /// Get the Position representing the end of the sequence
    pub fn max_pos(&self) -> BasePos {
        BasePos::new(self.len()).unwrap_or_default()
    }

    /// Does interval span a range that exists in this sequence. If sequence is empty no intervals are valid
    pub fn is_interval_valid(&self, interval: &BaseInterval) -> bool {
        match self.is_empty() {
            true => false,
            false => *interval.end() <= self.max_pos(),
        }
    }

    /// Does the sequence contain a particular position, or does it fall outside of the sequence
    /// length
    pub fn sequence_contains_position(&self, pos: BasePos) -> bool {
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
            Some(base) => base
                .try_chemical_class()
                .unwrap()
                .eq(&ChemClass::Pyrimidine),
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
             "
        )
    }

    /// Returns a borrowed slice of bases using Rust-style indices.
    ///
    /// This method follows standard Rust slicing semantics:
    /// - `start` is **0-based** and **inclusive**
    /// - `end` is **0-based** and **exclusive**
    /// - the returned slice is a **borrowed view** into the original sequence
    ///
    /// The returned slice:
    /// - contains bases `start..end`
    /// - performs **no allocation** and **no copying**
    /// - is valid only for the lifetime of `&self`
    ///
    /// # Errors
    ///
    /// Returns an error if `start > end` or if `end` is greater than the sequence length.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use seqlib::sequences::DnaSeq;
    ///
    /// let seq = DnaSeq::new("ACGTAC").unwrap();
    /// let view = seq.slice(1, 4).unwrap();
    ///
    /// // bases 1,2,3
    /// assert_eq!(view.len(), 3);
    /// ```
    pub fn slice(&self, start: usize, end: usize) -> Result<&[B]> {
        if start > end || end > self.len() {
            return Err(SequenceError::InvalidSlice {
                start,
                end,
                len: self.len(),
            });
        }

        Ok(&self.seq[start..end])
    }

    /// Returns a borrowed view of the subsequence defined by a [`Interval`].
    ///
    /// This method is the biologist-facing counterpart to [`Seq::slice`].
    /// It interprets `interval` using the coordinate contract of [`Interval`]
    /// (1-based coordinates with **both ends included**) and returns a
    /// **read-only, zero-copy** view into the sequence.
    ///
    /// The returned slice:
    /// - contains exactly the bases covered by the interval
    /// - does **not** allocate or copy
    /// - is tied to the lifetime of the original sequence
    ///
    /// # Errors
    ///
    /// Returns an error if the interval falls outside the bounds of the sequence.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use seqlib::{coords::{Pos, Interval}, sequences::{BaseSliceExt, DnaSeq}};
    ///
    /// let seq = DnaSeq::new("ACGTAC").unwrap();
    /// let interval = Interval::new(Pos::new(2).unwrap(), Pos::new(4).unwrap()).unwrap(); // 2..=4
    ///
    /// let slice = seq.subseq_slice(&interval).unwrap();
    /// assert_eq!(slice.to_string_upper(), "CGT");
    /// ```
    pub fn subseq_slice(&self, interval: &BaseInterval) -> Result<&[B]> {
        // Convert interval (1-based inclusive) to Rust indices (0-based, end-exclusive).
        let start = interval.start().as_0based_index();
        let end_exclusive = interval.end().as_0based_index() + 1;

        self.slice(start, end_exclusive)
    }

    /// Returns a borrowed view of the subsequence covered by an [`Interval`].
    ///
    /// Unlike [`Seq::subseq_slice`] this function never throws an error.
    /// If the interval extends beyond the bounds of the sequence, we just return whatever bases
    /// are covered.
    /// If the both the start and end of Interval are outsied the sequence bounds, we just return
    /// an empty slice.
    ///
    /// Returns a **read-only, zero-copy** view into the sequence plus a revised interval
    /// representing the part of the sequence that was actually returned (end clamped to sequence
    /// bounds)
    ///
    /// The returned slice:
    /// - contains the bases covered by the interval
    /// - does **not** allocate or copy
    /// - is tied to the lifetime of the original sequence
    /// - does **not** guarantee the length matches the interval size
    ///
    ///
    /// # Examples
    ///
    /// ```rust
    /// use seqlib::{coords::{Pos, Interval}, sequences::{BaseSliceExt, DnaSeq}, pos};
    ///
    /// let seq = DnaSeq::new("ACGTAC").unwrap();
    ///
    /// // Define interval that extends beyond the sequence
    /// let interval = Interval::new(Pos::new(2).unwrap(), Pos::new(100).unwrap()).unwrap(); // 2..=4
    ///
    /// // Grab the slice of sequence covered by the range and the corresponding clamped interval
    /// let (slice, clamped_interval) = seq.subseq_covered_slice(&interval);
    /// assert_eq!(slice.to_string_upper(), "CGTAC");
    /// assert_eq!(clamped_interval, Some(Interval::new(pos1!(2), pos1!(6)).unwrap()));
    /// ```
    pub fn subseq_covered_slice(&self, interval: &BaseInterval) -> (&[B], Option<BaseInterval>) {
        // Convert interval (1-based inclusive) to Rust indices (0-based, end-exclusive).
        let start = interval.start().as_0based_index();
        if start > self.len() - 1 {
            return (&self.seq[0..0], None);
        }

        let sequence_contains_end_position = self.sequence_contains_position(*interval.end());
        let end_exclusive = match sequence_contains_end_position {
            true => interval.end().as_0based_index() + 1,
            false => self.len(),
        };

        let new_interval = match sequence_contains_end_position {
            true => interval.clone(),
            false => BaseInterval::try_new(interval.start().to_owned(), self.max_pos().to_owned())
                .expect("Bug in subseq_covered_slice: creation of new end"),
        };

        let slice = self.slice(start, end_exclusive).expect("Bug in subseq_covered_slice: slicing should never fail because interval end should be clamped to seq size in above code");
        (slice, Some(new_interval))
    }

    /// Extracts a subsequence defined by an [`Interval`] as a new, independent [`Seq`].
    ///
    /// This is the classic “subsequence” operation for biologists:
    /// - `interval` uses the coordinate contract of [`Interval`]
    ///   (1-based coordinates with **both ends included**)
    /// - the result is an **owned** `Seq<B>` that does not borrow from the original
    ///
    /// The returned subsequence:
    /// - can be stored, returned, or mutated independently
    /// - does not change the original sequence
    ///
    /// # Errors
    ///
    /// Returns an error if the interval falls outside the bounds of the sequence.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use seqlib::{coords::{Pos, Interval}, sequences::DnaSeq};
    ///
    /// let seq = DnaSeq::new("ACGTAC").unwrap();
    /// let interval = Interval::new(Pos::new(2).unwrap(), Pos::new(4).unwrap()).unwrap(); // 2..=4
    ///
    /// let sub = seq.subseq(&interval).unwrap();
    /// assert_eq!(sub.to_string(), "CGT");
    /// assert_eq!(seq.to_string(), "ACGTAC"); // original unchanged
    /// ```
    pub fn subseq(&self, interval: &BaseInterval) -> Result<Seq<B>> {
        let slice = self.subseq_slice(interval)?;
        Ok(Seq {
            seq: slice.to_vec(),
        })
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
    pub fn format_with_highlight_pos(&self, pos: Option<BasePos>) -> String {
        let idx = pos.map(|position| position.as_0based_index());
        self.format_with_highlight_index(idx)
    }

    /// Highlight a series of bases using a interval. If interval end falls outside of sequence length
    /// it will be annotated with ]>EndPosition
    pub fn format_with_highlight_interval(&self, interval: Option<&BaseInterval>) -> String {
        if let Some(reg) = interval {
            let (start, end) = reg.as_0based_indices();
            let mut s = self.to_string();

            if self.sequence_contains_position(*reg.start()) {
                s.insert(start, '[');
            }

            if self.sequence_contains_position(*reg.end()) {
                s.insert(end + 1, ']');
            } else if !self.is_empty() {
                s.push_str(&format!("{}{}", "]>", reg.end()));
            };
            s
        } else {
            self.to_string()
        }
    }

    /// Format sequence as string, with ANSI color codes to get background highlights
    pub fn format_with_colour(&self) -> String {
        let mut out = String::new();

        for b in self.as_slice().iter() {
            out.push_str(&b.to_colourised_string());
        }

        out
    }

    /// Format the whole sequence using a single style.
    pub fn format_with_style(&self, style: &SeqStyler) -> String {
        style.paint(self.to_string_upper())
    }

    /// Format a sequence using one style inside `interval` and another
    /// style outside it.
    ///
    /// The interval uses 1-based inclusive coordinates. If the interval
    /// extends beyond the end of the sequence, it is clamped to the
    /// sequence length. If it starts beyond the sequence, the entire
    /// sequence receives `outside_style`.
    pub fn format_interval_with_styles(
        &self,
        interval: &BaseInterval,
        outside_style: &SeqStyler,
        inside_style: &SeqStyler,
    ) -> String {
        if self.is_empty() {
            return String::new();
        }

        // Interval start is 1-based inclusive.
        let start = interval.start().as_0based_index();

        // There is no overlap with the sequence.
        if start >= self.len() {
            return outside_style.paint(self);
        }

        // For a 1-based inclusive end coordinate, the numeric position is
        // also the corresponding 0-based exclusive slice boundary.
        //
        // For example:
        // interval 3..=5 -> Rust slice 2..5
        let end_exclusive = interval.end().get().min(self.len());

        let (before, remaining) = self.as_slice().split_at(start);
        let (within, after) = remaining.split_at(end_exclusive - start);

        let mut output = String::new();

        // Avoid emitting empty ANSI sections.
        if !before.is_empty() {
            output.push_str(&outside_style.paint(before.to_string_upper()));
        }

        if !within.is_empty() {
            output.push_str(&inside_style.paint(within.to_string_upper()));
        }

        if !after.is_empty() {
            output.push_str(&outside_style.paint(after.to_string_upper()));
        }

        output
    }

    /// Highlight an interval using the standard sequence styles.
    pub fn format_with_coloured_interval(&self, interval: &BaseInterval) -> String {
        self.format_interval_with_styles(interval, &SeqStyler::DIMMED, &SeqStyler::HIGHLIGHT)
    }

    /// Mutate a sequence changing an interval to a new sequence
    pub fn mutate(&mut self, interval: BaseInterval, new: &Seq<B>) -> Result<Self> {
        let mut cloned_seq = self.clone();
        cloned_seq.mutate_in_place(interval, new)?;

        Ok(cloned_seq)
    }

    pub fn mutate_in_place(&mut self, interval: BaseInterval, new: &Seq<B>) -> Result<()> {
        if !self.is_interval_valid(&interval) {
            return Err(SequenceError::FailedMutateInvalidInterval {
                interval,
                seqlength: self.len(),
            });
        }

        let start = interval.start().as_0based_index();
        let end = interval.end().get();

        // Replace range with expected values
        self.seq.splice(start..end, new.as_slice().iter().cloned());

        Ok(())
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

// Implement functions specific to sequences of concrete (unambiguous) bases
// See [`ConcreteBase`] for details on exactly what this means
impl<B: ConcreteBase> Seq<B> {
    /// Returns `true` if the sequence is palindromic
    ///
    /// A sequence is considered *palindromic* if it is identical to its
    /// **reverse complement**. A *biological palindrome*
    /// is **NOT** the same as the a *mirror palindrome* (like racecar).
    ///
    /// This method returns `true` **only if all of the following hold**:
    ///
    /// 1. **The sequence length is non-zero and even**
    ///    - Empty sequences are not considered palindromic.
    ///    - For DNA/RNA, no concrete base is self-complementary, so
    ///      odd-length sequences cannot form concrete palindromes.
    ///
    /// 2. **Each base matches the complement of its mirrored base**
    ///    - For every position `i` in the first half of the sequence,
    ///      `seq[i] == complement(seq[n - 1 - i])` must hold.
    ///
    /// In genetics, palindromic nucleotide sequences are of interest since they may form secondary structures like hairpins.
    /// They also alow homodimer enzymes to recognise to interact with recognition sequences symmetrically.
    /// This is why they make common  
    ///
    /// See also [`try_is_palindromic`] if you have a sequence of potentially ambiguous/degenerate
    /// bases
    //
    /// ## Examples
    ///
    /// ```rust
    /// use seqlib::sequences::{DnaSeq, RnaSeq};
    ///
    /// // A DNA biological palindrome
    /// assert!(DnaSeq::new("GAATTC").unwrap().is_palindromic());
    ///
    /// // An RNA biological palindrome
    /// assert!(RnaSeq::new("GAAUUC").unwrap().is_palindromic());
    ///
    /// // Returns false for mirror (non-genetic) palindrome
    /// assert!(!DnaSeq::new("ATTA").unwrap().is_palindromic());
    ///
    /// // Odd length → false
    /// assert!(!DnaSeq::new("AAA").unwrap().is_palindromic());
    /// ```
    pub fn is_palindromic(&self) -> bool {
        let n = self.len();

        // Empty sequences are not considered palindromes
        if n == 0 {
            return false;
        };

        // Only even numbered sequences can be palindromes
        if !n.is_multiple_of(2) {
            return false;
        }

        // Check palindrome status
        for i in 0..(n / 2) {
            if self.seq[i] != self.seq[n - 1 - i].complement() {
                return false;
            }
        }
        true
    }
}

// Implement functions specific to sequences of degenerate (potentially ambiguous) bases
// See [`ConcreteBase`] for details on exactly what this means
impl<B: DegenerateBase> Seq<B> {
    /// Check whether a sequence **palindromic**.
    ///
    /// A sequence is considered *palindromic* if it is identical to its
    /// **reverse complement** at the level of *concrete nucleotide bases*.
    /// This method only returns `Ok(true)` when palindromicity can be established
    /// with **100% certainty**.
    ///
    /// ## Certainty guarantees
    ///
    /// This method returns Ok(true) **only if all of the following hold**:
    ///
    /// 1. **The sequence contains no ambiguous bases**
    ///    - Any IUPAC ambiguity code (e.g. `N`, `R`, `Y`, `S`, etc.) makes it
    ///      impossible to determine palindromicity with certainty, because such
    ///      symbols represent multiple possible concrete bases.
    ///    - Sequences containing *any* ambiguous base always return an [`SequenceError::PalindromeError`].
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
    /// ## Examples
    ///
    /// ```rust
    /// use seqlib::sequences::IupacDnaSeq;
    ///
    /// // A concrete, unambiguous palindrome -> Ok(true)
    /// assert!(IupacDnaSeq::new("GAATTC").unwrap().is_palindromic_checked() == Ok(true));
    ///
    /// // Symbolically palindromic but ambiguous throws an error
    /// assert!(IupacDnaSeq::new("NNNNNN").unwrap().is_palindromic_checked().is_err());
    ///
    /// // Odd length -> Ok(false)
    /// assert!(IupacDnaSeq::new("AAA").unwrap().is_palindromic_checked() == Ok(false));
    /// ```
    pub fn is_palindromic_checked(&self) -> Result<bool> {
        // Any ambiguous characters make it impossible to identify palindromes with certainty.
        // So we start by converting our degenerate sequence to a concrete one
        let concrete_seq = match self.try_to_concrete() {
            Ok(c) => c,
            Err(_) => {
                return Err(SequenceError::AmbiguousPalindrome);
            }
        };

        // Then we can use the normal is_palindromic method
        Ok(concrete_seq.is_palindromic())
    }
    /// Attempts to convert a degenerate sequence into its concrete equivalent.
    ///
    /// This method narrows a sequence whose base alphabet may contain Iupac ambiguity
    /// symbols into the corresponding concrete base alphabet. For example:
    ///
    /// - `Seq<IupacDnaBase>` becomes `Seq<DnaBase>`
    /// - `Seq<IupacRnaBase>` becomes `Seq<RnaBase>`
    ///
    /// Conversion succeeds only if every base in the sequence has an unambiguous
    /// concrete representation. Bases such as `A`, `C`, `G`, and `T`/`U` can be
    /// converted, while ambiguity codes such as `N`, `R`, or `Y` cause conversion
    /// to fail.
    ///
    /// The returned sequence preserves the original base order and length.
    ///
    /// # Errors
    ///
    /// Returns [`SequenceError::CannotConvertDegenerateSequence`] if any base cannot be
    /// represented in the concrete alphabet. The error reports the first ambiguous
    /// base encountered, using a 1-based sequence position.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use seqlib::sequences::IupacDnaSeq;
    ///
    /// let seq = IupacDnaSeq::new("ACGT").unwrap();
    /// let concrete = seq.try_to_concrete().unwrap();
    ///
    /// assert_eq!(concrete.to_string(), "ACGT");
    /// ```
    ///
    /// ```rust
    /// use seqlib::sequences::IupacDnaSeq;
    ///
    /// let seq = IupacDnaSeq::new("ACNT").unwrap();
    ///
    /// assert!(seq.try_to_concrete().is_err());
    /// ```
    pub fn try_to_concrete(&self) -> Result<Seq<B::ConcreteEquivalent>> {
        let mut out = Vec::with_capacity(self.len());

        for (idx, base) in self.as_slice().iter().copied().enumerate() {
            let concrete = base.try_to_concrete().ok_or_else(|| {
                SequenceError::CannotConvertDegenerateSequence {
                    position: BasePos::new(idx + 1).unwrap(),
                    base: base.to_char(),
                }
            })?;

            out.push(concrete);
        }

        Ok(Seq { seq: out })
    }
}

pub trait BaseSliceExt<B: Base> {
    fn to_string_upper(&self) -> String;
}

impl<B: Base> BaseSliceExt<B> for [B] {
    fn to_string_upper(&self) -> String {
        self.iter().map(|b| b.to_char()).collect()
    }
}

pub const fn validate_dna_literal(s: &str) {
    let bytes = s.as_bytes();
    let mut i: usize = 0;

    while i < bytes.len() {
        let b: u8 = bytes[i];

        if IupacDnaBase::from_ascii_const(b).is_none() {
            // Render printable ASCII, otherwise show a placeholder.
            let shown: char = if b >= b' ' && b <= b'~' {
                b as char
            } else {
                '�'
            };

            const_panic::concat_panic!(
                "invalid DNA base at position ",
                i + 1,
                " (",
                shown,
                ")",
                ", byte: ",
                b,
                ". Allowed: A,C,G,T + IUPAC ambiguity codes; case-insensitive.",
            );
        }

        i += 1;
    }
}

#[macro_export]
macro_rules! dna {
    ($lit:literal) => {{
        // Force compile-time validation *at the call site*.
        const _: () = {
            $crate::sequences::validate_dna_literal($lit);
        };

        // Now construct using the existing, single runtime constructor.
        // If your const checker matches `try_from_ascii`, this unwrap is safe.
        $crate::sequences::DnaSeq::new($lit).unwrap()
    }};
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Helpers ---

    fn dna(s: &str) -> IupacDnaSeq {
        IupacDnaSeq::new(s).unwrap()
    }

    fn rna(s: &str) -> IupacRnaSeq {
        IupacRnaSeq::new(s).unwrap()
    }

    // --- Construction / strict alphabets ---

    #[test]
    fn new_rejects_invalid_characters_dna() {
        assert!(IupacDnaSeq::new("ACGTX").is_err());
    }

    #[test]
    fn new_rejects_u_in_dna_strict() {
        // Strict DNA: U is not allowed
        assert!(IupacDnaSeq::new("ACGU").is_err());
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
        assert_eq!(*s.middlebase().unwrap(), IupacDnaBase::A);
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
        assert!(dna("GAATTC").is_palindromic_checked().unwrap());
    }

    #[test]
    fn is_palindromic_false_for_non_palindrome() {
        assert!(!dna("AGACT").is_palindromic_checked().unwrap());
    }

    #[test]
    fn is_palindromic_handles_edge_cases() {
        // Empty: empty sequences are not considered palindromic
        assert!(!dna("").is_palindromic_checked().unwrap());

        // Length-1 DNA cannot be palindromic unless the base equals its own complement.
        // For DNA A<->T and C<->G, so no unambiguous base is self-complementary.
        assert!(!dna("A").is_palindromic_checked().unwrap());
        assert!(!dna("C").is_palindromic_checked().unwrap());

        // Ambiguous IUPAC symbols like  S (C/G) can never be classified as palindromic with
        // certainty since each S might be a different base! Should return error.
        assert!(dna("SAAS").is_palindromic_checked().is_err());

        // Odd length sequences are never palindromes since the middle base will always break the
        // palindrome (it cannot be identical when reverse complemented)
        assert!(!dna("AAA").is_palindromic_checked().unwrap());
    }

    // --- Subsequence methods ---

    #[test]
    fn slice_rust_indices_returns_expected_view() {
        let s = dna("ACGTAC");

        // 0-based, end-exclusive: bases 1,2,3 => C,G,T
        let view = s.slice(1, 4).unwrap();
        assert_eq!(view.to_string_upper(), "CGT");

        // empty slice is allowed when start == end
        let empty = s.slice(2, 2).unwrap();
        assert_eq!(empty.len(), 0);

        // full slice
        let full = s.slice(0, s.len()).unwrap();
        assert_eq!(full.to_string_upper(), "ACGTAC");
    }

    #[test]
    fn slice_rust_indices_errors_on_invalid_ranges() {
        let s = dna("ACGT");

        // start > end
        assert!(matches!(
            s.slice(3, 2),
            Err(SequenceError::InvalidSlice {
                start: 3,
                end: 2,
                len: 4
            })
        ));

        // end out of bounds
        assert!(matches!(
            s.slice(0, 10),
            Err(SequenceError::InvalidSlice {
                start: 0,
                end: 10,
                len: 4
            })
        ));
    }

    #[test]
    fn subseq_slice_by_interval_returns_expected_view_inclusive_1based() {
        let s = dna("ACGTAC");

        // Interval is 1-based inclusive: 2..=4 => C,G,T
        let interval =
            BaseInterval::try_new(BasePos::new(2).unwrap(), BasePos::new(4).unwrap()).unwrap();
        let view = s.subseq_slice(&interval).unwrap();
        assert_eq!(view.to_string_upper(), "CGT");

        // single base: 1..=1 => A
        let r1 = BaseInterval::try_new(BasePos::new(1).unwrap(), BasePos::new(1).unwrap()).unwrap();
        let one = s.subseq_slice(&r1).unwrap();
        assert_eq!(one.to_string_upper(), "A");

        // last base: 6..=6 => C
        let rlast =
            BaseInterval::try_new(BasePos::new(6).unwrap(), BasePos::new(6).unwrap()).unwrap();
        let last = s.subseq_slice(&rlast).unwrap();
        assert_eq!(last.to_string_upper(), "C");
    }

    #[test]
    fn subseq_by_interval_returns_expected_owned_seq_and_does_not_modify_original() {
        let s = dna("ACGTAC");
        let interval =
            BaseInterval::try_new(BasePos::new(2).unwrap(), BasePos::new(4).unwrap()).unwrap();

        let sub = s.subseq(&interval).unwrap();
        assert_eq!(sub.to_string_upper(), "CGT");

        // original unchanged
        assert_eq!(s.to_string_upper(), "ACGTAC");
    }

    #[test]
    fn subseq_owned_is_independent_of_original() {
        let s = dna("ACGTAC");
        let interval =
            BaseInterval::try_new(BasePos::new(2).unwrap(), BasePos::new(4).unwrap()).unwrap();

        let mut sub = s.subseq(&interval).unwrap();
        sub.rev_in_place();

        // subseq changed
        assert_eq!(sub.to_string_upper(), "TGC");

        // original unchanged
        assert_eq!(s.to_string_upper(), "ACGTAC");
    }

    #[test]
    fn subseq_and_subseq_slice_agree_on_content() {
        let s = dna("ACGTAC");
        let interval =
            BaseInterval::try_new(BasePos::new(2).unwrap(), BasePos::new(5).unwrap()).unwrap(); // 2..=5 => CGTA

        let view = s.subseq_slice(&interval).unwrap();
        let owned = s.subseq(&interval).unwrap();

        assert_eq!(view.to_string_upper(), owned.to_string_upper());
    }

    #[test]
    fn subseq_slice_errors_when_interval_out_of_bounds() {
        let s = dna("ACGTAC");

        // End beyond sequence length (len=6). Interval 1..=7 should fail.
        let interval =
            BaseInterval::try_new(BasePos::new(1).unwrap(), BasePos::new(7).unwrap()).unwrap();
        assert!(s.subseq_slice(&interval).is_err());
        assert!(s.subseq(&interval).is_err());
    }

    #[test]
    fn subseq_methods_work_for_rna_too() {
        let s = rna("ACGUAC"); // length 6

        // Interval 2..=4 => C,G,U
        let interval =
            BaseInterval::try_new(BasePos::new(2).unwrap(), BasePos::new(4).unwrap()).unwrap();

        let view = s.subseq_slice(&interval).unwrap();
        assert_eq!(view.to_string_upper(), "CGU");

        let owned = s.subseq(&interval).unwrap();
        assert_eq!(owned.to_string_upper(), "CGU");

        // Rust slice 1..4 => C,G,U
        let rust_view = s.slice(1, 4).unwrap();
        assert_eq!(rust_view.to_string_upper(), "CGU");
    }
}
