use core::fmt;

use crate::{
    base::{Base, ChemClass, DnaBase, RnaBase},
    coord::Pos,
    sequence::Seq,
};

pub type DnaSmallMutation = SmallMutation<DnaBase>;
pub type RnaSmallMutation = SmallMutation<RnaBase>;

/// A small mutation (SNV/MNV/indel) over a specific nucleotide alphabet `B`.
///
/// - `SmallMutation<DnaBase>` is DNA
/// - `SmallMutation<RnaBase>` is RNA
///
/// ## Coordinate and allele semantics
/// - `position` is **1-based** (VCF-style) and refers to the **start** position.
/// - `reference` and `alternative` are stored **as provided** by the caller.
///   No left/right trimming, normalization, or decomposition is performed.
///
/// ## Multi-allelic sites
/// - `multiallelic` indicates that this mutation originated from a site with multiple ALT
///   alleles (e.g. a VCF record with `ALT=A,C`). It is purely metadata for downstream
///   handling; it does not change coordinate or allele semantics.
///
/// ## Filter / QC status (`pass`)
/// - `pass` indicates whether the source record **passed upstream filtering**.
///   In VCF terms this typically corresponds to the `FILTER` field being `PASS`
///   (or `.` depending on your chosen convention).
/// - `pass` is **metadata only**:
///   - it does *not* imply biological validity,
///   - it does *not* change classification (`SNV`/`INDEL`/etc),
///   - and it should not be silently acted on by this type.
///
/// ## Context
/// - `context` is optional sequence context (e.g. trinucleotide context) if already
///   computed by the caller. This is commonly added later after reference lookup.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmallMutation<B: Base> {
    chromosome: String,
    position: Pos, // 1-based position (start)
    reference: Seq<B>,
    alternative: Seq<B>,
    multiallelic: bool,
    pass: bool,
    context: Option<Seq<B>>,
}

// Implement the `fmt::Display` trait for `Point`.
impl<B: Base> fmt::Display for SmallMutation<B> {
    /// Render a compact, human-readable representation of the mutation.
    ///
    /// Format:
    /// `chrom:pos REF>ALT (delta: D; class: C; multiallelic: M; pass: P)`
    ///
    /// Intended for logging / CLI output rather than stable serialization.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use the write! macro to format the output.
        write!(
            f,
            "{}:{} {}>{} (delta: {}; class: {}; multiallelic:{}; pass:{})",
            self.chromosome,
            self.position,
            self.reference,
            self.alternative,
            self.delta(),
            self.class(),
            self.multiallelic,
            self.pass
        )
    }
}

impl<B: Base> SmallMutation<B> {
    /// Construct a new [`SmallMutation`].
    ///
    /// This constructor does **not** attempt to validate biological correctness beyond
    /// what is already guaranteed by [`Seq<B>`] (i.e. sequences conform to the alphabet).
    ///
    /// In particular, this type:
    /// - assumes `position` is **1-based** (VCF-style coordinates)
    /// - stores `reference` and `alternative` as provided (no trimming/normalization)
    /// - treats `multiallelic` as an external flag (e.g. derived from a multi-ALT record)
    /// - treats `pass` as an external flag describing upstream filtering/QC outcome
    /// - accepts an optional `context` sequence if the caller has already computed it
    ///
    /// If you need allele normalization (left/right trimming of shared prefix/suffix),
    /// do it before constructing this type.
    ///
    /// # Parameters
    /// - `chromosome`: reference sequence / contig name (e.g. `"chr1"`)
    /// - `position`: 1-based start coordinate
    /// - `reference`: reference allele sequence
    /// - `alternative`: alternative allele sequence
    /// - `multiallelic`: whether the originating site had multiple ALT alleles
    /// - `pass`: whether the originating record passed upstream filters/QC
    /// - `context`: optional context sequence (e.g. trinucleotide context)
    pub fn new(
        chromosome: String,
        position: Pos,
        reference: Seq<B>,
        alternative: Seq<B>,
        multiallelic: bool,
        pass: bool,
        context: Option<Seq<B>>,
    ) -> Self {
        Self {
            chromosome,
            position,
            reference,
            alternative,
            multiallelic,
            pass,
            context,
        }
    }

    // --- Accessors (read-only) ---

    /// Returns the chromosome / contig name (e.g. `"chr1"`).
    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }

    /// Returns the 1-based start position of the mutation.
    ///
    /// This follows VCF conventions: the coordinate refers to the first base of `reference`.
    pub fn position(&self) -> Pos {
        self.position
    }

    /// Returns the reference allele sequence.
    pub fn reference(&self) -> &Seq<B> {
        &self.reference
    }

    /// Returns the alternative allele sequence.
    pub fn alternative(&self) -> &Seq<B> {
        &self.alternative
    }

    /// Returns whether this mutation originated from a multi-allelic site.
    ///
    /// This is metadata (e.g. a VCF record with multiple ALT alleles) and does not change
    /// coordinate or allele semantics.
    pub fn is_multiallelic(&self) -> bool {
        self.multiallelic
    }

    /// Returns whether this mutation passed upstream filtering/QC.
    ///
    /// In VCF terms, this commonly corresponds to the `FILTER` field being `PASS`
    /// (and sometimes `.` depending on the caller’s convention).
    ///
    /// This is metadata only; downstream code should decide how to handle non-passing
    /// mutations explicitly.
    pub fn is_pass(&self) -> bool {
        self.pass
    }

    /// Returns the optional sequence context for this mutation.
    ///
    /// This is typically used for contexts like trinucleotide sequence around the
    /// mutation site (computed later using a reference genome).
    ///
    /// Returns `None` if no context has been attached.
    pub fn context(&self) -> Option<&Seq<B>> {
        self.context.as_ref()
    }

    // --- Context mutation (write) ---

    /// Set (or overwrite) the context sequence for this mutation.
    ///
    /// This is intended for pipelines where the mutation is parsed first and context
    /// is computed later (e.g. reference lookup).
    ///
    /// This overwrites any existing context.
    pub fn set_context(&mut self, seq: Seq<B>) {
        self.context = Some(seq);
    }

    /// Clear any attached context sequence.
    pub fn clear_context(&mut self) {
        self.context = None;
    }

    // --- Computed Properties (read-only) ---

    /// Return the length of the reference allele in bases.
    ///
    /// This is a convenience wrapper around [`Seq::len`].
    pub fn reflen(&self) -> usize {
        self.reference.len()
    }

    /// Return the length of the alternative allele in bases.
    ///
    /// This is a convenience wrapper around [`Seq::len`].
    pub fn altlen(&self) -> usize {
        self.alternative.len()
    }
    /// Return the signed size change implied by this mutation.
    ///
    /// Defined as:
    /// ```text
    /// delta = alt_length - ref_length
    /// ```
    ///
    /// Interpretation:
    /// - `delta == 0` → equal-length substitution (SNV/DOUBLET/MNV)
    /// - `delta > 0`  → insertion (net gain of bases)
    /// - `delta < 0`  → deletion (net loss of bases)
    ///
    /// This is purely length-based and does not depend on allele normalization.
    pub fn delta(&self) -> i64 {
        self.altlen() as i64 - self.reflen() as i64
    }
    /// Return the mutation class derived from allele lengths.
    ///
    /// This is computed on demand using [`SmallMutationType::from_lengths`]
    /// and is therefore always consistent with the stored alleles.
    pub fn class(&self) -> SmallMutationType {
        SmallMutationType::from_lengths(self.reflen(), self.altlen())
    }

    /// Compute transition/transversion classification for this mutation.
    ///
    /// This classification is only defined for **single-nucleotide substitutions**
    /// ([`SmallMutationType::SNV`]). For all other mutation types, this returns `None`.
    ///
    /// This method is conservative in the presence of ambiguity:
    /// - If either base has ambiguous chemical class (e.g. `N`, `S`, `W`), returns `None`.
    ///
    /// # Returns
    /// - `Some(TiTv::Transition)` for A↔G or C↔T substitutions (including unambiguous
    ///   IUPAC codes that resolve to a single chemical class).
    /// - `Some(TiTv::Transversion)` for purine↔pyrimidine substitutions.
    /// - `None` if not an SNV or if ambiguity prevents a confident classification.
    pub fn titv(&self) -> Option<TiTv> {
        if self.class() != SmallMutationType::SNV {
            return None;
        }

        let r = self.reference.as_slice().first()?;
        let a = self.alternative.as_slice().first()?;

        TiTv::from_chemical_class(r.chemical_class(), a.chemical_class())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmallMutationType {
    SNV,
    DOUBLET,
    MNV,
    INSERTION,
    DELETION,
}

impl fmt::Display for SmallMutationType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            SmallMutationType::SNV => "SNV",
            SmallMutationType::DOUBLET => "DOUBLET",
            SmallMutationType::MNV => "MNV",
            SmallMutationType::INSERTION => "INSERTION",
            SmallMutationType::DELETION => "DELETION",
        };
        write!(f, "{s}")
    }
}

impl SmallMutationType {
    /// Classify a small variant based on reference and alternative allele lengths.
    ///
    /// Classification rules:
    /// - If `altlen > reflen` → [`SmallMutationType::INSERTION`]
    /// - If `altlen < reflen` → [`SmallMutationType::DELETION`]
    /// - If lengths are equal:
    ///   - `reflen == 1` → [`SmallMutationType::SNV`]
    ///   - `reflen == 2` → [`SmallMutationType::DOUBLET`]
    ///   - `reflen >= 3` → [`SmallMutationType::MNV`]
    ///
    ///
    /// # Notes
    /// - This function is purely *shape-based* and does not inspect sequence content.
    /// - A length of `0` is invalid for VCF alleles; this function currently maps
    ///   `reflen == 0 && altlen == 0` to [`SmallMutationType::MNV`] and assumes such cases are rejected
    ///   elsewhere.
    pub fn from_lengths(reflen: usize, altlen: usize) -> Self {
        match altlen.cmp(&reflen) {
            std::cmp::Ordering::Greater => Self::INSERTION,
            std::cmp::Ordering::Less => Self::DELETION,
            std::cmp::Ordering::Equal => match reflen {
                0 => Self::MNV, // OR HANDLE AS ERROR ELSEWHERE (0-LENGTH ALLELES ARE INVALID FOR VCF)
                1 => Self::SNV,
                2 => Self::DOUBLET,
                _ => Self::MNV,
            },
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TiTv {
    Transition,
    Transversion,
}

impl TiTv {
    /// Classify a substitution as a transition or transversion based on chemical class.
    ///
    /// Returns:
    /// - [`Some(TiTv::Transition)`] if both bases are certainly purines (A↔G) or both
    ///   are certainly pyrimidines (C↔T)
    /// - [`Some(TiTv::Transversion)`] if the substitution is purine↔pyrimidine
    /// - `None` if either input class is ambiguous (e.g. derived from an ambiguous base)
    ///
    /// This conservative behavior matches the library philosophy: ambiguous inputs
    /// should not be silently guessed.
    pub fn from_chemical_class(reference: ChemClass, alternative: ChemClass) -> Option<TiTv> {
        match (reference, alternative) {
            (ChemClass::Purine, ChemClass::Purine) => Some(TiTv::Transition),
            (ChemClass::Pyrimidine, ChemClass::Pyrimidine) => Some(TiTv::Transition),
            (ChemClass::Purine, ChemClass::Pyrimidine) => Some(TiTv::Transversion),
            (ChemClass::Pyrimidine, ChemClass::Purine) => Some(TiTv::Transversion),
            _ => None, // If either chemical class is ambiguous, return None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{ChemClass, DnaBase, RnaBase};
    use crate::sequence::Seq;

    // --- Helpers ---

    fn dna(s: &str) -> Seq<DnaBase> {
        Seq::<DnaBase>::new(s).unwrap()
    }

    fn rna(s: &str) -> Seq<RnaBase> {
        Seq::<RnaBase>::new(s).unwrap()
    }

    fn dna_mut(ref_allele: &str, alt_allele: &str) -> DnaSmallMutation {
        SmallMutation::new(
            "chr1".to_string(),
            Pos::new(123).unwrap(),
            dna(ref_allele),
            dna(alt_allele),
            false,
            true,
            None,
        )
    }

    fn rna_mut(ref_allele: &str, alt_allele: &str) -> RnaSmallMutation {
        SmallMutation::new(
            "tx1".to_string(),
            Pos::new(7).unwrap(),
            rna(ref_allele),
            rna(alt_allele),
            false,
            true,
            None,
        )
    }

    // --- Construction / context ---

    #[test]
    fn new_sets_fields_and_lengths_dna() {
        let m = dna_mut("A", "G");
        assert_eq!(m.reflen(), 1);
        assert_eq!(m.altlen(), 1);
        assert_eq!(m.delta(), 0);
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), Some(TiTv::Transition));
    }

    #[test]
    fn set_context_sets_context() {
        let mut m = dna_mut("A", "G");
        assert_eq!(m.context, None);

        m.set_context(dna("TCA"));
        assert_eq!(m.context, Some(dna("TCA")));

        // overwrite
        m.set_context(dna("AAA"));
        assert_eq!(m.context, Some(dna("AAA")));
    }

    // --- SmallMutationType::from_lengths ---

    #[test]
    fn from_lengths_classifies_equal_length_substitutions() {
        assert_eq!(
            SmallMutationType::from_lengths(1, 1),
            SmallMutationType::SNV
        );
        assert_eq!(
            SmallMutationType::from_lengths(2, 2),
            SmallMutationType::DOUBLET
        );
        assert_eq!(
            SmallMutationType::from_lengths(3, 3),
            SmallMutationType::MNV
        );
        assert_eq!(
            SmallMutationType::from_lengths(10, 10),
            SmallMutationType::MNV
        );
    }

    #[test]
    fn from_lengths_classifies_indels() {
        assert_eq!(
            SmallMutationType::from_lengths(1, 2),
            SmallMutationType::INSERTION
        );
        assert_eq!(
            SmallMutationType::from_lengths(2, 1),
            SmallMutationType::DELETION
        );
        assert_eq!(
            SmallMutationType::from_lengths(5, 9),
            SmallMutationType::INSERTION
        );
        assert_eq!(
            SmallMutationType::from_lengths(9, 5),
            SmallMutationType::DELETION
        );
    }

    // --- delta / class integration tests ---

    #[test]
    fn class_and_delta_match_expected_for_snv_doublet_mnv() {
        let snv = dna_mut("A", "C");
        assert_eq!(snv.class(), SmallMutationType::SNV);
        assert_eq!(snv.delta(), 0);

        let dbl = dna_mut("AC", "GT");
        assert_eq!(dbl.class(), SmallMutationType::DOUBLET);
        assert_eq!(dbl.delta(), 0);

        let mnv = dna_mut("ACG", "TTA");
        assert_eq!(mnv.class(), SmallMutationType::MNV);
        assert_eq!(mnv.delta(), 0);
    }

    #[test]
    fn class_and_delta_match_expected_for_insertion_and_deletion() {
        let ins = dna_mut("A", "AT");
        assert_eq!(ins.class(), SmallMutationType::INSERTION);
        assert_eq!(ins.delta(), 1);

        let del = dna_mut("AT", "A");
        assert_eq!(del.class(), SmallMutationType::DELETION);
        assert_eq!(del.delta(), -1);
    }

    // --- TiTv::from_chemical_class ---

    #[test]
    fn titv_from_chemical_class_transition_and_transversion() {
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Purine, ChemClass::Purine),
            Some(TiTv::Transition)
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Pyrimidine, ChemClass::Pyrimidine),
            Some(TiTv::Transition)
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Purine, ChemClass::Pyrimidine),
            Some(TiTv::Transversion)
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Pyrimidine, ChemClass::Purine),
            Some(TiTv::Transversion)
        );
    }

    #[test]
    fn titv_from_chemical_class_ambiguous_returns_none() {
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Ambiguous, ChemClass::Purine),
            None
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Pyrimidine, ChemClass::Ambiguous),
            None
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Ambiguous, ChemClass::Ambiguous),
            None
        );
    }

    // --- SmallMutation::titv ---

    #[test]
    fn titv_only_defined_for_snvs() {
        let mnv = dna_mut("AC", "GT");
        assert_eq!(mnv.class(), SmallMutationType::DOUBLET);
        assert_eq!(mnv.titv(), None);

        let ins = dna_mut("A", "AT");
        assert_eq!(ins.class(), SmallMutationType::INSERTION);
        assert_eq!(ins.titv(), None);

        let del = dna_mut("AT", "A");
        assert_eq!(del.class(), SmallMutationType::DELETION);
        assert_eq!(del.titv(), None);
    }

    #[test]
    fn titv_transition_examples_dna() {
        // A <-> G is a transition (purine <-> purine)
        let m = dna_mut("A", "G");
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), Some(TiTv::Transition));

        // C <-> T is a transition (pyrimidine <-> pyrimidine)
        let m2 = dna_mut("C", "T");
        assert_eq!(m2.titv(), Some(TiTv::Transition));
    }

    #[test]
    fn titv_transversion_examples_dna() {
        // A <-> C is a transversion (purine <-> pyrimidine)
        let m = dna_mut("A", "C");
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), Some(TiTv::Transversion));

        // G <-> T is a transversion
        let m2 = dna_mut("G", "T");
        assert_eq!(m2.titv(), Some(TiTv::Transversion));
    }

    #[test]
    fn titv_returns_none_when_ambiguous_base_present() {
        // N is ChemClass::Ambiguous in your design
        let m = dna_mut("N", "A");
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), None);

        let m2 = dna_mut("A", "N");
        assert_eq!(m2.class(), SmallMutationType::SNV);
        assert_eq!(m2.titv(), None);
    }

    // --- Display ---

    #[test]
    fn display_includes_core_fields() {
        let m = dna_mut("A", "G");
        let s = m.to_string();

        // Keep this intentionally loose so formatting tweaks don't require rewrites.
        assert!(s.contains("chr1:123"));
        assert!(s.contains("A>G"));
        assert!(s.contains("delta: 0"));
        assert!(s.contains("class: SNV"));
        assert!(s.contains("multiallelic:false"));
    }

    // --- RNA smoke tests (generic over Base works) ---

    #[test]
    fn rna_small_mutation_works_and_classifies() {
        let m = rna_mut("A", "G");
        assert_eq!(m.reflen(), 1);
        assert_eq!(m.altlen(), 1);
        assert_eq!(m.delta(), 0);
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), Some(TiTv::Transition)); // purine<->purine is transition

        let ins = rna_mut("A", "AU");
        assert_eq!(ins.class(), SmallMutationType::INSERTION);
        assert_eq!(ins.delta(), 1);
        assert_eq!(ins.titv(), None);
    }

    #[test]
    fn rna_ambiguous_titv_none() {
        // 'N' exists for RNA alphabet in your base impl
        let m = rna_mut("N", "A");
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), None);
    }
}
