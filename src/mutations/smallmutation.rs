use core::fmt;

use crate::{
    base::{Base, ChemClass, DnaBase, IupacDnaBase, IupacRnaBase, RnaBase},
    coords::{BasePos, InterbaseInterval, InterbasePos, Region, Strand},
    error::MutationError as Error,
    sequences::Seq,
};

pub(crate) type Result<T> = std::result::Result<T, Error>;

pub type IupacDnaSmallMutation = SmallMutation<IupacDnaBase>;
pub type IupacRnaSmallMutation = SmallMutation<IupacRnaBase>;
pub type DnaSmallMutation = SmallMutation<DnaBase>;
pub type RnaSmallMutation = SmallMutation<RnaBase>;

/// A small mutation (SNV/MNV/indel) over a specific nucleotide alphabet `B`.
/// - `SmallMutation<DnaBase>` is DNA
/// - `SmallMutation<RnaBase>` is RNA
///
/// ## Coordinate and allele semantics
/// - `location` in reference sequence affected by this mutation. **0-start interbase coordinates** (VRS-style).
/// - `reference` and `alternative` are stored **as provided** by the caller.
///   No left/right trimming, normalization, or decomposition is performed.
///
/// ## Invariants (enforced by constructor)
/// Length of reference sequence must be equal to the size of affected region according to
/// `location`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmallMutation<B: Base> {
    /// 0-based interbase interval describing the set of bases
    location: Region,

    /// Original base sequence in reference. Length must be equal to
    reference: Seq<B>,

    /// Mutation of sequence
    alternative: Seq<B>,

    /// If mutation was on a double stranded molecule like DNA
    /// then field indicates the strand of the reference genome the reference allele was on.
    /// If mutation comes from an unstranded molecule (single stranded DNA / RNA)
    /// or strand is not known can be set to None.
    strand: Option<Strand>,
}

// Implement the `fmt::Display` trait for `SmallMutation`.
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
            "{}:{}-{} {}>{} (strand: {}, delta: {}; class: {})",
            self.chromosome(),
            self.start(),
            self.end(),
            self.reference,
            self.alternative,
            self.strand
                .map(|opt| opt.to_string())
                .unwrap_or("None".to_string()),
            self.delta(),
            self.class(),
        )
    }
}

impl<B: Base> SmallMutation<B> {
    /// Creates a [`SmallMutation`]
    ///
    /// # Errors
    /// Returns Error::MutationReferenceLengthDoesNotMatchInterval if length of interval in location
    /// is not equal to the size of the reference allele sequence.
    ///
    /// # Panics
    /// Any panics in this function represent a bug in the codebase. Please create a new github issue.
    ///
    pub fn new(
        location: Region,
        reference: Seq<B>,
        alternative: Seq<B>,
        strand: Option<Strand>,
    ) -> Result<Self> {
        // Return Error if interval length does not match reference sequence length
        if location.interval().len() != reference.len() {
            return Err(Error::MutationReferenceLengthDoesNotMatchInterval {
                interval_length: location.interval().len(),
                reference_length: reference.len(),
                variant_id: format!(
                    "{}:{}-{} {}>{}",
                    location.name(),
                    location.interval().start(),
                    location.interval().end(),
                    reference,
                    alternative
                ),
            });
        }

        // Otherwise create mutation object
        Ok(Self {
            location,
            reference,
            alternative,
            strand,
        })
    }

    /// Construct a new [`SmallMutation`].
    ///
    /// This constructor does **not** attempt to validate biological correctness beyond
    /// what is already guaranteed by [`Seq<B>`] (i.e. sequences conform to the alphabet).
    ///
    /// In particular, this type:
    /// - assumes `position` is **1-based** (VCF-style coordinates)
    /// - stores `reference` and `alternative` as provided (no trimming/normalization)
    ///
    /// If you need allele normalization (left/right trimming of shared prefix/suffix),
    /// do it before constructing this type.
    ///
    /// # Parameters
    /// - `chromosome`: reference sequence / contig name (e.g. `"chr1"`)
    /// - `position`: 1-based start coordinate ([`BasePos`])
    /// - `reference`: reference allele sequence.
    /// - `alternative`: alternative allele sequence
    /// - `strand`: Strand
    ///
    /// # Errors
    /// None. Correctness is guaranteed by argument types.
    ///
    /// # Panics
    /// All panics represent bugs in seqlib itself and should be reported as github issues.
    pub fn new_from_vcflike(
        chromosome: String,
        position: BasePos,
        reference: Seq<B>,
        alternative: Seq<B>,
        strand: Option<Strand>,
    ) -> Result<Self> {
        let start = usize::from(position) - 1usize; // Subtraction should never fail since BasePos is a NonZeroUsize
        let end = start + reference.len();

        let interval =
            InterbaseInterval::new(InterbasePos::new(start), InterbasePos::new(end))
            .unwrap_or_else(|err| {
                panic!(
                "Unexpected seqlib failure converting a VCF-like mutation. \
                end position <= start position - his should not be possible and represents a bug in seqlib. \n\
                Error: {err}")
            });
        let region = Region::new(chromosome, interval);

        Ok(
            Self::new(region, reference, alternative, strand).unwrap_or_else(|err| {
                panic!(
                    "Unexpected seqlib failure converting a VCF-like mutation. \
                     Please report this as a GitHub issue.\n\
                     Error: {err}"
                )
            }),
        )
    }

    // --- Accessors (read-only) ---
    /// Returns the chromosome / contig name (e.g. `"chr1"`).
    pub fn chromosome(&self) -> &str {
        self.location.name()
    }

    /// Returns the reference allele sequence.
    pub fn reference(&self) -> &Seq<B> {
        &self.reference
    }

    /// Returns the alternative allele sequence.
    pub fn alternative(&self) -> &Seq<B> {
        &self.alternative
    }

    /// Returns the strand of the source molecule that contains
    /// the reference allele of this mutation (or None if from  
    /// a single stranded molecule like RNA
    pub fn strand(&self) -> Option<&Strand> {
        self.strand.as_ref()
    }

    /// Returns the start position of the mutated interval
    pub fn start(&self) -> &InterbasePos {
        self.location.interval().start()
    }

    /// Returns the end position of the mutated interval
    pub fn end(&self) -> &InterbasePos {
        self.location.interval().end()
    }

    pub fn location(&self) -> &Region {
        &self.location
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

        let r_chemical_class = r.try_chemical_class().ok()?;
        let a_chemical_class = a.try_chemical_class().ok()?;

        let titv = TiTv::from_chemical_class(r_chemical_class, a_chemical_class);
        Some(titv)
    }

    /// Render a compact `chrom:pos REF>ALT` representation of the mutation.
    ///
    /// This is a convenience formatter intended for **human-readable output**
    /// (e.g. logs, debugging, CLI display). It mirrors the common genomics shorthand
    /// used in VCF-style contexts, but is **not** a stable serialization format.
    ///
    /// ## Format
    ///
    /// ```text
    /// chromosome:position REF>ALT
    /// ```
    ///
    /// Where:
    /// - `chromosome` is the contig / reference name (e.g. `"chr1"`)
    /// - `position` is the **1-based** start coordinate (VCF-style)
    /// - `REF` is the reference allele sequence
    /// - `ALT` is the alternative allele sequence
    ///
    /// Alleles are rendered using their sequence display implementations
    /// (uppercase IUPAC symbols).
    ///
    /// ## Examples
    ///
    /// ```rust
    /// use seqlib::mutations::DnaSmallMutation;
    /// use seqlib::sequences::DnaSeq;
    /// use seqlib::coords::{Region, Strand};
    ///
    /// let m = DnaSmallMutation::new(
    ///     Region::new("chr1", 122..123),
    ///     DnaSeq::new("A").unwrap(),
    ///     DnaSeq::new("G").unwrap(),
    ///     Some(Strand::Positive),
    /// ).unwrap();
    ///
    /// assert_eq!(m.chrom_start_end_ref_alt(), "chr1:122-123 A>G");
    /// ```
    pub fn chrom_start_end_ref_alt(&self) -> String {
        format!(
            "{}:{}-{} {}>{}",
            self.chromosome(),
            self.start(),
            self.end(),
            self.reference(),
            self.alternative()
        )
    }

    /// Render `chromosome:position REF>ALT` with alleles coloured based on their default ANSI
    /// style (based on sequence)
    pub fn format_with_colour(&self) -> String {
        format!(
            "{}:{}-{} {}>{}",
            self.chromosome(),
            self.start(),
            self.end(),
            self.reference().format_with_colour(),
            self.alternative().format_with_colour()
        )
    }

    /// Reverse complement mutation
    /// Reverse complements reference and alternative sequence and flips strand
    /// field if present.
    pub fn reverse_complement(&self) -> Self {
        // Flip Strand
        let newstrand = self.strand.map(|strand| strand.flip());

        Self {
            location: self.location().clone(),
            reference: self.reference.reverse_complement(),
            alternative: self.alternative.reverse_complement(),
            strand: newstrand,
        }
    }

    /// Ensure middle base of reference sequence is a pyrimidine
    /// This is accomplished by reverse complementing the mutation
    /// if there is a purine centered.
    ///
    /// # Errors
    /// If the middlebase is ambiguous, or the reference sequence is even in length (has no
    /// middlebase) will return an Error
    pub fn try_pyrimidine_center(&self) -> Result<Self> {
        let refseq = &self.reference;

        let Some(middlebase) = refseq.middlebase() else {
            return Err(Error::MissingMiddleBase);
        };

        let middle_chemclass =
            middlebase
                .try_chemical_class()
                .map_err(|_| Error::AmbiguousMiddleBase {
                    base: middlebase.to_char(),
                })?;

        // Flip Chemclass
        match middle_chemclass {
            ChemClass::Purine => Ok(self.reverse_complement()),
            ChemClass::Pyrimidine => Ok(self.clone()),
        }
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
    /// - This function is purely *length-based* and does not inspect sequence content.
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
    /// Purines -> Pyrimidine = Transversion
    /// Pyrimidine -> Purine = Transversion
    /// Purine -> Purine = Transition
    /// Pyrimidine -> Pyrimidine = Transition
    pub fn from_chemical_class(reference: ChemClass, alternative: ChemClass) -> TiTv {
        match (reference, alternative) {
            (ChemClass::Purine, ChemClass::Purine) => TiTv::Transition,
            (ChemClass::Pyrimidine, ChemClass::Pyrimidine) => TiTv::Transition,
            (ChemClass::Purine, ChemClass::Pyrimidine) => TiTv::Transversion,
            (ChemClass::Pyrimidine, ChemClass::Purine) => TiTv::Transversion,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{ChemClass, IupacDnaBase, IupacRnaBase};
    use crate::sequences::Seq;

    // --- Helpers ---

    fn dna_iupac(s: &str) -> Seq<IupacDnaBase> {
        Seq::<IupacDnaBase>::new(s).unwrap()
    }

    fn _dna(s: &str) -> Seq<DnaBase> {
        Seq::<DnaBase>::new(s).unwrap()
    }

    fn rna_iupac(s: &str) -> Seq<IupacRnaBase> {
        Seq::<IupacRnaBase>::new(s).unwrap()
    }

    fn dna_mut(ref_allele: &str, alt_allele: &str) -> IupacDnaSmallMutation {
        let reference = dna_iupac(ref_allele);
        IupacDnaSmallMutation::new(
            Region::new("chr1", 122..122 + reference.len()),
            reference,
            dna_iupac(alt_allele),
            Some(Strand::Positive),
        )
        .unwrap()
    }

    fn rna_mut(ref_allele: &str, alt_allele: &str) -> IupacRnaSmallMutation {
        let reference = rna_iupac(ref_allele);
        IupacRnaSmallMutation::new(
            Region::new("tx1", 0..reference.len()),
            reference,
            rna_iupac(alt_allele),
            Some(Strand::Positive),
        )
        .unwrap()
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
    fn new_preserves_region_and_interbase_accessors() {
        let m = dna_mut("ACG", "TTA");
        let region = Region::new("chr1", 122..125);

        assert_eq!(m.location(), &region);
        assert_eq!(m.chromosome(), "chr1");
        assert_eq!(*m.start(), InterbasePos::new(122));
        assert_eq!(*m.end(), InterbasePos::new(125));
        assert_eq!(m.reference(), &dna_iupac("ACG"));
        assert_eq!(m.alternative(), &dna_iupac("TTA"));
        assert_eq!(m.strand(), Some(&Strand::Positive));
    }

    #[test]
    fn new_rejects_reference_and_interval_length_mismatches() {
        for end in [123, 125] {
            let region = Region::new("chr1", 122..end);
            let err = IupacDnaSmallMutation::new(region, dna_iupac("AC"), dna_iupac("GT"), None)
                .unwrap_err();

            assert_eq!(
                err,
                Error::MutationReferenceLengthDoesNotMatchInterval {
                    interval_length: end - 122,
                    reference_length: 2,
                    variant_id: format!("chr1:122-{end} AC>GT"),
                }
            );
        }
    }

    #[test]
    fn new_accepts_empty_reference_at_interbase_boundary() {
        let region = Region::new("chr1", 0..0);
        let m = IupacDnaSmallMutation::new(region.clone(), dna_iupac(""), dna_iupac("A"), None)
            .unwrap();

        assert_eq!(m.location(), &region);
        assert_eq!(m.reflen(), 0);
        assert_eq!(m.altlen(), 1);
        assert_eq!(m.class(), SmallMutationType::INSERTION);
        assert_eq!(m.chrom_start_end_ref_alt(), "chr1:0-0 >A");
    }

    #[test]
    fn new_from_vcflike_converts_start_and_reference_span() {
        let first = IupacDnaSmallMutation::new_from_vcflike(
            "chr1".to_string(),
            BasePos::new(1).unwrap(),
            dna_iupac("A"),
            dna_iupac("G"),
            None,
        )
        .unwrap();
        assert_eq!(first.location(), &Region::new("chr1", 0..1));

        let multibase = IupacDnaSmallMutation::new_from_vcflike(
            "chr2".to_string(),
            BasePos::new(123).unwrap(),
            dna_iupac("ACG"),
            dna_iupac("TTA"),
            Some(Strand::Negative),
        )
        .unwrap();
        assert_eq!(multibase.location(), &Region::new("chr2", 122..125));
        assert_eq!(multibase.reference(), &dna_iupac("ACG"));
        assert_eq!(multibase.alternative(), &dna_iupac("TTA"));
        assert_eq!(multibase.strand(), Some(&Strand::Negative));
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
            TiTv::Transition
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Pyrimidine, ChemClass::Pyrimidine),
            TiTv::Transition
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Purine, ChemClass::Pyrimidine),
            TiTv::Transversion
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Pyrimidine, ChemClass::Purine),
            TiTv::Transversion
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
        assert!(s.contains("chr1:122-123"));
        assert!(s.contains("A>G"));
        assert!(s.contains("delta: 0"));
        assert_eq!(m.chrom_start_end_ref_alt(), "chr1:122-123 A>G");
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
