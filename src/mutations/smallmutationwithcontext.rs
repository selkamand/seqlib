use crate::{
    base::Base,
    coords::{BaseInterval, BasePos},
    error::MutationError,
    mutations::SmallMutation,
    render::SeqStyler,
    sequences::{Seq, SourcedSeq},
};

pub(crate) type Result<T> = std::result::Result<T, MutationError>;

/// A small mutation annotated with the reference sequence from which it was created
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MutationWithContext<B: Base> {
    mutation: SmallMutation<B>,
    context: SourcedSeq<B>,
}

impl<B: Base> std::fmt::Display for MutationWithContext<B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "-- Mutation --")?;
        writeln!(f, "{}", self.mutation)?;
        writeln!(f)?;
        writeln!(f, "-- Context --")?;
        writeln!(f, "{}", self.context)
    }
}

impl<B: Base> MutationWithContext<B> {
    /// Construct a new [`MutationWithContext`]
    ///
    /// Describing a [`SmallMutation`] and its local sequence context.
    /// Mutation must occur within the region described by context.
    /// If strand is supplied, it must be supplied for both the mutation and the context.
    ///
    /// Note if strand is supplied for both but the value itself differs, the context will be
    /// reverse complemented so they match.
    ///
    /// ## Errors
    /// - Returns a [`MutationError::MismatchedStrandOption`] error if a strand is specified for mutation or context but not both.
    /// - Returns a [`MutationError::MismatchedChromosomeName`] if chromosome names don't match.
    /// - Returns a [`MutationError::MutationPositionOutsideInterval`] if mutation position is outside the region described by the context.
    pub fn new(mutation: SmallMutation<B>, context: SourcedSeq<B>) -> Result<Self> {
        // Strand comparison
        let mutstrand = mutation.strand();
        let contextstrand = context.strand();

        //TODO: Because this function takes ownership of context, we could use
        //reverse_complement_in_place to make more memory efficient but because
        //contexts will usually be quite short and very few in memory at any one time we will
        //keep the copy-on-modify approach here for now
        let normalised_context = match (mutstrand, contextstrand) {
            (None, None) => context,

            // If strand is given for mutation but not for context sequence or vice-versa throw an error
            (None, Some(_)) => return Err(MutationError::MismatchedStrandOption),
            (Some(_), None) => return Err(MutationError::MismatchedStrandOption),

            // Reverse complement context to ensure strand of context and mutation always match
            (Some(strand1), Some(strand2)) => match strand1 == strand2 {
                true => context,
                false => context.reverse_complement(),
            },
        };

        // Create struct to gain access to helpers to check validity
        let mutwithcontext = Self {
            mutation,
            context: normalised_context,
        };

        // Enforce invariants
        mutwithcontext.check_chromosome_names_match()?;
        mutwithcontext.check_mutation_position_within_sequence_context_interval()?;
        mutwithcontext.check_reference_bases_viable()?;

        // If all pass return context
        Ok(mutwithcontext)
    }

    // < Validity Checks>

    /// Returns an [`MutationError::MismatchedChromosomeName`] if chromosome names don't match
    fn check_chromosome_names_match(&self) -> Result<()> {
        let mutated_chromosome = self.mutation().chromosome();
        let context_chromosome = self.context().region().name();

        if mutated_chromosome != context_chromosome {
            return Err(MutationError::MismatchedChromosomeName {
                mutated_chromosome: mutated_chromosome.into(),
                context_chromosome: context_chromosome.into(),
            });
        }

        Ok(())
    }

    /// Returns an [`MutationError::MutationPositionOutsideInterval`] if mutation position is
    /// outside the region described by the context
    fn check_mutation_position_within_sequence_context_interval(&self) -> Result<()> {
        let mutation_position = self.mutation().position();
        let region = self.context().region();
        let context_start = region.interval().start();
        let context_end = region.interval().end();
        // If mutation is outside the bounds of the sequence context return an error
        if mutation_position < *context_start || mutation_position > *context_end {
            return Err(MutationError::MutationPositionOutsideInterval {
                position: mutation_position.get(),
                start: context_start.get(),
                end: context_end.get(),
            });
        }

        Ok(())
    }

    fn check_reference_bases_viable(&self) -> Result<()> {
        if !self.reference_bases_viable() {
            return Err(MutationError::MismatchedReferenceAlleleAndContextSeq {
                mutation: self.mutation().chrom_pos_ref_alt(),
                context: self
                    .context()
                    .seq()
                    .format_with_highlighted_base_interval(Some(&self.mutated_interval())),
            });
        }
        Ok(())
    }
    // < Getters >

    /// Get context sequence local 1-based start position of mutation.
    ///
    /// For example, anchor of Pos(1) means the mutation affect the first base of sequence context.
    /// See [`MutationWithContext::mutated_interval`] to get the full interval of bases in sequence
    /// context that are mutated.
    ///
    /// # Panics
    /// Panics only if the invariants established by [`MutationWithContext::new`]
    /// have been violated, which indicates an internal bug.
    pub fn anchor(&self) -> BasePos {
        let mutation_position = self.mutation().position().get();
        let context_start = self.context().region().interval().start().get();

        let local_position = mutation_position
            .checked_sub(context_start)
            .and_then(|offset| offset.checked_add(1))
            .expect(
                "MutationWithContext invariant violated: \
             mutation position must lie within the context region",
            );

        BasePos::new(local_position).expect(
            "MutationWithContext invariant violated: \
         a validated local mutation position must be non-zero",
        )
    }

    /// Get the mutation
    pub fn mutation(&self) -> &SmallMutation<B> {
        &self.mutation
    }

    /// Get the context of the mutation
    pub fn context(&self) -> &SourcedSeq<B> {
        &self.context
    }

    /// Get the interval of context affected by the mutation
    ///
    ///
    /// Returns the Interval representing the mutated position (with respect to the context
    /// sequence)
    ///
    /// # Panics:
    /// Panics only if the invariants established by [`MutationWithContext::new`]
    /// have been violated or mutation position + reflength offset creates an overflow error, which indicate internal bugs.
    pub fn mutated_interval(&self) -> BaseInterval {
        let start = self.anchor();

        let offset = self
            .mutation()
            .reflen()
            .checked_sub(1)
            .expect("reference alleles must contain at least one base");

        let end = start
            .try_add(offset)
            .expect("mutation interval position overflowed usize");

        BaseInterval::new(start, end)
            .expect("start plus a non-negative offset must form a valid interval")
    }

    /// Is the context sequence viable given the mutation reference allele and anchor position.
    /// Returns true if the context sequence matchs the mutation reference sequence at the expected region
    /// (based on anchor position and reference size)
    /// Returns false in all other cases
    fn reference_bases_viable(&self) -> bool {
        let interval = self.mutated_interval();

        let Ok(mutated_bases) = self.context().seq().subseq_by_base_interval(&interval) else {
            return false;
        };

        mutated_bases == *self.mutation().reference()
    }

    /// Apply a mutation to a sequence context
    ///
    /// # Panics
    /// Panics only if the invariants established by [`MutationWithContext::new`] which indicates
    /// internal bugs.
    pub fn apply_mutation(&self) -> Seq<B> {
        let mut seq = self.context().seq().clone();

        let interval = self.mutated_interval();

        seq.mutate(interval, self.mutation().alternative()).expect(
            "any mutation in mutation_with_context to be possible to apply to the reference seq",
        )
    }

    /// Apply a mutation to a sequence context and return the interval occupied by ALT bases.
    ///
    /// The returned interval uses 1-based coordinates local to the mutated sequence.
    /// Mutations with an empty ALT allele return `None` because no ALT bases are present in the
    /// mutated sequence.
    ///
    /// # Panics
    /// Panics only if the invariants established by [`MutationWithContext::new`] which indicates
    /// internal bugs.
    pub fn apply_mutation_with_alt_interval(&self) -> (Seq<B>, Option<BaseInterval>) {
        let mut seq = self.context().seq().clone();

        let interval = self.mutated_interval();

        let mutated = seq.mutate(interval, self.mutation().alternative()).expect(
            "any mutation in mutation_with_context to be possible to apply to the reference seq",
        );

        let alt_interval = if self.mutation().altlen() == 0 {
            None
        } else {
            let start = self.anchor();
            let end = start
                .try_add(
                    self.mutation()
                        .altlen()
                        .checked_sub(1)
                        .expect("non-empty alternative alleles must contain at least one base"),
                )
                .expect("mutated interval position overflowed usize");

            Some(
                BaseInterval::new(start, end)
                    .expect("start plus a non-negative offset must form a valid interval"),
            )
        };

        (mutated, alt_interval)
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
    /// - there is insufficient flanking sequence on either side of the anchor
    ///
    /// ## Notes
    /// - No mutation-type checks are performed (e.g. SNV vs indel)
    /// - No strand normalization or reverse-complementing is applied
    /// - The returned slice is borrowed; no allocation or copying occurs
    pub fn kmer_centered_on_anchor(&self, k: usize) -> Option<&[B]> {
        if k == 0 || k.is_multiple_of(2) {
            return None;
        }
        let center = self.anchor().as_0based_index();
        let half = k / 2;
        let start = center.checked_sub(half)?;
        let end = center + half + 1;
        if end > self.context().seq().len() {
            return None;
        }
        // subseq_slice returns Result<&[B]>; convert to Option
        self.context().seq().slice(start, end).ok()
    }

    /// Visualise the mutation in context
    ///
    /// creates a string with two lines. One with reference sequence (brakets indicating mutated
    /// interval/position) and one with the alternative sequence
    ///
    /// If there were errors applying the mutation to the sequence context, alternative will just
    /// display the error message.
    pub fn to_difference_string(&self) -> String {
        // Grab context
        let ctx = self.context();

        // Create formatted reference string
        let refstring = ctx
            .seq()
            .format_with_highlighted_base_interval(Some(&self.mutated_interval()));

        // Create formatted altstring
        let altstring = self.apply_mutation().to_string();

        format!("{refstring}\n{altstring}")
    }

    /// Create a string representing context sequence before mutation, with mutated region coloured
    /// using ansi codes.
    pub fn format_context_sequence_and_highlight_mutated_bases(&self) -> String {
        self.context()
            .seq()
            .format_with_coloured_interval(&self.mutated_interval())
    }
    pub fn format_mutated_sequence_and_highlight_mutated_bases(&self) -> String {
        let (mutated, alt_interval) = self.apply_mutation_with_alt_interval();

        if let Some(interval) = alt_interval {
            mutated.format_with_coloured_interval(&interval)
        } else {
            mutated.format_with_style(&SeqStyler::DIMMED)
        }
    }

    /// Render the reference context and mutated context as two ANSI-coloured lines.
    pub fn format_with_colour(&self) -> String {
        format!(
            "Mutation [{}]\nreference: {}\nmutated:   {}",
            self.mutation.format_with_colour(),
            self.format_context_sequence_and_highlight_mutated_bases(),
            self.format_mutated_sequence_and_highlight_mutated_bases()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        base::DnaBase,
        coords::{Region, Strand},
        sequences::DnaSeq,
    };

    fn dna_mutation_with_context(
        position: usize,
        reference: &str,
        alternative: &str,
    ) -> MutationWithContext<DnaBase> {
        let mutation = SmallMutation::new(
            "chr1".to_string(),
            BasePos::new(position).unwrap(),
            DnaSeq::new(reference).unwrap(),
            DnaSeq::new(alternative).unwrap(),
            Some(Strand::Positive),
        );
        let context = SourcedSeq::new(
            DnaSeq::new("ACGTACGT").unwrap(),
            Region::new(
                "chr1",
                BaseInterval::new(BasePos::new(100).unwrap(), BasePos::new(107).unwrap()).unwrap(),
            ),
            Some(Strand::Positive),
        );

        MutationWithContext::new(mutation, context).unwrap()
    }

    #[test]
    fn format_mutated_sequence_highlights_snv_alt_base() {
        let mutation = dna_mutation_with_context(103, "T", "A");
        let mutated = DnaSeq::new("ACGAACGT").unwrap();
        let interval =
            BaseInterval::new(BasePos::new(4).unwrap(), BasePos::new(4).unwrap()).unwrap();
        let (applied, alt_interval) = mutation.apply_mutation_with_alt_interval();

        assert_eq!(applied, mutated);
        assert_eq!(alt_interval, Some(interval.clone()));

        assert_eq!(
            mutation.format_mutated_sequence_and_highlight_mutated_bases(),
            mutated.format_with_coloured_interval(&interval)
        );
    }

    #[test]
    fn format_mutated_sequence_highlights_mnv_alt_bases() {
        let mutation = dna_mutation_with_context(103, "TAC", "GCA");
        let mutated = DnaSeq::new("ACGGCAGT").unwrap();
        let interval =
            BaseInterval::new(BasePos::new(4).unwrap(), BasePos::new(6).unwrap()).unwrap();
        let (applied, alt_interval) = mutation.apply_mutation_with_alt_interval();

        assert_eq!(applied, mutated);
        assert_eq!(alt_interval, Some(interval.clone()));

        assert_eq!(
            mutation.format_mutated_sequence_and_highlight_mutated_bases(),
            mutated.format_with_coloured_interval(&interval)
        );
    }

    #[test]
    fn format_mutated_sequence_highlights_entire_insertion_alt_allele() {
        let mutation = dna_mutation_with_context(103, "T", "TGG");
        let mutated = DnaSeq::new("ACGTGGACGT").unwrap();
        let interval =
            BaseInterval::new(BasePos::new(4).unwrap(), BasePos::new(6).unwrap()).unwrap();
        let (applied, alt_interval) = mutation.apply_mutation_with_alt_interval();

        assert_eq!(applied, mutated);
        assert_eq!(alt_interval, Some(interval.clone()));

        assert_eq!(
            mutation.format_mutated_sequence_and_highlight_mutated_bases(),
            mutated.format_with_coloured_interval(&interval)
        );
    }

    #[test]
    fn format_mutated_sequence_highlights_non_empty_deletion_alt_allele() {
        let mutation = dna_mutation_with_context(103, "TAC", "T");
        let mutated = DnaSeq::new("ACGTGT").unwrap();
        let interval =
            BaseInterval::new(BasePos::new(4).unwrap(), BasePos::new(4).unwrap()).unwrap();
        let (applied, alt_interval) = mutation.apply_mutation_with_alt_interval();

        assert_eq!(applied, mutated);
        assert_eq!(alt_interval, Some(interval.clone()));

        assert_eq!(
            mutation.format_mutated_sequence_and_highlight_mutated_bases(),
            mutated.format_with_coloured_interval(&interval)
        );
    }

    #[test]
    fn format_mutated_sequence_dims_zero_length_alt_without_highlighted_bases() {
        let mutation = dna_mutation_with_context(103, "TAC", "");
        let mutated = DnaSeq::new("ACGGT").unwrap();
        let (applied, alt_interval) = mutation.apply_mutation_with_alt_interval();

        assert_eq!(applied, mutated);
        assert_eq!(alt_interval, None);

        assert_eq!(
            mutation.format_mutated_sequence_and_highlight_mutated_bases(),
            mutated.format_with_style(&SeqStyler::DIMMED)
        );
    }

    #[test]
    fn format_mutated_sequence_highlights_shorter_replacement_alt_allele() {
        let mutation = dna_mutation_with_context(103, "TAC", "G");
        let mutated = DnaSeq::new("ACGGGT").unwrap();
        let interval =
            BaseInterval::new(BasePos::new(4).unwrap(), BasePos::new(4).unwrap()).unwrap();
        let (applied, alt_interval) = mutation.apply_mutation_with_alt_interval();

        assert_eq!(applied, mutated);
        assert_eq!(alt_interval, Some(interval.clone()));

        assert_eq!(
            mutation.format_mutated_sequence_and_highlight_mutated_bases(),
            mutated.format_with_coloured_interval(&interval)
        );
    }

    #[test]
    fn format_mutated_sequence_highlights_first_base_alt_allele() {
        let mutation = dna_mutation_with_context(100, "A", "G");
        let mutated = DnaSeq::new("GCGTACGT").unwrap();
        let interval =
            BaseInterval::new(BasePos::new(1).unwrap(), BasePos::new(1).unwrap()).unwrap();
        let (applied, alt_interval) = mutation.apply_mutation_with_alt_interval();

        assert_eq!(applied, mutated);
        assert_eq!(alt_interval, Some(interval.clone()));

        assert_eq!(
            mutation.format_mutated_sequence_and_highlight_mutated_bases(),
            mutated.format_with_coloured_interval(&interval)
        );
    }

    #[test]
    fn format_mutated_sequence_highlights_last_base_alt_allele() {
        let mutation = dna_mutation_with_context(107, "T", "A");
        let mutated = DnaSeq::new("ACGTACGA").unwrap();
        let interval =
            BaseInterval::new(BasePos::new(8).unwrap(), BasePos::new(8).unwrap()).unwrap();
        let (applied, alt_interval) = mutation.apply_mutation_with_alt_interval();

        assert_eq!(applied, mutated);
        assert_eq!(alt_interval, Some(interval.clone()));

        assert_eq!(
            mutation.format_mutated_sequence_and_highlight_mutated_bases(),
            mutated.format_with_coloured_interval(&interval)
        );
    }
}
