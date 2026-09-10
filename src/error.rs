use crate::{
    base::Alphabet,
    coords::{Interval1, Pos1},
};

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum BaseError {
    #[error("invalid {alphabet} base: '{invalid}'")]
    InvalidCharacter { alphabet: Alphabet, invalid: char },

    #[error("invalid {alphabet} byte value: {invalid}; expected an ASCII nucleotide symbol")]
    InvalidByte { alphabet: Alphabet, invalid: u8 },

    #[error("base '{base}' has ambiguous chemical class")]
    AmbiguousChemicalClass { base: char },
}

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum CoordError {
    #[error("position must be 1-based so 0 is not a valid position")]
    PositionIsZero,

    #[error("position {value} cannot be represented on this platform; max allowed position: {max}")]
    PositionOverflowU64 { value: u64, max: Pos1 },

    #[error("position {value} cannot be represented on this platform; max allowed position: {max}")]
    PositionOverflowU32 { value: u32, max: Pos1 },

    #[error("position underflow: {lhs} - {rhs} would be < {min}")]
    PositionUnderflow { lhs: usize, rhs: usize, min: usize },

    #[error("position overflow: {lhs} + {rhs} would exceed {max}")]
    PositionOverflowAdd { lhs: usize, rhs: usize, max: usize },

    #[error("start position of interval ({start}) must be <= end position ({end})")]
    InvalidIntervalCoords { start: usize, end: usize },

    #[error("end position of range ({end}) cannot be less than start position ({start})")]
    RangeEndTooSmall { start: usize, end: usize },
}

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum SequenceError {
    #[error(transparent)]
    InvalidBase(#[from] BaseError),

    #[error(
        "invalid subsequence coordinates: requested range [{start}, {end}) on a sequence of length {len}"
    )]
    InvalidSlice {
        start: usize,
        end: usize,
        len: usize,
    },

    #[error("cannot mutate interval: {interval}; it spans beyond sequence length {seqlength}")]
    FailedMutateInvalidInterval {
        interval: Interval1,
        seqlength: usize,
    },

    #[error("cannot determine palindrome for a sequence containing ambiguous bases")]
    AmbiguousPalindrome,

    #[error(
        "cannot convert degenerate sequence to concrete sequence: ambiguous base '{base}' at position {position}"
    )]
    CannotConvertDegenerateSequence { position: Pos1, base: char },
}

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum MutationError {
    #[error(transparent)]
    Coord(#[from] CoordError),

    #[error("cannot pyrimidine-center mutation because the reference sequence has no middle base")]
    MissingMiddleBase,

    #[error(
        "cannot pyrimidine-center mutation because middle base '{base}' has ambiguous chemical class"
    )]
    AmbiguousMiddleBase { base: char },

    #[error("strand field of mutation and context must both be either Some or None.")]
    MismatchedStrandOption,

    #[error("failed to apply mutation {mutation} to sequence context")]
    FailedToApplyMutationToContext {
        mutation: String,

        #[source]
        source: SequenceError,
    },

    #[error(
        "reference allele of mutation does NOT match context sequence at the expected position. Problematic mutation: {mutation}. Problematic context: {context}"
    )]
    MismatchedReferenceAlleleAndContextSeq { mutation: String, context: String },

    #[error(
        "name of the mutated chromosome ({mutated_chromosome}) is different to chromosome the sequence context originated from ({context_chromosome})"
    )]
    MismatchedChromosomeName {
        mutated_chromosome: String,
        context_chromosome: String,
    },
    #[error(
        "mutation position ({position}) is outside the interval {start}-{end} (1-based inclusive)"
    )]
    MutationPositionOutsideInterval {
        position: usize,
        start: usize,
        end: usize,
    },
}

// #[derive(thiserror::Error, Debug, PartialEq, Eq)]
// pub enum Error {
//     #[error(transparent)]
//     Base(#[from] BaseError),
//
//     #[error(transparent)]
//     Coord(#[from] CoordError),
//
//     #[error(transparent)]
//     Sequence(#[from] SequenceError),
//
//     #[error(transparent)]
//     Mutation(#[from] MutationError),
// }

// pub type Result<T> = std::result::Result<T, Error>;
