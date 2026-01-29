use crate::base::Alphabet;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error(
        "Invalid {alphabet} base: '{invalid}'. \
         Allowed symbols are standard bases plus IUPAC ambiguity codes for {alphabet}."
    )]
    InvalidCharacter { alphabet: Alphabet, invalid: char },

    #[error(
        "Invalid {alphabet} byte value: {invalid}. \
         Expected an ASCII letter representing a nucleotide (e.g. A,C,G,T/U,N)."
    )]
    InvalidByte { alphabet: Alphabet, invalid: u8 },

    #[error(
        "Invalid subsequence coordinates: requested range [{start}, {end}) on a sequence of length {len}. \
         Indices are 0-based and the end position is exclusive (like Rust slicing)."
    )]
    InvalidSlice {
        start: usize,
        end: usize,
        len: usize,
    },

    #[error("Position must be 1-based so 0 is not a valid position")]
    PositionIsZero,
}

pub type Result<T> = std::result::Result<T, Error>;
