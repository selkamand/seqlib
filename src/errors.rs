use crate::base::Alphabet;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SeqError {
    #[error("sequence contains invalid characters for {alphabet:?}: {invalid}")]
    InvalidCharacters { alphabet: Alphabet, invalid: String },
    #[error("Invalid character for {alphabet:?}: {invalid}")]
    InvalidCharacter { alphabet: Alphabet, invalid: char },
    #[error("Invalid byte for {alphabet:?}: {invalid}")]
    InvalidByte { alphabet: Alphabet, invalid: u8 },
}
