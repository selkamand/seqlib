use crate::errors::SeqError;
use core::fmt;
use std::collections::BTreeSet;

/// A biological sequence with an associated alphabet.
///
/// `Seq` represents a validated nucleotide sequence (DNA or RNA).
/// The alphabet determines which characters are considered valid and
/// influences downstream operations (e.g. reverse complementing).
///
/// This struct is not appropriate for larger-than-memory sequences.
#[derive(Debug)]
pub struct Seq {
    seq: String,
    alphabet: Alphabet,
}

impl fmt::Display for Seq {
    /// Formats the sequence as its raw string representation.
    ///
    /// This prints only the underlying sequence characters, without
    /// additional metadata.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.seq)
    }
}

impl Seq {
    /// Returns the length of the sequence in characters.
    ///
    /// This corresponds to the number of bases in the sequence.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Returns `true` if the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Constructs a new `Seq` from a string slice and an alphabet.
    ///
    /// The sequence is validated against the provided alphabet. If any
    /// invalid characters are found, construction fails with a
    /// [`SeqError::InvalidCharacters`] error containing the unique set
    /// of offending characters.
    ///
    /// # Errors
    ///
    /// Returns `SeqError::InvalidCharacters` if the sequence contains
    /// characters not permitted by the given alphabet.
    pub fn new(seq: &str, alphabet: Alphabet) -> Result<Self, SeqError> {
        if seq.chars().any(|c| !alphabet.is_valid_char(c)) {
            let invalid: BTreeSet<char> = seq
                .chars()
                .filter(|&c| !alphabet.is_valid_char(c))
                .collect();

            if !invalid.is_empty() {
                let invalid = invalid.into_iter().collect::<String>();
                return Err(SeqError::InvalidCharacters { alphabet, invalid });
            }
        }

        Ok(Self {
            seq: seq.to_owned(),
            alphabet,
        })
    }

    /// Returns the length of the sequence in characters.
    ///
    /// This is an alias for [`Seq::len`].
    pub fn length(&self) -> usize {
        self.seq.len()
    }

    /// Returns a reference to the sequence’s alphabet.
    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }
}

/// The alphabet describing which nucleotide symbols are valid.
#[derive(Debug)]
pub enum Alphabet {
    /// DNA alphabet (A, C, G, T, N, and IUPAC ambiguity codes).
    DNA,
    /// RNA alphabet (A, C, G, U, N, and IUPAC ambiguity codes).
    RNA,
}

impl Alphabet {
    /// Returns `true` if the given character is valid for this alphabet.
    ///
    /// This method accepts both uppercase and lowercase characters and
    /// includes IUPAC ambiguity codes.
    fn is_valid_char(&self, c: char) -> bool {
        match self {
            Alphabet::DNA => matches!(
                c.to_ascii_uppercase(),
                'A' | 'C'
                    | 'G'
                    | 'T'
                    | 'R'
                    | 'Y'
                    | 'S'
                    | 'W'
                    | 'K'
                    | 'M'
                    | 'B'
                    | 'D'
                    | 'H'
                    | 'V'
                    | 'N'
            ),

            Alphabet::RNA => matches!(
                c.to_ascii_uppercase(),
                'A' | 'C'
                    | 'G'
                    | 'U'
                    | 'R'
                    | 'Y'
                    | 'S'
                    | 'W'
                    | 'K'
                    | 'M'
                    | 'B'
                    | 'D'
                    | 'H'
                    | 'V'
                    | 'N'
            ),
        }
    }
    fn complement_char(&self, c: char) -> Result<char, SeqError> {
        let is_lower = c.is_ascii_lowercase();
        let u = c.to_ascii_uppercase();
        let comp = match self {
            Alphabet::DNA => match u {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                'U' => 'A', // tolerate U in DNA input (optional; remove if you prefer strictness)
                'R' => 'Y',
                'Y' => 'R',
                'S' => 'S',
                'W' => 'W',
                'K' => 'M',
                'M' => 'K',
                'B' => 'V',
                'V' => 'B',
                'D' => 'H',
                'H' => 'D',
                'N' => 'N',
                _ => {
                    return Err(SeqError::InvalidCharacters {
                        alphabet: Alphabet::DNA,
                        invalid: u.to_string(),
                    });
                }
            },
            Alphabet::RNA => match u {
                'A' => 'U',
                'C' => 'G',
                'G' => 'C',
                'U' => 'A',
                'T' => 'A', // tolerate T in RNA input (optional; remove if you prefer strictness)
                'R' => 'Y',
                'Y' => 'R',
                'S' => 'S',
                'W' => 'W',
                'K' => 'M',
                'M' => 'K',
                'B' => 'V',
                'V' => 'B',
                'D' => 'H',
                'H' => 'D',
                'N' => 'N',
                _ => {
                    return Err(SeqError::InvalidCharacters {
                        alphabet: Alphabet::RNA,
                        invalid: u.to_string(),
                    });
                }
            },
        };

        Ok(if is_lower {
            comp.to_ascii_lowercase()
        } else {
            comp
        })
    }

    fn complement(&self, seq: &str) -> Result<String, SeqError> {
        seq.chars()
            .map(|c| self.complement_char(c)) // Iterator<Item = Result<char, SeqError>>
            .collect::<Result<String, SeqError>>() // early returns on first Err
    }
}
