// Data structures for representing and transforming individual nucleotides.
// We currently implement DnaBase and RnaBase, which both implement the Base trait.
// The Base trait defines the small set of operations any "nucleotide alphabet" must support.

use core::fmt;

use crate::errors::SeqError;

/// A minimal interface for “a single nucleotide symbol” (DNA, RNA, or a future alphabet).
///
/// Think of this trait as: *"What can I do with one base?"*
///
/// If you implement `Base` for an enum, you get:
/// - a **complement** mapping (A<->T/U, ambiguity codes too)
/// - conversion to **ASCII** bytes for printing / writing to files
/// - parsing from an ASCII byte with helpful errors
/// - a flag telling you if the base is ambiguous (e.g. R means A/G)
///
/// Why ASCII bytes (`u8`) instead of `char`?
/// - Real sequence files are bytes.
/// - DNA/RNA/IUPAC codes are ASCII characters.
/// - Using bytes is fast and avoids Unicode complexity.
///
/// The contract is:
/// - `complement` is infallible (it always returns something)
/// - `try_from_ascii` is the validation gate (it may fail)
pub trait Base: Copy + Eq + fmt::Debug + fmt::Display {
    /// Name of the alphabet
    const ALPHABET: Alphabet;

    /// Return the complement of this base.
    ///
    /// Examples (DNA):
    /// - A ↔ T
    /// - C ↔ G
    /// - R (A/G) ↔ Y (C/T)
    ///
    /// This must be infallible: every valid base must have a defined complement.
    fn complement(self) -> Self;

    /// Convert this base to an uppercase ASCII byte (e.g. `b'A'`).
    ///
    /// This is useful for:
    /// - printing
    /// - writing FASTA/FASTQ
    /// - building `String`s
    fn to_ascii(self) -> u8;

    fn to_char(self) -> char {
        self.to_ascii() as char
    }

    /// Convert this base to a lowercase ASCII byte (e.g. `b'a'`).
    ///
    /// Lowercase bases are often used for “soft-masking” (e.g. low-confidence regions),
    /// even though the biological base is the same.
    fn to_ascii_lower(self) -> u8;

    /// Check if base only represents one possible nucleotide (e.g. A/C/T/U/G).
    /// Bases like `R` can represent multiple possible nucleotides (in this case, A or G)
    fn is_unambiguous(self) -> bool {
        !self.is_ambiguous()
    }

    /// Parse a single ASCII byte into a base.
    ///
    /// This is the “gatekeeper” function: it checks whether an input symbol is allowed.
    ///
    /// Rules:
    /// - input must be ASCII (bytes 0–127)
    /// - parsing is case-insensitive (both `b'a'` and `b'A'` work)
    /// - invalid input returns a `SeqError` describing what went wrong
    fn try_from_ascii(b: u8) -> Result<Self, SeqError>;

    /// Returns `true` if this symbol can represent more than one concrete nucleotide.
    ///
    /// Examples:
    /// - `A` is not ambiguous (it always means A)
    /// - `R` is ambiguous (it means A **or** G)
    /// - `N` is ambiguous/unknown (it means any base)
    ///
    /// This is useful because some operations (like translation) usually require
    /// unambiguous sequences.
    fn is_ambiguous(self) -> bool;

    /// Classify the type of Base (Purine Vs Pyrimidine vs Uncertain)
    fn chemical_class(self) -> Option<ChemClass>;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum ChemClass {
    Purine,
    Pyrimidine,
}

/// DNA nucleotide symbols including IUPAC ambiguity codes.
///
/// Stored as an enum so Rust can enforce that functions like `complement` handle
/// every possible symbol (no missing cases).
#[repr(u8)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum DnaBase {
    A,
    C,
    G,
    T,
    /// Unknown / any base
    N,
    /// A or G
    R,
    /// C or T
    Y,
    /// G or C
    S,
    /// A or T
    W,
    /// G or T
    K,
    /// A or C
    M,
    /// C or G or T
    B,
    /// A or G or T
    D,
    /// A or C or T
    H,
    /// A or C or G
    V,
}

/// RNA nucleotide symbols including IUPAC ambiguity codes.
#[repr(u8)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum RnaBase {
    A,
    C,
    G,
    U,
    /// Unknown / any base
    N,
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
}

impl Base for DnaBase {
    const ALPHABET: Alphabet = Alphabet::DNA;

    fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::N => Self::N,
            Self::R => Self::Y,
            Self::Y => Self::R,
            Self::S => Self::S,
            Self::W => Self::W,
            Self::K => Self::M,
            Self::M => Self::K,
            Self::B => Self::V,
            Self::V => Self::B,
            Self::D => Self::H,
            Self::H => Self::D,
        }
    }

    fn try_from_ascii(b: u8) -> Result<Self, SeqError> {
        // If the byte is not ASCII, it cannot be a nucleotide symbol.
        if !b.is_ascii() {
            return Err(SeqError::InvalidByte {
                alphabet: Alphabet::DNA,
                invalid: b,
            });
        }

        match b.to_ascii_uppercase() {
            b'A' => Ok(DnaBase::A),
            b'C' => Ok(DnaBase::C),
            b'G' => Ok(DnaBase::G),
            b'T' => Ok(DnaBase::T),
            b'N' => Ok(DnaBase::N),
            b'R' => Ok(DnaBase::R),
            b'Y' => Ok(DnaBase::Y),
            b'S' => Ok(DnaBase::S),
            b'W' => Ok(DnaBase::W),
            b'K' => Ok(DnaBase::K),
            b'M' => Ok(DnaBase::M),
            b'B' => Ok(DnaBase::B),
            b'D' => Ok(DnaBase::D),
            b'H' => Ok(DnaBase::H),
            b'V' => Ok(DnaBase::V),
            _ => Err(SeqError::InvalidCharacter {
                alphabet: Alphabet::DNA,
                invalid: b as char,
            }),
        }
    }

    fn to_ascii(self) -> u8 {
        match self {
            DnaBase::A => b'A',
            DnaBase::C => b'C',
            DnaBase::G => b'G',
            DnaBase::T => b'T',
            DnaBase::N => b'N',
            DnaBase::R => b'R',
            DnaBase::Y => b'Y',
            DnaBase::S => b'S',
            DnaBase::W => b'W',
            DnaBase::K => b'K',
            DnaBase::M => b'M',
            DnaBase::B => b'B',
            DnaBase::D => b'D',
            DnaBase::H => b'H',
            DnaBase::V => b'V',
        }
    }

    fn to_ascii_lower(self) -> u8 {
        match self {
            DnaBase::A => b'a',
            DnaBase::C => b'c',
            DnaBase::G => b'g',
            DnaBase::T => b't',
            DnaBase::N => b'n',
            DnaBase::R => b'r',
            DnaBase::Y => b'y',
            DnaBase::S => b's',
            DnaBase::W => b'w',
            DnaBase::K => b'k',
            DnaBase::M => b'm',
            DnaBase::B => b'b',
            DnaBase::D => b'd',
            DnaBase::H => b'h',
            DnaBase::V => b'v',
        }
    }

    fn is_ambiguous(self) -> bool {
        !matches!(self, DnaBase::A | DnaBase::C | DnaBase::G | DnaBase::T)
    }

    fn chemical_class(self) -> Option<ChemClass> {
        match self {
            DnaBase::A => Some(ChemClass::Purine),
            DnaBase::G => Some(ChemClass::Purine),
            DnaBase::C => Some(ChemClass::Pyrimidine),
            DnaBase::T => Some(ChemClass::Pyrimidine),
            DnaBase::N => None,
            DnaBase::R => Some(ChemClass::Purine),
            DnaBase::Y => Some(ChemClass::Pyrimidine),
            DnaBase::S => None,
            DnaBase::W => None,
            DnaBase::K => None,
            DnaBase::M => None,
            DnaBase::B => None,
            DnaBase::D => None,
            DnaBase::H => None,
            DnaBase::V => None,
        }
    }
}

impl fmt::Display for DnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use core::fmt::Write;
        f.write_char((*self).to_char())
    }
}

impl Base for RnaBase {
    const ALPHABET: Alphabet = Alphabet::RNA;

    fn try_from_ascii(b: u8) -> Result<Self, SeqError> {
        if !b.is_ascii() {
            return Err(SeqError::InvalidByte {
                alphabet: Alphabet::RNA,
                invalid: b,
            });
        }

        match b.to_ascii_uppercase() {
            b'A' => Ok(RnaBase::A),
            b'C' => Ok(RnaBase::C),
            b'G' => Ok(RnaBase::G),
            b'U' => Ok(RnaBase::U),
            b'N' => Ok(RnaBase::N),
            b'R' => Ok(RnaBase::R),
            b'Y' => Ok(RnaBase::Y),
            b'S' => Ok(RnaBase::S),
            b'W' => Ok(RnaBase::W),
            b'K' => Ok(RnaBase::K),
            b'M' => Ok(RnaBase::M),
            b'B' => Ok(RnaBase::B),
            b'D' => Ok(RnaBase::D),
            b'H' => Ok(RnaBase::H),
            b'V' => Ok(RnaBase::V),
            _ => Err(SeqError::InvalidCharacter {
                alphabet: Alphabet::RNA,
                invalid: b as char,
            }),
        }
    }

    fn complement(self) -> Self {
        match self {
            Self::A => Self::U,
            Self::U => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::R => Self::Y,
            Self::Y => Self::R,
            Self::S => Self::S,
            Self::W => Self::W,
            Self::K => Self::M,
            Self::M => Self::K,
            Self::B => Self::V,
            Self::V => Self::B,
            Self::D => Self::H,
            Self::H => Self::D,
            Self::N => Self::N,
        }
    }

    fn to_ascii(self) -> u8 {
        match self {
            RnaBase::A => b'A',
            RnaBase::C => b'C',
            RnaBase::G => b'G',
            RnaBase::U => b'U',
            RnaBase::N => b'N',
            RnaBase::R => b'R',
            RnaBase::Y => b'Y',
            RnaBase::S => b'S',
            RnaBase::W => b'W',
            RnaBase::K => b'K',
            RnaBase::M => b'M',
            RnaBase::B => b'B',
            RnaBase::D => b'D',
            RnaBase::H => b'H',
            RnaBase::V => b'V',
        }
    }

    fn to_ascii_lower(self) -> u8 {
        match self {
            RnaBase::A => b'a',
            RnaBase::C => b'c',
            RnaBase::G => b'g',
            RnaBase::U => b'u',
            RnaBase::N => b'n',
            RnaBase::R => b'r',
            RnaBase::Y => b'y',
            RnaBase::S => b's',
            RnaBase::W => b'w',
            RnaBase::K => b'k',
            RnaBase::M => b'm',
            RnaBase::B => b'b',
            RnaBase::D => b'd',
            RnaBase::H => b'h',
            RnaBase::V => b'v',
        }
    }

    fn is_ambiguous(self) -> bool {
        !matches!(self, RnaBase::A | RnaBase::C | RnaBase::G | RnaBase::U)
    }

    fn chemical_class(self) -> Option<ChemClass> {
        match self {
            RnaBase::A => Some(ChemClass::Purine),
            RnaBase::G => Some(ChemClass::Purine),
            RnaBase::C => Some(ChemClass::Pyrimidine),
            RnaBase::U => Some(ChemClass::Pyrimidine),
            RnaBase::N => None,
            RnaBase::R => Some(ChemClass::Purine),
            RnaBase::Y => Some(ChemClass::Pyrimidine),
            RnaBase::S => None,
            RnaBase::W => None,
            RnaBase::K => None,
            RnaBase::M => None,
            RnaBase::B => None,
            RnaBase::D => None,
            RnaBase::H => None,
            RnaBase::V => None,
        }
    }
}

impl fmt::Display for RnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use core::fmt::Write;
        f.write_char((*self).to_char())
    }
}

/// Names of supported alphabets.
///
/// This is mainly used for error reporting (e.g. “invalid character for DNA”).
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Alphabet {
    DNA,
    RNA,
}

impl fmt::Display for Alphabet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Alphabet::DNA => write!(f, "DNA"),
            Alphabet::RNA => write!(f, "RNA"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_try_from_ascii_accepts_case_insensitive() {
        assert_eq!(DnaBase::try_from_ascii(b'a').unwrap(), DnaBase::A);
        assert_eq!(DnaBase::try_from_ascii(b'A').unwrap(), DnaBase::A);
        assert_eq!(DnaBase::try_from_ascii(b't').unwrap(), DnaBase::T);
    }

    #[test]
    fn rna_try_from_ascii_accepts_case_insensitive() {
        assert_eq!(RnaBase::try_from_ascii(b'u').unwrap(), RnaBase::U);
        assert_eq!(RnaBase::try_from_ascii(b'U').unwrap(), RnaBase::U);
    }

    #[test]
    fn complement_is_involutive_for_some_representative_bases() {
        // "Involutive" means comp(comp(x)) == x.
        // We only test a few representative bases to keep this minimal.
        let reps_dna = [DnaBase::A, DnaBase::C, DnaBase::R, DnaBase::B, DnaBase::N];
        for b in reps_dna {
            assert_eq!(b.complement().complement(), b);
        }

        let reps_rna = [RnaBase::A, RnaBase::C, RnaBase::R, RnaBase::B, RnaBase::N];
        for b in reps_rna {
            assert_eq!(b.complement().complement(), b);
        }
    }

    #[test]
    fn ascii_rendering_is_consistent() {
        // Upper then lower should match ASCII casing expectations.
        let b = DnaBase::G;
        assert_eq!(b.to_ascii(), b'G');
        assert_eq!(b.to_ascii_lower(), b'g');

        let r = RnaBase::U;
        assert_eq!(r.to_ascii(), b'U');
        assert_eq!(r.to_ascii_lower(), b'u');
    }

    #[test]
    fn ambiguity_flag_matches_expectations() {
        assert!(!DnaBase::A.is_ambiguous());
        assert!(DnaBase::N.is_ambiguous());
        assert!(DnaBase::R.is_ambiguous());

        assert!(!RnaBase::U.is_ambiguous());
        assert!(RnaBase::N.is_ambiguous());
    }
}
