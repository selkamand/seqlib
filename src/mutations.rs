use core::fmt;

use crate::{
    base::{Base, ChemClass},
    sequences::{DnaSeq, Seq},
};

// Small Variant Class
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmallMutation {
    chromosome: String,
    position: u64, // 1-based position (start)
    reference: DnaSeq,
    alternative: DnaSeq,
    multiallelic: bool,
    context: DnaSeq,
}

// Implement the `fmt::Display` trait for `Point`.
impl fmt::Display for SmallMutation {
    // This trait requires the `fmt` function with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use the write! macro to format the output.
        write!(
            f,
            "{}:{} {}>{} (delta: {}; class: {}; multiallelic:{})",
            self.chromosome,
            self.position,
            self.reference,
            self.alternative,
            self.delta(),
            self.class(),
            self.multiallelic
        )
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
    pub fn from_chemical_bases(reference: ChemClass, alternative: ChemClass) -> TiTv {
        match reference == alternative {
            true => TiTv::Transition,
            false => TiTv::Transversion,
        }
    }
}

impl SmallMutation {
    pub fn reflen(&self) -> usize {
        self.reference.len()
    }
    pub fn altlen(&self) -> usize {
        self.alternative.len()
    }

    pub fn delta(&self) -> i64 {
        self.altlen() as i64 - self.reflen() as i64
    }

    pub fn class(&self) -> SmallMutationType {
        SmallMutationType::from_lengths(self.reflen(), self.altlen())
    }

    pub fn class_titv(&self) -> Option<TiTv> {
        match self.class() {
            SmallMutationType::SNV => Some(TiTv::from_chemical_bases(
                self.reference.middlebase()?.chemical_class(),
                self.alternative.middlebase()?.chemical_class(),
            )),
            SmallMutationType::DOUBLET => None,
            SmallMutationType::MNV => None,
            SmallMutationType::INSERTION => None,
            SmallMutationType::DELETION => None,
        }
    }
}
