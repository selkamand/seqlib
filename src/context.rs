use crate::error::{Error, Result};
use crate::{base::Base, sequence::Seq};
use std::num::NonZeroUsize;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContextWindow<B: Base> {
    seq: Seq<B>,
    start: u64,
    anchor: u64, // Position in self.seq about which the context was derived
    orientation: Orientation,
}

impl<B: Base> std::fmt::Display for ContextWindow<B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let end = if self.seq.is_empty() {
            self.start
        } else {
            self.start + self.seq.len() as u64 - 1
        };

        write!(
            f,
            "Sequence: {} | from {}-{} (both 1-based) | Orientation: {}",
            self.seq, self.start, end, self.orientation
        )
    }
}

impl<B: Base> ContextWindow<B> {
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Orientation {
    Forward, // as provided by the reference sequence in this coord space
    Reverse, // if you ever store reverse-oriented windows (I’d avoid this if possible)
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Orientation::Forward => "Forward",
            Orientation::Reverse => "Reverse",
        };
        write!(f, "{s}")
    }
}

pub struct Pos(core::num::NonZeroUsize);

impl std::fmt::Display for Pos {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Pos {
    pub fn new(position: usize) -> Result<Self> {
        if let Some(m) = NonZeroUsize::new(position) {
            Ok(Self(m))
        } else {
            Err(Error::PositionIsZero)
        }
    }
}
