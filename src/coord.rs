use crate::error::{Error, Result};
use std::num::NonZeroUsize;

/// A non-zero 1-based coordinate.
///
/// This type is intended for biological coordinate systems that are conventionally
/// 1-based (e.g. VCF POS). It prevents accidental construction of an invalid `0`
/// coordinate, which helps avoid off-by-one bugs when converting to 0-based indices
/// for slicing.
///
/// # Invariants
/// - Always `>= 1`.
///
/// # Notes
/// - This type stores a `NonZeroUsize`. On both 32 and 64-bit platforms this comfortably fits common
///   genome/transcript coordinate ranges (32 bits: 4,294,967,295, 64 bits: 18,446,744,073,709,551,615).
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Pos(NonZeroUsize);

impl Pos {
    /// Maximum Allowed Position Value
    pub const MAX: Self = Self(NonZeroUsize::MAX);

    /// Minimum Allowed Position Value
    pub const MIN: Self = Self(NonZeroUsize::MIN);

    /// Construct a 1-based position.
    ///
    /// # Errors
    /// Returns [`Error::PositionIsZero`] if `position == 0`.
    pub fn new(position: usize) -> Result<Self> {
        NonZeroUsize::new(position)
            .map(Pos)
            .ok_or(Error::PositionIsZero)
    }

    /// Return the underlying 1-based coordinate as a `usize`.
    pub fn get(self) -> usize {
        self.0.get()
    }

    // Position Shifting

    /// Add an offset to this position.
    ///
    /// Returns `None` if the result would overflow `usize`.
    pub fn checked_add(self, offset: usize) -> Option<Self> {
        let v = self.get().checked_add(offset)?;
        Pos::new(v).ok()
    }

    /// Add an offset, saturating at `Pos::MAX` on overflow.
    pub fn saturating_add(self, offset: usize) -> Self {
        let v = self.get().saturating_add(offset);
        // `v` is never 0 here, so `new` cannot fail.
        // But we still avoid unwrap by falling back to MAX defensively.
        Pos::new(v).unwrap_or(Pos::MAX)
    }

    /// Subtract an offset, saturating at `Pos::MIN` (Position 1) on underflow.
    pub fn saturating_sub(self, offset: usize) -> Self {
        let v = self.get().saturating_sub(offset);
        Pos::new(v).unwrap_or(Pos::MIN)
    }

    /// Add an offset to this position.
    ///
    /// # Errors
    /// Returns [`Error::PositionOverflowAdd`] if `self + offset` cannot be represented
    /// on this platform.
    pub fn try_add(self, offset: usize) -> Result<Self> {
        match self.get().checked_add(offset) {
            Some(v) => Pos::new(v).map_err(|_| Error::PositionOverflowAdd {
                lhs: self,
                rhs: offset,
                max: Pos::MAX,
            }),
            None => Err(Error::PositionOverflowAdd {
                lhs: self,
                rhs: offset,
                max: Pos::MAX,
            }),
        }
    }

    /// Subtract an offset from this position.
    ///
    /// # Errors
    /// Returns [`Error::PositionUnderflow`] if `self - offset` would be < 1.
    pub fn try_sub(self, offset: usize) -> Result<Self> {
        match self.get().checked_sub(offset) {
            Some(v) => Pos::new(v).map_err(|_| Error::PositionUnderflow {
                lhs: self,
                rhs: offset,
            }),
            None => Err(Error::PositionUnderflow {
                lhs: self,
                rhs: offset,
            }),
        }
    }
}

impl core::fmt::Display for Pos {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.get())
    }
}
impl Default for Pos {
    fn default() -> Self {
        Pos::MIN
    }
}

impl TryFrom<u64> for Pos {
    type Error = Error;

    /// Fallibly convert a `u64` into a 1-based [`Pos`].
    ///
    /// This conversion is **platform dependent** because [`Pos`] stores a `NonZeroUsize`.
    /// On targets where `usize` is smaller than `u64` (e.g. 32-bit or 16-bit), large
    /// values may not be representable and will be rejected.
    ///
    /// # Errors
    /// - [`Error::PositionIsZero`] if `value == 0`.
    /// - [`Error::PositionOverflowU64`] if `value` cannot be represented as a `usize`
    ///   on the current platform.
    fn try_from(value: u64) -> Result<Self> {
        if value == 0 {
            return Err(Error::PositionIsZero);
        }

        // Fail on 32-bit (or any platform) if it doesn't fit in usize.
        let as_usize = usize::try_from(value).map_err(|_| Error::PositionOverflowU64 {
            value,
            max: Pos::MAX,
        })?;

        // as_usize is non-zero because value != 0
        Pos::new(as_usize)
    }
}

impl TryFrom<u32> for Pos {
    type Error = Error;

    /// Fallibly convert a `u32` into a 1-based [`Pos`].
    ///
    /// This conversion is always safe on 32-bit and 64-bit targets, but may fail on
    /// narrower targets (e.g. 16-bit) where `usize::MAX < u32::MAX`.
    ///
    /// # Errors
    /// - [`Error::PositionIsZero`] if `value == 0`.
    /// - [`Error::PositionOverflowU32`] if `value` cannot be represented as a `usize`
    ///   on the current platform.
    fn try_from(value: u32) -> Result<Self> {
        if value == 0 {
            return Err(Error::PositionIsZero);
        }

        let as_usize = usize::try_from(value).map_err(|_| Error::PositionOverflowU32 {
            value,
            max: Pos::MAX,
        })?;

        Pos::new(as_usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::convert::TryFrom;

    #[test]
    fn new_rejects_zero() {
        let err = Pos::new(0).unwrap_err();
        assert_eq!(err, Error::PositionIsZero);
    }

    #[test]
    fn new_accepts_one_and_get_roundtrips() {
        let p = Pos::new(1).unwrap();
        assert_eq!(p.get(), 1);

        let p2 = Pos::new(42).unwrap();
        assert_eq!(p2.get(), 42);
    }

    #[test]
    fn min_and_max_constants_are_sane() {
        assert_eq!(Pos::MIN.get(), 1);
        assert_eq!(Pos::MAX.get(), usize::MAX);
        assert!(Pos::MAX.get() >= Pos::MIN.get());
    }

    #[test]
    fn display_prints_numeric_value() {
        let p = Pos::new(123).unwrap();
        assert_eq!(p.to_string(), "123");
    }

    #[test]
    fn try_from_u32_rejects_zero() {
        let err = Pos::try_from(0_u32).unwrap_err();
        assert_eq!(err, Error::PositionIsZero);
    }

    #[test]
    fn try_from_u64_rejects_zero() {
        let err = Pos::try_from(0_u64).unwrap_err();
        assert_eq!(err, Error::PositionIsZero);
    }

    #[test]
    fn try_from_u32_accepts_nonzero() {
        let p = Pos::try_from(1_u32).unwrap();
        assert_eq!(p.get(), 1);

        let p2 = Pos::try_from(123_u32).unwrap();
        assert_eq!(p2.get(), 123);
    }

    #[test]
    fn try_from_u64_accepts_nonzero_that_fits() {
        let p = Pos::try_from(1_u64).unwrap();
        assert_eq!(p.get(), 1);

        let p2 = Pos::try_from(123_u64).unwrap();
        assert_eq!(p2.get(), 123);
    }

    // --- overflow behavior depends on pointer width ---

    #[cfg(target_pointer_width = "16")]
    #[test]
    fn try_from_u32_rejects_values_that_do_not_fit_on_16bit() {
        // On 16-bit: usize::MAX is 65535, so 65536 should overflow.
        let v: u32 = (u16::MAX as u32) + 1;

        let err = Pos::try_from(v).unwrap_err();
        match err {
            Error::PositionOverflowU32 { value, max } => {
                assert_eq!(value, v);
                assert_eq!(max, Pos::MAX);
            }
            other => panic!("expected PositionOverflowU32, got {other:?}"),
        }
    }

    #[cfg(any(target_pointer_width = "16", target_pointer_width = "32"))]
    #[test]
    fn try_from_u64_rejects_values_that_do_not_fit_on_non_64bit() {
        // On 16/32-bit: pick a value > usize::MAX.
        let v: u64 = (usize::MAX as u64) + 1;

        let err = Pos::try_from(v).unwrap_err();
        match err {
            Error::PositionOverflowU64 { value, max } => {
                assert_eq!(value, v);
                assert_eq!(max, Pos::MAX);
            }
            other => panic!("expected PositionOverflowU64, got {other:?}"),
        }
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn try_from_u64_accepts_large_values_on_64bit() {
        // On 64-bit platforms, any non-zero u64 should fit into usize? Not quite:
        // usize::MAX == u64::MAX on 64-bit, so yes, all non-zero u64 fit.
        let v: u64 = u64::MAX;
        let p = Pos::try_from(v).unwrap();
        assert_eq!(p.get() as u64, v);
    }
}

#[cfg(test)]
mod pos_arith_tests {
    use super::*;

    #[test]
    fn try_add_ok() {
        let p = Pos::new(10).unwrap();
        let q = p.try_add(5).unwrap();
        assert_eq!(q.get(), 15);
    }

    #[test]
    fn try_sub_ok() {
        let p = Pos::new(10).unwrap();
        let q = p.try_sub(3).unwrap();
        assert_eq!(q.get(), 7);
    }

    #[test]
    fn try_sub_underflow_to_zero_errors() {
        let p = Pos::new(1).unwrap();
        let err = p.try_sub(1).unwrap_err();
        match err {
            Error::PositionUnderflow { lhs, rhs } => {
                assert_eq!(lhs, p);
                assert_eq!(rhs, 1);
            }
            other => panic!("expected PositionUnderflow, got {other:?}"),
        }
    }

    #[test]
    fn try_sub_underflow_below_zero_errors() {
        let p = Pos::new(1).unwrap();
        let err = p.try_sub(2).unwrap_err();
        match err {
            Error::PositionUnderflow { lhs, rhs } => {
                assert_eq!(lhs, p);
                assert_eq!(rhs, 2);
            }
            other => panic!("expected PositionUnderflow, got {other:?}"),
        }
    }

    #[test]
    fn try_add_overflow_errors() {
        let p = Pos::MAX;
        let err = p.try_add(1).unwrap_err();
        match err {
            Error::PositionOverflowAdd { lhs, rhs, max } => {
                assert_eq!(lhs, p);
                assert_eq!(rhs, 1);
                assert_eq!(max, Pos::MAX);
            }
            other => panic!("expected PositionOverflowAdd, got {other:?}"),
        }
    }
}
