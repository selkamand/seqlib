//! Definition of a Position (Pos) struct. Is 1-based and non-zero (enforced by typesystem).
use const_panic::concat_panic;

use crate::error::CoordError as Error;
pub(crate) type Result<T> = std::result::Result<T, Error>;

use std::num::NonZeroUsize;

/// A position in a zero-base inter-base coordinate system
///
/// Numbers are assigned to the space between bases (starting at 0).
///
/// ```{text}
///   A   T   A   C   G
/// 0   1   2   3   4   5
/// ```
///
/// So an [`InterbasePos`] of `4` doesn't mean much by itself, but packaged into an [`InterbaseInterval`] (e.g. `1-4`) unambiguously describe the sequence (`TAC`).
///
/// Interbase coordinate systems are also great for unambiguosly describing mutated sequences (including insertions) (which happen between bases).
/// This is why (inspired by the GA4GH Variant Representation Specification) `seqlib` mutation data types use interbase coordinates.
///
/// # Example
/// ```
/// use seqlib::coords::InterbasePos;
///
/// // Define first position
/// let position = InterbasePos::from(0usize);
///
/// // Define second position
/// let position = InterbasePos::from(1usize);
/// ```
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InterbasePos(usize);

impl InterbasePos {
    /// Maximum Allowed Position Value
    pub const MAX: Self = Self(usize::MAX);

    /// Minimum Allowed Position Value
    pub const MIN: Self = Self(usize::MIN);

    /// Get position as a usize
    pub fn get(&self) -> usize {
        self.0
    }

    /// Create an [`InterbasePos`] from a usize (infallable as any valid usize is a valid position)
    pub const fn new(position: usize) -> Self {
        Self(position)
    }

    // <- Position Shifting ->

    /// Add an offset to this position.
    ///
    /// Returns `None` if the result would overflow `usize`.
    pub fn checked_add(self, offset: usize) -> Option<Self> {
        let v = self.get().checked_add(offset)?;
        Some(Self::from(v))
    }

    /// Add an offset, saturating at `InterbasePos::MAX` on overflow.
    pub fn saturating_add(self, offset: usize) -> Self {
        // Pos::MAX is usize::MAX so we can just use usize saturating_add
        Self::from(self.get().saturating_add(offset))
    }

    /// Subtract an offset, saturating at `InterbasePos::MIN` on underflow.
    pub fn saturating_sub(self, offset: usize) -> Self {
        Self::from(self.get().saturating_sub(offset))
    }

    /// Add an offset to this position.
    ///
    /// # Errors
    /// Returns [`Error::PositionOverflowAdd`] if `self + offset` cannot be represented
    /// on this platform.
    pub fn try_add(self, offset: usize) -> Result<Self> {
        match self.get().checked_add(offset) {
            Some(p) => Ok(Self::from(p)),
            None => Err(Error::PositionOverflowAdd {
                lhs: self.into(),
                rhs: offset,
                max: Self::MAX.into(),
            }),
        }
    }

    /// Subtract an offset from this position.
    ///
    /// # Errors
    /// Returns [`Error::PositionUnderflow`] if `self - offset` would be < 1.
    pub fn try_sub(self, offset: usize) -> Result<Self> {
        match self.get().checked_sub(offset) {
            Some(p) => Ok(InterbasePos::from(p)),
            None => Err(Error::PositionUnderflow {
                lhs: self.into(),
                rhs: offset,
                min: Self::MIN.into(),
            }),
        }
    }
}

// Conversions into InterbasePos
impl From<usize> for InterbasePos {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<BasePos> for InterbasePos {
    fn from(value: BasePos) -> Self {
        let x = value.get() - 1; // Can never fail (go below zero) because BasePos is a NonZeroUsize
        Self(x)
    }
}

// Converions out of InterbasePos
impl From<InterbasePos> for usize {
    fn from(value: InterbasePos) -> Self {
        value.get()
    }
}

impl std::fmt::Display for InterbasePos {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.get())
    }
}

/// A position in a 1-based in-base coordinate system.
///
/// This type is intended for biological coordinate systems that are conventionally
/// 1-based (e.g. VCF POS). It prevents accidental construction of an invalid `0`
/// coordinate, which helps avoid off-by-one bugs when converting to 0-based indices
/// for slicing.
///
///
///
/// Each base is numbered starting at 1
///
/// ```{text}
/// A T A C G
/// 1 2 3 4 5
/// ```
///
/// So a [`BasePos`] of `4` refers to a `C`
/// and an [`BaseInterval`] of 2-4 refers to the 3bp sequence `TAC` (both-end inclusive)
///
///
/// # Invariants
/// - Always `>= 1`.
///
/// # Notes
/// - This type stores a `NonZeroUsize`. On both 32 and 64-bit platforms this comfortably fits common
///   genome/transcript coordinate ranges (32 bits: 4,294,967,295, 64 bits: 18,446,744,073,709,551,615).
#[repr(transparent)]
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BasePos(NonZeroUsize);

impl BasePos {
    /// Maximum Allowed Position Value
    pub const MAX: Self = Self(NonZeroUsize::MAX);

    /// Minimum Allowed Position Value
    pub const MIN: Self = Self(NonZeroUsize::MIN);

    /// Construct a 1-based [`BasePos`].
    ///
    /// # Errors
    /// Returns [`Error::PositionIsZero`] if `position == 0`.
    pub fn new(position: usize) -> Result<Self> {
        match NonZeroUsize::new(position) {
            Some(validpos) => Ok(Self(validpos)),
            None => Err(Error::PositionIsZero),
        }
    }

    /// Create a Pos - compile time panic if it fails. Powers the pos! macro
    pub const fn new_panic(position: usize) -> Self {
        match NonZeroUsize::new(position) {
            Some(validpos) => Self(validpos),
            None => concat_panic!(
                "Failed to create position from [",
                position,
                "]. Must be a usize from ",
                NonZeroUsize::MIN,
                "-",
                NonZeroUsize::MAX
            ),
        }
    }

    /// An unchecked constructor that works because all NonZeroUsize values are valid positions.
    /// Powers the pos! macro
    pub const fn new_unchecked(position: NonZeroUsize) -> Self {
        BasePos(position)
    }
    /// Return the underlying 1-based coordinate as a `usize`.
    pub fn get(self) -> usize {
        self.0.get()
    }

    /// Return the position as a 0based index (e.g. for indexnig into a `Seq` object)
    pub fn as_0based_index(&self) -> usize {
        self.get().saturating_sub(1)
    }

    // <- Position Shifting ->

    /// Add an offset to this position.
    ///
    /// Returns `None` if the result would overflow `usize`.
    pub fn checked_add(self, offset: usize) -> Option<Self> {
        let v = self.get().checked_add(offset)?;
        BasePos::new(v).ok()
    }

    /// Add an offset, saturating at `Pos::MAX` on overflow.
    pub fn saturating_add(self, offset: usize) -> Self {
        let v = self.get().saturating_add(offset);
        // `v` is never 0 here, so `new` cannot fail.
        // But we still avoid unwrap by falling back to MAX defensively.
        BasePos::new(v).unwrap_or(BasePos::MAX)
    }

    /// Subtract an offset, saturating at `Pos::MIN` (Position 1) on underflow.
    pub fn saturating_sub(self, offset: usize) -> Self {
        let v = self.get().saturating_sub(offset);
        BasePos::new(v).unwrap_or(BasePos::MIN)
    }

    /// Add an offset to this position.
    ///
    /// # Errors
    /// Returns [`Error::PositionOverflowAdd`] if `self + offset` cannot be represented
    /// on this platform.
    pub fn try_add(self, offset: usize) -> Result<Self> {
        match self.get().checked_add(offset) {
            Some(v) => BasePos::new(v).map_err(|_| Error::PositionOverflowAdd {
                lhs: self.into(),
                rhs: offset,
                max: Self::MAX.into(),
            }),
            None => Err(Error::PositionOverflowAdd {
                lhs: self.into(),
                rhs: offset,
                max: Self::MAX.into(),
            }),
        }
    }

    /// Subtract an offset from this position.
    ///
    /// # Errors
    /// Returns [`Error::PositionUnderflow`] if `self - offset` would be < 1.
    pub fn try_sub(self, offset: usize) -> Result<Self> {
        match self.get().checked_sub(offset) {
            Some(v) => BasePos::new(v).map_err(|_| Error::PositionUnderflow {
                lhs: self.into(),
                rhs: offset,
                min: Self::MIN.into(),
            }),
            None => Err(Error::PositionUnderflow {
                lhs: self.into(),
                rhs: offset,
                min: Self::MIN.into(),
            }),
        }
    }

    /// Construct a [`Pos`] from a [`NonZeroUsize`]
    ///
    /// This is infallable because every `NonZeroUsize` is a valid 1-based position
    pub fn from_nonzero(position: NonZeroUsize) -> Self {
        Self(position)
    }

    /// Return the underlying non-zero 1-based coordinate.
    pub const fn as_nonzero(self) -> NonZeroUsize {
        self.0
    }
}

impl core::fmt::Display for BasePos {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.get())
    }
}
impl Default for BasePos {
    fn default() -> Self {
        BasePos::MIN
    }
}

impl TryFrom<u64> for BasePos {
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
            max: BasePos::MAX,
        })?;

        // as_usize is non-zero because value != 0
        BasePos::new(as_usize)
    }
}

impl TryFrom<u32> for BasePos {
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
            max: BasePos::MAX,
        })?;

        BasePos::new(as_usize)
    }
}

impl From<NonZeroUsize> for BasePos {
    fn from(value: NonZeroUsize) -> Self {
        Self::from_nonzero(value)
    }
}

impl From<BasePos> for NonZeroUsize {
    fn from(value: BasePos) -> Self {
        value.0
    }
}

impl From<BasePos> for usize {
    fn from(value: BasePos) -> Self {
        value.get()
    }
}
/// Construct a [`BasePos`] from a **compile-time** integer literal.
///
/// This macro is intended for constant contexts and test code where the position
/// is known at compile time.
///
/// - `basepos!(1)` expands to a [`BasePos`] representing 1.
/// - `basepos!(0)` is rejected (a [`BasePos`] is always >= 1).
///
/// # Failure mode
/// This macro does **not** introduce a runtime panic in normal use:
/// it evaluates the invariant at compile time for literal inputs.
/// If the value is invalid (e.g. `0`), compilation fails.
///
/// For dynamic values (runtime variables), use [`BasePos::new`] instead.
///
/// # Examples
/// ```
/// use seqlib::coords::{BasePos};
/// use seqlib::basepos;
///
/// const P: BasePos = basepos!(123);
/// assert_eq!(P.get(), 123);
/// ```
#[macro_export]
macro_rules! basepos {
    ($lit:literal) => {{
        const P: BasePos = BasePos::new_panic($lit);
        P

        // Safe because 0 is handled above.
        // const P: core::num::NonZeroUsize = match core::num::NonZeroUsize::new($lit) {
        // Some(v) => v,
        // None => compile_error!("Cannot create valid position from ($lit)"),
        // $crate::coords::Pos::new_unchecked(P)
    }};
}

/// Construct an [`InterbasePos`] from a **compile-time** integer literal.
///
/// This macro is intended for constant contexts and test code where the position
/// is known at compile time.
///
/// - `interbasepos!(0)` expands to an [`InterbasePos`] representing 0.
///
/// # Failure mode
/// This macro does **not** introduce a runtime panic in normal use:
/// it evaluates the invariant at compile time for literal inputs.
///
/// For dynamic values (runtime variables), use [`InterbasePos::from`] or [`InterbasePos::new`] instead.
///
/// # Examples
/// ```
/// use seqlib::coords::{BasePos};
/// use seqlib::interbasepos;
///
/// const P: InterbasePos = interbasepos!(123);
/// assert_eq!(P.get(), 123);
/// ```
#[macro_export]
macro_rules! interbasepos {
    ($lit:literal) => {{
        const P: InterbasePos = InterbasePos::new_panic($lit);
        P
    }};
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::convert::TryFrom;

    #[test]
    fn new_rejects_zero() {
        let err = BasePos::new(0).unwrap_err();
        assert_eq!(err, Error::PositionIsZero);
    }

    #[test]
    fn new_accepts_one_and_get_roundtrips() {
        let p = BasePos::new(1).unwrap();
        assert_eq!(p.get(), 1);

        let p2 = BasePos::new(42).unwrap();
        assert_eq!(p2.get(), 42);
    }

    #[test]
    fn min_and_max_constants_are_sane() {
        assert_eq!(BasePos::MIN.get(), 1);
        assert_eq!(BasePos::MAX.get(), usize::MAX);
        assert!(BasePos::MAX.get() >= BasePos::MIN.get());
    }

    #[test]
    fn display_prints_numeric_value() {
        let p = BasePos::new(123).unwrap();
        assert_eq!(p.to_string(), "123");
    }

    #[test]
    fn try_from_u32_rejects_zero() {
        let err = BasePos::try_from(0_u32).unwrap_err();
        assert_eq!(err, Error::PositionIsZero);
    }

    #[test]
    fn try_from_u64_rejects_zero() {
        let err = BasePos::try_from(0_u64).unwrap_err();
        assert_eq!(err, Error::PositionIsZero);
    }

    #[test]
    fn try_from_u32_accepts_nonzero() {
        let p = BasePos::try_from(1_u32).unwrap();
        assert_eq!(p.get(), 1);

        let p2 = BasePos::try_from(123_u32).unwrap();
        assert_eq!(p2.get(), 123);
    }

    #[test]
    fn try_from_u64_accepts_nonzero_that_fits() {
        let p = BasePos::try_from(1_u64).unwrap();
        assert_eq!(p.get(), 1);

        let p2 = BasePos::try_from(123_u64).unwrap();
        assert_eq!(p2.get(), 123);
    }

    // --- overflow behavior depends on pointer width ---

    #[cfg(target_pointer_width = "16")]
    #[test]
    fn try_from_u32_rejects_values_that_do_not_fit_on_16bit() {
        // On 16-bit: usize::MAX is 65535, so 65536 should overflow.
        let v: u32 = (u16::MAX as u32) + 1;

        let err = BasePos::try_from(v).unwrap_err();
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

        let err = BasePos::try_from(v).unwrap_err();
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
        let p = BasePos::try_from(v).unwrap();
        assert_eq!(p.get() as u64, v);
    }
}

#[cfg(test)]
mod pos_arith_tests {
    use super::*;

    #[test]
    fn try_add_ok() {
        let p = BasePos::new(10).unwrap();
        let q = p.try_add(5).unwrap();
        assert_eq!(q.get(), 15);
    }

    #[test]
    fn try_sub_ok() {
        let p = BasePos::new(10).unwrap();
        let q = p.try_sub(3).unwrap();
        assert_eq!(q.get(), 7);
    }

    #[test]
    fn try_sub_underflow_to_zero_errors() {
        let p = BasePos::new(1).unwrap();
        let err = p.try_sub(1).unwrap_err();
        match err {
            Error::PositionUnderflow { lhs, rhs, min } => {
                assert_eq!(lhs, p.get());
                assert_eq!(rhs, 1);
                assert_eq!(min, BasePos::MIN.get());
            }
            other => panic!("expected PositionUnderflow, got {other:?}"),
        }
    }

    #[test]
    fn try_sub_underflow_below_zero_errors() {
        let p = BasePos::new(1).unwrap();
        let err = p.try_sub(2).unwrap_err();
        match err {
            Error::PositionUnderflow { lhs, rhs, min } => {
                assert_eq!(lhs, p.get());
                assert_eq!(rhs, 2);
                assert_eq!(min, BasePos::MIN.get());
            }
            other => panic!("expected PositionUnderflow, got {other:?}"),
        }
    }

    #[test]
    fn try_add_overflow_errors() {
        let p = BasePos::MAX;
        let err = p.try_add(1).unwrap_err();
        match err {
            Error::PositionOverflowAdd { lhs, rhs, max } => {
                assert_eq!(lhs, p.get());
                assert_eq!(rhs, 1);
                assert_eq!(max, BasePos::MAX.get());
            }
            other => panic!("expected PositionOverflowAdd, got {other:?}"),
        }
    }
}
