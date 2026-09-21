//! Interval Type (simple range defined by start and end position)

use std::num::NonZeroUsize;

use crate::coords::{BasePos, InterbasePos};
use crate::error::CoordError as Error;
pub(crate) type Result<T> = std::result::Result<T, Error>;

/// A 0-based inter-residue interval  
///
/// 0 is the position before the first residue in a sequence
///
/// For the numbering of a the 3 base sequence:
///
/// ```text
///  A C T
/// 0 1 2 3
/// ```
///
/// The interval describing the full sequence is 0-3
///
/// Empty intervals are valid and represent a boundary between residues. This is useful
/// for describing insertions, where the reference span has length 0.
///
/// # Examples
/// ```
/// use seqlib::coords::{InterbaseInterval, InterbasePos};
/// use std::ops::Range;
///
/// // Create interval from two InterbasePositions
/// let i = InterbaseInterval::new(InterbasePos::from(0), InterbasePos::from(3)).expect("valid interval");
///
/// // Create equivalent interval from a rust Range
/// let i2 = InterbaseInterval::from(0..3);
///
/// // Note if creating from a rust Range, if Start > End the conversion automatically flips the
/// // coords so Start <= End. However when using [`Interbase::new()`] this case would return an explicit error
/// let i3 = InterbaseInterval::from(3..0);
///
/// // Ensure all InterbaseIntervals are the same
/// assert_eq!(i, i2);
/// assert_eq!(i, i3);
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct InterbaseInterval {
    start: InterbasePos,
    end: InterbasePos,
}

impl InterbaseInterval {
    /// Create a new interbase, 0-based [`InterbaseInterval`].
    ///
    /// ```text
    ///       A   T   A   C
    ///     0   1   2   3   4
    ///    └─────┘
    /// ```
    ///
    ///
    /// ```text
    ///       A   T   A   C
    ///     0   1   2   3   4
    ///    └─────┘
    /// ```
    ///
    /// # Invariants
    /// Start must be less than or equal to end. Empty intervals are allowed.
    ///
    pub fn new(start: InterbasePos, end: InterbasePos) -> Result<Self> {
        if start > end {
            return Err(Error::InvalidIntervalCoords {
                start: start.get(),
                end: end.get(),
            });
        }

        // Return a valid interval
        Ok(Self { start, end })
    }

    /// Create a new inter-residue, zero-based [`Interval`]
    ///
    /// Skips check  `start <= end`. Used in macros.
    ///
    ///
    /// ```text
    ///       A   T   A   C
    ///     0   1   2   3   4
    ///    └─────┘
    /// ```
    pub(crate) fn new_unchecked(start: InterbasePos, end: InterbasePos) -> Self {
        Self { start, end }
    }

    /// Creates an interval around `pos` with `left` bases before and `right` bases after it.
    ///
    /// The bounds saturate at [`InterbasePos::MIN`] and [`InterbasePos::MAX`] rather than failing on
    /// underflow or overflow.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{InterbaseInterval, InterbasePos};
    ///
    /// let interval = InterbaseInterval::around_position(InterbasePos::from(10usize), 2, 3);
    ///
    /// assert_eq!(*interval.start(), InterbasePos::from(8));
    /// assert_eq!(*interval.end(), InterbasePos::from(13));
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn around_position(pos: InterbasePos, left: usize, right: usize) -> Self {
        Self {
            start: pos.saturating_sub(left),
            end: pos.saturating_add(right),
        }
    }

    /// Returns the start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{InterbaseInterval, InterbasePos};
    ///
    /// let interval = InterbaseInterval::new(InterbasePos::from(3), InterbasePos::from(7))?;
    ///
    /// assert_eq!(*interval.start(), InterbasePos::from(3));
    ///
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn start(&self) -> &InterbasePos {
        &self.start
    }

    /// Returns the end position.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{InterbaseInterval, InterbasePos};
    ///
    /// let interval = InterbaseInterval::new(InterbasePos::from(2), InterbasePos::from(7)).unwrap();
    ///
    /// assert_eq!(*interval.end(), InterbasePos::from(7));
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn end(&self) -> &InterbasePos {
        &self.end
    }

    /// Check if interval is empty.
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }

    /// Returns the number of bases spanned by the interval.
    ///
    /// Because intervals are interbase, `0-1` has length 1.
    ///
    /// ```text
    ///       A   T   A   C
    ///     0   1   2   3   4
    ///    └─────┘
    /// ```
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{InterbaseInterval, InterbasePos};
    ///
    /// let interval = InterbaseInterval::new(InterbasePos::from(2), InterbasePos::from(4))?;
    ///
    /// assert_eq!(interval.len(), 2);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn len(&self) -> usize {
        // this cannot overflow because during construction we ensure end >= start
        self.end().get() - self.start().get()
    }

    /// Converts a feature interval (`original`) to offsets within the region `self`.
    ///
    /// Use this to locate a variant or annotation in a sequence extracted from `self`.
    /// Both inputs use parent-sequence coordinates; the returned interval treats
    /// the start of `self` as position `0`. (See example for clearer examples)
    ///
    /// Returns `None` unless `original` is fully contained within `self`.
    ///
    /// # Examples
    ///  
    /// ```text
    ///                        A   T   A   C   G
    /// Parent coordinates:   0   1   2   3   4   5
    /// Region (self):                └───────────┘    2-5
    /// Feature (original):               └───┘        3-4
    /// Local coordinates:            0   1   2   3
    /// Local interval:                   └───┘        1-2 <- Returned
    /// ```
    ///
    /// ```
    /// use seqlib::coords::{InterbaseInterval, InterbasePos};
    ///
    /// // A window described in chromosome coordinates, and its extracted bases.
    /// let window = InterbaseInterval::new(InterbasePos::from(2), InterbasePos::from(5))?;
    ///
    /// // The variant's reference interval, also in chromosome coordinates.
    /// let variant = InterbaseInterval::new(InterbasePos::from(3), InterbasePos::from(4))?;
    ///
    /// // To describe the variant location within `window` interval (e.g. turn 3–4 into 2–3)
    /// // we use the local_interval method
    /// let local = window.local_interval(variant).unwrap();
    ///
    /// assert_eq!(local, InterbaseInterval::new(InterbasePos::from(1), InterbasePos::from(2))?);
    ///
    /// # Ok::<(), seqlib::error::CoordError>(())    
    /// ```
    ///
    pub fn local_interval(&self, original: InterbaseInterval) -> Option<InterbaseInterval> {
        // Check supplied interval is contained within self (recall that our constructor guarantees
        // start < end so we only need to do these two checks
        if original.start < self.start || original.end > self.end {
            return None;
        }

        let local_start = original.start.get() - self.start.get();
        let local_end = local_start + original.len();

        let local_interval = InterbaseInterval::new_unchecked(
            InterbasePos::new(local_start),
            InterbasePos::new(local_end),
        );

        Some(local_interval)
    }
}

/// Converts an [`InterbaseInterval`] into the equivalent zero-based,
/// half open suitable for indexing into a rust array. Note rust Range is effectively equivalent to
/// InterbaseIntervals, so conversion is very straightforward
///
/// ```text
/// Seq:                    A C T
/// InterbasePos:          0 1 2 3
/// Rust:                  0 1 2 3
/// InterbaseInterval 1-3:    ---
/// Rust Range 1..3:          ---
/// ```
///
impl From<InterbaseInterval> for std::ops::Range<usize> {
    fn from(val: InterbaseInterval) -> Self {
        val.start.get()..val.end.get()
    }
}

impl From<&InterbaseInterval> for std::ops::Range<usize> {
    fn from(val: &InterbaseInterval) -> Self {
        val.start.get()..val.end.get()
    }
}

/// Converts a rust [`std::ops::Range`] into an equivalent [`InterbaseInterval`].
///
/// It convenient to think of rust ranges as interbase, 0-start inter-element coordinates equivalent
/// to InterbaseIntervals. The key difference is that `Range` types can have start > end but
/// InterbaseIntervals can not. To deal with this, when Range end > start we simply switch them when
/// loading into an InterbaseInterval.
///
// ```text
/// Seq:                    A C T
/// InterbasePos:          0 1 2 3
/// Rust:                  0 1 2 3
/// InterbaseInterval 1-3:    ---
/// Rust Range 1..3:          ---
/// ```
///
/// ```
/// use seqlib::coords::{InterbaseInterval, InterbasePos};
/// use std::ops::Range;
///
/// let interval = InterbaseInterval::from(1..3);
/// ```
impl From<std::ops::Range<usize>> for InterbaseInterval {
    fn from(val: std::ops::Range<usize>) -> Self {
        if val.start <= val.end {
            Self::new(InterbasePos::new(val.start), InterbasePos::new(val.end))
                .expect("Bug in conversion: Creating an InterbaseInterval from a Range<usize> should never error. Please report this bug")
        } else {
            Self::new(InterbasePos::new(val.end), InterbasePos::new(val.start))
                .expect("Bug in conversion: Creating an InterbaseInterval from a Range<usize> should never error. Please report this bug")
        }
    }
}

/// A genomic interval (Start & End)
/// Both are 1-based and both-end inclusive
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseInterval {
    start: BasePos,
    end: BasePos,
}

impl std::fmt::Display for BaseInterval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", self.start, self.end)
    }
}

impl std::fmt::Display for InterbaseInterval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", self.start, self.end)
    }
}

impl BaseInterval {
    /// Creates a 1-based both-end inclusive residue interval from `start` to `end`
    ///
    /// For the numbering of a the 3 base sequence:
    ///  A C T
    ///  1 2 3
    ///
    /// The BaseInterval describing the full sequence is 1-3
    ///
    /// # Errors
    ///
    /// Returns [`Error::RangeEndTooSmall`] if `end` is less than `start`.
    //
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    ///
    /// let interval = BaseInterval::new(BasePos::new(2usize)?, BasePos::new(5)?)?;
    ///
    /// assert_eq!(*interval.start(), BasePos::new(2)?);
    /// assert_eq!(*interval.end(), BasePos::new(5)?);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    ///
    pub fn new(start: BasePos, end: BasePos) -> Result<Self> {
        if end < start {
            return Err(Error::RangeEndTooSmall {
                start: start.into(),
                end: end.into(),
            });
        }
        Ok(Self { start, end })
    }

    /// Creates an interval around `pos` with `left` bases before and `right` bases after it.
    ///
    /// The bounds saturate at [`BasePos::MIN`] and [`BasePos::MAX`] rather than failing on
    /// underflow or overflow.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    ///
    /// let interval = BaseInterval::around_position(BasePos::new(10)?, 2, 3);
    ///
    /// assert_eq!(*interval.start(), BasePos::new(8)?);
    /// assert_eq!(*interval.end(), BasePos::new(13)?);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn around_position(pos: BasePos, left: usize, right: usize) -> Self {
        Self {
            start: pos.saturating_sub(left),
            end: pos.saturating_add(right),
        }
    }

    /// Returns the inclusive start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    ///
    /// let interval = BaseInterval::new(BasePos::new(3)?, BasePos::new(7)?)?;
    ///
    /// assert_eq!(*interval.start(), BasePos::new(3)?);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn start(&self) -> &BasePos {
        &self.start
    }

    /// Returns the inclusive end position.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    /// use seqlib::basepos;
    ///
    /// let interval = BaseInterval::new(basepos!(2), basepos!(7)).unwrap();
    ///
    /// assert_eq!(*interval.end(), basepos!(7));
    /// ```
    pub fn end(&self) -> &BasePos {
        &self.end
    }

    /// Check if the interval is empty.
    ///
    /// Always returns false because base intervals are non-empty by construction.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Returns the number of positions spanned by the interval.
    ///
    /// Because intervals are inclusive, `1-1` has length 1.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    ///
    /// let interval = BaseInterval::new(BasePos::new(2)?, BasePos::new(5)?)?;
    ///
    /// assert_eq!(interval.len(), 4);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn len(&self) -> usize {
        // this cannot overflow because during construction we ensure end >= start
        self.end().get() - self.start().get() + 1
    }

    /// Returns the number of positions spanned by the interval.
    ///
    /// Because intervals are inclusive, `1-1` has length 1.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    /// use std::num::NonZeroUsize;
    /// let interval = BaseInterval::new(BasePos::new(2usize)?, BasePos::new(2usize)?).unwrap();
    ///
    /// assert_eq!(interval.len_nonzero(), NonZeroUsize::new(1).unwrap());
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ````
    ///
    pub fn len_nonzero(&self) -> NonZeroUsize {
        match NonZeroUsize::try_from(self.len()) {
            Ok(val) => val,
            Err(_) => unreachable!(
                "Implementation mistake: len_nonzero method of interval should never error because len() of a 1-based inclusive interval is always >=1 so long as constructor properly asserts end >= start and len() method calculates length correctly. Please report this error message on this repos github"
            ),
        }
    }

    /// Returns the interval as 0-based half-open indices.
    ///
    /// The returned `(start, end)` pair is suitable for Rust slicing, where
    /// `start` is inclusive and `end` is exclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    ///
    /// let interval = BaseInterval::new(BasePos::new(2)?, BasePos::new(5)?)?;
    ///
    /// assert_eq!(interval.as_0based_indices(), (1, 5));
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn as_0based_indices(&self) -> (usize, usize) {
        (self.start.as_0based_index(), self.end.as_0based_index() + 1)
    }

    /// Returns the 1-based local position of `pos` within the interval.
    ///
    /// The returned [`BasePos`] is suitable for APIs that expect a sequence-local
    /// coordinate, such as the `anchor` argument of
    /// [`MutationWithContext::new`](crate::mutations::MutationWithContext::new).
    /// Returns `None` if `pos` is outside the interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    ///
    /// let interval = BaseInterval::new(BasePos::new(8)?, BasePos::new(13)?)?;
    ///
    /// assert_eq!(interval.local_position(BasePos::new(10)?), Some(BasePos::new(3)?));
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    ///
    /// Saturated intervals still return the observed local position.
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    ///
    /// let interval = BaseInterval::around_position(BasePos::new(2)?, 5, 4);
    ///
    /// assert_eq!(*interval.start(), BasePos::new(1)?);
    /// assert_eq!(*interval.end(), BasePos::new(6)?);
    /// assert_eq!(interval.local_position(BasePos::new(2)?), Some(BasePos::new(2)?));
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    ///
    /// Positions outside the interval return `None`.
    ///
    /// ```
    /// use seqlib::coords::{BaseInterval, BasePos};
    ///
    /// let interval = BaseInterval::new(BasePos::new(8)?, BasePos::new(13)?)?;
    ///
    /// assert_eq!(interval.local_position(BasePos::new(7)?), None);
    /// assert_eq!(interval.local_position(BasePos::new(14)?), None);
    ///
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn local_position(&self, pos: BasePos) -> Option<BasePos> {
        if pos < self.start || pos > self.end {
            return None;
        }

        let local_position = pos.get() - self.start.get() + 1;
        BasePos::new(local_position).ok()
    }
}

impl Default for BaseInterval {
    fn default() -> Self {
        Self {
            start: BasePos::MIN,
            end: BasePos::MIN,
        }
    }
}

/// Converts an [`BaseInterval`] into the equivalent zero-based,
/// half-open interval suitable for indexing into a rust array.
///
/// ```text
/// Seq:                 A C T
/// BasePos:             1 2 3
/// Rust:                0 1 2
/// BaseInterval 2-3:      ---
/// Rust Range 1..3:       ---
/// ```
///
impl From<BaseInterval> for std::ops::Range<usize> {
    fn from(value: BaseInterval) -> Self {
        value.start.get() - 1..value.end.get()
    }
}

impl From<&BaseInterval> for std::ops::Range<usize> {
    fn from(value: &BaseInterval) -> Self {
        value.start.get() - 1..value.end.get()
    }
}
