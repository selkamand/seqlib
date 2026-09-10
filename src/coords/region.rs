//!  Region type

use crate::coords::BaseInterval;

/// A region (name of originating sequence + 1-based inclusive interval)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    name: String,
    interval: BaseInterval,
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{}",
            self.name(),
            self.interval().start(),
            self.interval().end()
        )
    }
}

impl Region {
    /// Create a new [`Region`] from a sequence **name** and [`Interval`]
    pub fn new<N, I>(name: N, interval: I) -> Self
    where
        N: Into<String>,
        I: Into<BaseInterval>,
    {
        Self {
            name: name.into(),
            interval: interval.into(),
        }
    }

    /// Get region [`Interval`]
    pub fn interval(&self) -> &BaseInterval {
        &self.interval
    }

    /// Name of sequence interval is from (e.g. chromosome / contig /etc)
    pub fn name(&self) -> &str {
        &self.name
    }
}
