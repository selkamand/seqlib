//!  Region type

use crate::coords::InterbaseInterval;

/// A region (name of originating sequence + [`InterbaseInterval()`])
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    name: String,
    interval: InterbaseInterval,
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
    /// Create a new [`Region`] from a sequence **name** and [`InterbaseInterval`]
    pub fn new<N, I>(name: N, interval: I) -> Self
    where
        N: Into<String>,
        I: Into<InterbaseInterval>,
    {
        Self {
            name: name.into(),
            interval: interval.into(),
        }
    }

    /// Get region [`Interval`]
    pub fn interval(&self) -> &InterbaseInterval {
        &self.interval
    }

    /// Name of sequence interval is from (e.g. chromosome / contig /etc)
    pub fn name(&self) -> &str {
        &self.name
    }
}
