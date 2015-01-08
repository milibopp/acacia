//! Abstraction of spatial partitioning schemes

use std::vec::IntoIter;


/// A type describing a partition of some space
pub trait Partition<T>: Sized {
    // TODO: uncomment the following, as soon as this bug is fixed.
    // link: https://github.com/rust-lang/rust/issues/20551

    // The type of subsequent partitions
    // type Subpartition: Partition<T> = Self;

    /// The type of iterator used for subdivision
    type Iter: Iterator<Item=Self>;

    /// Subdivide into smaller partitions
    fn subdivide(&self) -> Self::Iter;

    /// Does the partition contain an element?
    fn contains(&self, elem: &T) -> bool;

    /// Totality proposition
    ///
    /// A partition is required to be total, i.e. any element contained in the
    /// partition will be contained in exactly one of its subdivided partitions.
    /// At the same time, an element not contained in the partition can not be
    /// contained by any subdivided partition.
    fn prop_is_total(&self, elem: &T) -> bool {
        if self.contains(elem) {
            self.subdivide().filter(|sub| sub.contains(elem)).count() == 1
        }
        else {
            self.subdivide().all(|sub| !sub.contains(elem))
        }
    }
}


/// A half-open interval [a, b) between two f64
#[derive(Copy)]
pub struct Interval {
    start: f64,
    end: f64
}

impl Interval {
    /// Create a new interval given lower and upper bound
    ///
    /// This constructor dynamically asserts that `start < end`.
    pub fn new(start: f64, end: f64) -> Interval {
        assert!(start < end);
        Interval { start: start, end: end }
    }
}

impl Partition<f64> for Interval {
    type Iter = IntoIter<Interval>;

    fn subdivide(&self) -> IntoIter<Interval> {
        let mid = (self.start + self.end) / 2.0;
        vec![
            Interval::new(self.start, mid),
            Interval::new(mid, self.end)
        ].into_iter()
    }

    fn contains(&self, elem: &f64) -> bool {
        (self.start <= *elem) && (*elem < self.end)
    }
}


#[cfg(test)]
mod test {
    use super::{Interval, Partition};
    use quickcheck::{quickcheck, TestResult};

    #[test]
    fn interval_total() {
        fn interval_total((a, b): (f64, f64), c: f64) -> TestResult {
            if b < a { TestResult::discard() }
            else { TestResult::from_bool({
                Interval::new(a, b).prop_is_total(&c)
            })}
        }
        quickcheck(interval_total as fn((f64, f64), f64) -> TestResult);
    }
}
