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


/// The notion of a mid point between two inputs
trait Mid {
    /// Return the mid between this point and another
    fn mid(self, other: Self) -> Self;
}

impl Mid for f64 {
    fn mid(self, other: f64) -> f64 { (self + other) / 2.0 }
}

impl Mid for f32 {
    fn mid(self, other: f32) -> f32 { (self + other) / 2.0 }
}


/// A half-open interval [a, b) between two points a and b
#[derive(Copy)]
pub struct Interval<T> {
    start: T,
    end: T,
}

impl<T: PartialOrd> Interval<T> {
    /// Create a new interval given lower and upper bound
    ///
    /// This constructor dynamically asserts that `start < end`.
    pub fn new(start: T, end: T) -> Interval<T> {
        assert!(start < end);
        Interval { start: start, end: end }
    }
}

impl<T: Mid + PartialOrd + Copy> Partition<T> for Interval<T> {
    type Iter = IntoIter<Interval<T>>;

    fn subdivide(&self) -> IntoIter<Interval<T>> {
        let mid = self.start.mid(self.end);
        vec![
            Interval { start: self.start, end: mid },
            Interval { start: mid, end: self.end },
        ].into_iter()
    }

    fn contains(&self, elem: &T) -> bool {
        (self.start <= *elem) && (*elem < self.end)
    }
}


#[cfg(test)]
mod test {
    use super::{Interval, Partition};
    use quickcheck::{quickcheck, TestResult};

    #[test]
    fn interval_total_f64() {
        fn interval_total_f64((a, b): (f64, f64), c: f64) -> TestResult {
            if b < a { TestResult::discard() }
            else { TestResult::from_bool({
                Interval::new(a, b).prop_is_total(&c)
            })}
        }
        quickcheck(interval_total_f64 as fn((f64, f64), f64) -> TestResult);
    }
}
