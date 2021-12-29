//! Interval partition

#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::{Arbitrary, Gen};
use super::{Partition, Subdivide, Mid};


/// A half-open interval [a, b) between two points a and b
#[derive(Copy, Clone, Debug)]
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

impl<T: Mid + Copy> Subdivide for Interval<T> {
    fn subdivide(&self) -> Vec<Interval<T>> {
        let mid = self.start.mid(&self.end);
        vec![
            Interval { start: self.start, end: mid },
            Interval { start: mid, end: self.end },
        ]
    }
}

impl<T: Mid + PartialOrd + Copy> Partition<T> for Interval<T> {
    fn contains(&self, elem: &T) -> bool {
        (self.start <= *elem) && (*elem < self.end)
    }
}

#[cfg(any(test, feature = "arbitrary"))]
impl<T: PartialOrd + Arbitrary> Arbitrary for Interval<T> {
    fn arbitrary(g: &mut Gen) -> Interval<T> {
        use std::iter::repeat;
        let a: T = Arbitrary::arbitrary(g);
        let b = repeat(())
            .map(|_| Arbitrary::arbitrary(g))
            .filter(|b: &T| b > &a)
            .next()
            .unwrap();
        Interval::new(a, b)
    }
}


#[cfg(test)]
mod test {
    use super::Interval;

    partition_quickcheck!(interval_f32_partition, Interval<f32>, f32);
    partition_quickcheck!(interval_f64_partition, Interval<f64>, f64);
}
