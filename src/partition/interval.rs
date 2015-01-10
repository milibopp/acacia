//! Interval partition

use std::iter::repeat;
#[cfg(test)]
use quickcheck::{Arbitrary, Gen};
use partition::{Partition, Mid};


/// A half-open interval [a, b) between two points a and b
#[derive(Copy, Clone, Show)]
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
    fn subdivide(&self) -> Vec<Interval<T>> {
        let mid = self.start.mid(&self.end);
        vec![
            Interval { start: self.start, end: mid },
            Interval { start: mid, end: self.end },
        ]
    }

    fn contains(&self, elem: &T) -> bool {
        (self.start <= *elem) && (*elem < self.end)
    }
}

#[cfg(test)]
impl<T: PartialOrd + Arbitrary> Arbitrary for Interval<T> {
    fn arbitrary<G: Gen>(g: &mut G) -> Interval<T> {
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
    use partition::{Partition, prop_is_total};
    use quickcheck::{quickcheck, TestResult};

    #[test]
    fn interval_total_f32() {
        quickcheck(prop_is_total as fn(Interval<f32>, f32) -> bool);
    }

    #[test]
    fn interval_total_f64() {
        quickcheck(prop_is_total as fn(Interval<f64>, f64) -> bool);
    }
}
