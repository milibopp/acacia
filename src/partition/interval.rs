//! Interval partition

use partition::{Partition, Mid};


/// A half-open interval [a, b) between two points a and b
#[derive(Copy, Clone)]
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
mod test {
    use super::Interval;
    use partition::Partition;
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
