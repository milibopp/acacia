//! Abstraction of spatial partitioning schemes

use nalgebra::{Vec2, Vec3};


/// A type describing a partition of some space
pub trait Partition<T>: Sized {
    // TODO: uncomment the following, as soon as this bug is fixed.
    // link: https://github.com/rust-lang/rust/issues/20551

    // The type of subsequent partitions
    // type Subpartition: Partition<T> = Self;

    /// Subdivide into smaller partitions
    fn subdivide(&self) -> Vec<Self>;

    /// Does the partition contain an element?
    fn contains(&self, elem: &T) -> bool;

    /// Totality proposition
    ///
    /// A partition is required to be total, i.e. any element contained in the
    /// partition will be contained in exactly one of its subdivided partitions.
    /// At the same time, an element not contained in the partition can not be
    /// contained by any subdivided partition.
    fn prop_is_total(&self, elem: &T) -> bool {
        let subs = self.subdivide();
        if self.contains(elem) {
            subs.iter().filter(|sub| sub.contains(elem)).count() == 1
        }
        else {
            subs.iter().all(|sub| !sub.contains(elem))
        }
    }
}


/// The notion of a mid point between two inputs
trait Mid {
    /// Return the mid between this point and another
    fn mid(&self, other: &Self) -> Self;
}

impl Mid for f64 {
    fn mid(&self, other: &f64) -> f64 { (*self + *other) / 2.0 }
}

impl Mid for f32 {
    fn mid(&self, other: &f32) -> f32 { (*self + *other) / 2.0 }
}


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


macro_rules! impl_box {
    ($b: ident, $($param: ident),*) => (
        impl<T> $b<T> {
            /// Create a new box from intervals
            pub fn new($($param: Interval<T>),*) -> $b<T> {
                $b { $($param: $param),* }
            }
        }
    )
}

macro_rules! impl_partition_for_box {
    ($b: ident, $v: ident, $($param: ident),*) => (
        impl<T: Mid + PartialOrd + Copy> Partition<$v<T>> for $b<T> {
            fn subdivide(&self) -> Vec<$b<T>> {
                let mut subs = vec![];
                $(let $param = self.$param.subdivide();)*
                subs.extend(
                    iproduct!($($param.iter()),*)
                        .map(|($(&$param),*)| $b::new($($param),*))
                );
                subs
            }

            fn contains(&self, elem: &$v<T>) -> bool {
                true $(&& self.$param.contains(&elem.$param))*
            }
        }
    )
}


/// A 2d box of intervals
#[derive(Copy, Clone)]
pub struct Box2<T> {
    x: Interval<T>,
    y: Interval<T>,
}

impl_box!(Box2, x, y);
impl_partition_for_box!(Box2, Vec2, x, y);


/// A 3d box of intervals
#[derive(Copy, Clone)]
pub struct Box3<T> {
    x: Interval<T>,
    y: Interval<T>,
    z: Interval<T>,
}

impl_box!(Box3, x, y, z);
impl_partition_for_box!(Box3, Vec3, x, y, z);


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
