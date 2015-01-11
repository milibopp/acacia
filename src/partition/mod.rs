//! Abstraction of spatial partitioning schemes

pub use partition::interval::Interval;
pub use partition::boxes::{Box2, Box3};
pub use partition::ncube::Ncube;
pub use partition::unitquad::UnitQuad;

#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::TestResult;

/// A type describing a partition of some space
///
/// In addition to this trait signature an implementation is required to satisfy
/// `prop_is_total`.
pub trait Partition<T>: Sized {
    // TODO: uncomment the following, as soon as this bug is fixed.
    // link: https://github.com/rust-lang/rust/issues/20551

    // The type of subsequent partitions
    // type Subpartition: Partition<T> = Self;

    /// Subdivide into smaller partitions
    fn subdivide(&self) -> Vec<Self>;

    /// Does the partition contain an element?
    fn contains(&self, elem: &T) -> bool;

    /// Dispatch an element to the correct subpartition
    ///
    /// The default implementation works, if the totality proposition is
    /// fulfilled. However, note that its performance is not optimal, as it
    /// checks for all subpartitions whether they contain the element, until one
    /// is found that does.
    fn dispatch(&self, elem: &T) -> usize {
        for (i, part) in self.subdivide().iter().enumerate() {
            if part.contains(elem) {
                return i;
            }
        }
        panic!("partition dispatch impossible");
    }
}


/// Totality proposition (contained)
///
/// Any element contained in the partition will be contained in exactly one of
/// its subpartitions.
#[cfg(any(test, feature = "arbitrary"))]
pub fn prop_is_total_elem<P: Partition<T>, T>(partition: P, elem: T) -> TestResult {
    if partition.contains(&elem) {
        TestResult::from_bool(
            partition.subdivide().iter()
                .filter(|sub| sub.contains(&elem)).count()
            == 1
        )
    }
    else {
        TestResult::discard()
    }
}


/// Totality proposition (not contained)
///
/// If a partition does not contain an element, none of its subpartitions
/// contain it either.
#[cfg(any(test, feature = "arbitrary"))]
pub fn prop_is_total_non_elem<P: Partition<T>, T>(partition: P, elem: T) -> TestResult {
    if partition.contains(&elem) {
        TestResult::discard()
    }
    else {
        TestResult::from_bool(
            partition.subdivide().iter()
                .all(|sub| !sub.contains(&elem))
        )
    }
}


/// Consistency of dispatch mechanism
///
/// An element contained in the partition must be dispatched consistently, i.e.
/// the subpartition to which it is dispatched must contain it.
#[cfg(any(test, feature = "arbitrary"))]
pub fn prop_consistent_dispatch<P: Partition<T>, T>(partition: P, elem: T) -> TestResult {
    if partition.contains(&elem) {
        TestResult::from_bool(partition.subdivide()[partition.dispatch(&elem)].contains(&elem))
    }
    else {
        TestResult::discard()
    }
}


#[macro_escape]
macro_rules! partition_quickcheck (
    ($testfn: ident, $p: ty, $t: ty) => (
        mod $testfn {
            use $crate::partition::{
                prop_is_total_elem, prop_is_total_non_elem,
                prop_consistent_dispatch
            };
            use quickcheck::{quickcheck, TestResult};
            use super::*;

            #[test]
            fn is_total_elem() {
                quickcheck(prop_is_total_elem as fn($p, $t) -> TestResult);
            }

            #[test]
            fn is_total_non_elem() {
                quickcheck(prop_is_total_non_elem as fn($p, $t) -> TestResult);
            }

            #[test]
            fn consistent_dispatch() {
                quickcheck(prop_consistent_dispatch as fn($p, $t) -> TestResult);
            }
        }
    )
);


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


mod interval;
mod boxes;
mod ncube;
mod unitquad;
