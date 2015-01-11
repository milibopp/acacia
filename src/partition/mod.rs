//! Abstraction of spatial partitioning schemes

pub use partition::interval::Interval;
pub use partition::boxes::{Box2, Box3};
pub use partition::ncube::Ncube;
pub use partition::unitquad::UnitQuad;


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


/// Totality proposition
///
/// A partition is required to be total, i.e. any element contained in the
/// partition will be contained in exactly one of its subdivided partitions.
/// At the same time, an element not contained in the partition can not be
/// contained by any subdivided partition.
pub fn prop_is_total<P: Partition<T>, T>(partition: P, elem: T) -> bool {
    let subs = partition.subdivide();
    if partition.contains(&elem) {
        subs.iter().filter(|sub| sub.contains(&elem)).count() == 1
    }
    else {
        subs.iter().all(|sub| !sub.contains(&elem))
    }
}


/// Consistency of dispatch mechanism
///
/// An element contained in the partition must be dispatched consistently, i.e.
/// the subpartition to which it is dispatched must contain it.
pub fn prop_consistent_dispatch<P: Partition<T>, T>(partition: P, elem: T) -> bool {
    if partition.contains(&elem) {
        partition.subdivide()[partition.dispatch(&elem)].contains(&elem)
    }
    else {
        true
    }
}


#[macro_escape]
macro_rules! partition_quickcheck (
    ($testfn: ident, $p: ty, $t: ty) => (
        mod $testfn {
            use $crate::partition::{prop_is_total, prop_consistent_dispatch};
            use quickcheck::quickcheck;
            use super::*;

            #[test]
            fn is_total() {
                quickcheck(prop_is_total as fn($p, $t) -> bool);
            }

            #[test]
            fn consistent_dispatch() {
                quickcheck(prop_consistent_dispatch as fn($p, $t) -> bool);
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
