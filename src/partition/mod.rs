//! Abstraction of spatial partitioning schemes

pub use partition::interval::Interval;
pub use partition::boxes::{Box2, Box3};
pub use partition::ncube::Ncube;

mod interval;
mod boxes;
mod ncube;


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
