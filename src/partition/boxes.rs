//! Box partitioning

use nalgebra::{Vector2, Vector3};
#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::{Arbitrary, Gen};
use partition::{Partition, Subdivide, Interval, Mid};
use std::fmt::Debug;


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
        impl<T: Mid + PartialOrd + Copy> Subdivide for $b<T> {
            fn subdivide(&self) -> Vec<$b<T>> {
                let mut subs = vec![];
                $(let $param = self.$param.subdivide();)*
                subs.extend(
                    iproduct!($($param.iter()),*)
                        .map(|($(&$param),*)| $b::new($($param),*))
                );
                subs
            }
        }

        impl<T: Debug + Mid + PartialOrd + Copy + 'static> Partition<$v<T>> for $b<T> {
            fn contains(&self, elem: &$v<T>) -> bool {
                true $(&& self.$param.contains(&elem.$param))*
            }
        }
    )
}

macro_rules! impl_arb_for_box {
    ($b: ident, $($param: ident),*) => (
        #[cfg(any(test, feature = "arbitrary"))]
        impl<T: PartialOrd + Arbitrary> Arbitrary for $b<T> {
            fn arbitrary<G: Gen>(g: &mut G) -> $b<T> {
                $b { $($param: Arbitrary::arbitrary(g)),* }
            }
        }
    )
}


/// A 2d box of intervals
#[derive(Copy, Clone, Debug)]
pub struct Box2<T> {
    x: Interval<T>,
    y: Interval<T>,
}

impl_box!(Box2, x, y);
impl_partition_for_box!(Box2, Vector2, x, y);
impl_arb_for_box!(Box2, x, y);


/// A 3d box of intervals
#[derive(Copy, Clone, Debug)]
pub struct Box3<T> {
    x: Interval<T>,
    y: Interval<T>,
    z: Interval<T>,
}

impl_box!(Box3, x, y, z);
impl_partition_for_box!(Box3, Vector3, x, y, z);
impl_arb_for_box!(Box3, x, y, z);


#[cfg(test)]
mod test {
    pub use nalgebra::{Vector2, Vector3};
    pub use super::*;

    partition_quickcheck!(box2_f32_partition, Box2<f32>, Vector2<f32>);
    partition_quickcheck!(box2_f64_partition, Box2<f64>, Vector2<f64>);
    partition_quickcheck!(box3_f32_partition, Box3<f32>, Vector3<f32>);
    partition_quickcheck!(box3_f64_partition, Box3<f64>, Vector3<f64>);
}
