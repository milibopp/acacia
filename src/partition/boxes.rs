//! Box partitioning

use nalgebra::{Vec2, Vec3};
use partition::{Partition, Interval, Mid};


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
