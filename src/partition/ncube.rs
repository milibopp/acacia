//! N-cube or hypercube partitioning scheme

use std::ops::{Index, IndexMut};
use std::num::{Int, cast};
use std::cmp::PartialOrd;
use std::iter::AdditiveIterator;
use nalgebra::{Dim, BaseFloat, Zero, zero};
#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::{Arbitrary, Gen};
use partition::{Partition, Subdivide};


/// An N-cube based partitioning scheme
#[derive(Copy, Clone, Debug)]
pub struct Ncube<P, S> {
    center: P,
    width: S,
}

impl<P, S: PartialOrd + Zero> Ncube<P, S> {
    /// Create a new N-cube given its center and width
    pub fn new(center: P, width: S) -> Ncube<P, S> {
        assert!(width > zero());
        Ncube { center: center, width: width }
    }
}

impl<P, S: Clone> Ncube<P, S> {
    /// The width of the N-cube
    pub fn width(&self) -> S { self.width.clone() }
}

impl<P: Clone, S> Ncube<P, S> {
    /// The center of the N-cube
    pub fn center(&self) -> P { self.center.clone() }
}

impl<P, S> Subdivide for Ncube<P, S>
    where P: Dim + Index<usize, Output=S> + IndexMut<usize, Output=S> + Copy,
          S: BaseFloat + PartialOrd,
{
    fn subdivide(&self) -> Vec<Ncube<P, S>> {
        let _2 = cast(2.0f64).unwrap();
        let dim = Dim::dim(None::<P>);
        let new_width = self.width / _2;
        (0..2.pow(dim as u32))
            .map(|n: i32| {
                let mut new_center = self.center;
                let dx = new_width / _2;
                for i in 0..dim {
                    new_center[i] = new_center[i] + match n / 2.pow(i as u32) % 2 {
                        0 => -dx,
                        1 => dx,
                        _ => unreachable!(),
                    };
                }
                Ncube {
                    center: new_center,
                    width: new_width,
                }
            })
        .collect()
    }
}

impl<P, S> Partition<P> for Ncube<P, S>
    where P: Dim + Index<usize, Output=S> + IndexMut<usize, Output=S> + Copy,
          S: BaseFloat + PartialOrd,
{
    fn contains(&self, elem: &P) -> bool {
        let _2 = cast(2.0f64).unwrap();
        (0..Dim::dim(None::<P>))
            .all(|i| {
                let off = (self.center[i] - elem[i]) * _2;
                (-self.width <= off) && (off < self.width)
            })
    }

    fn dispatch(&self, elem: &P) -> usize {
        (0..Dim::dim(None::<P>))
            .map(|k| if elem[k] < self.center[k] {0} else {1 << k})
            .sum()
    }
}

#[cfg(any(test, feature = "arbitrary"))]
impl<P: Arbitrary, S: PartialOrd + Zero + Arbitrary> Arbitrary for Ncube<P, S> {
    fn arbitrary<G: Gen>(g: &mut G) -> Ncube<P, S> {
        use std::iter::repeat;
        Ncube::new(
            Arbitrary::arbitrary(g),
            repeat(())
                .map(|_| Arbitrary::arbitrary(g))
                .filter(|w: &S| w > &zero())
                .next()
                .unwrap()
        )
    }
}


#[cfg(test)]
mod test {
    pub use nalgebra::Pnt2;
    pub use super::*;

    partition_quickcheck!(ncube_pnt2_f32_partition, Ncube<Pnt2<f32>, f32>, Pnt2<f32>);
}
