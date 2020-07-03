//! N-cube or hypercube partitioning scheme

use nalgebra::base::{allocator::Allocator, default_allocator::DefaultAllocator};
use nalgebra::{zero, DimName, Point, Scalar};
use num_traits::{Float, NumCast, PrimInt};
use partition::{Partition, Subdivide};
#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::{Arbitrary, Gen};
use std::cmp::PartialOrd;
use std::fmt::Debug;

/// An N-cube based partitioning scheme
#[derive(Clone, Debug)]
pub struct Ncube<D: DimName, S: Scalar + Copy + Debug>
where
    DefaultAllocator: Allocator<S, D>,
{
    center: Point<S, D>,
    width: S,
}

impl<D: DimName, S: Scalar + Copy> Copy for Ncube<D, S>
where
    Point<S, D>: Copy,
    DefaultAllocator: Allocator<S, D>,
{
}

impl<D: DimName, S: Debug + PartialOrd + Float + Copy + 'static> Ncube<D, S>
where
    Point<S, D>: Copy,
    DefaultAllocator: Allocator<S, D>,
{
    /// Create a new N-cube given its center and width
    pub fn new(center: Point<S, D>, width: S) -> Ncube<D, S> {
        assert!(width > zero());
        Ncube {
            center: center,
            width: width,
        }
    }
}

impl<D: DimName, S: Debug + PartialEq + Copy + 'static> Ncube<D, S>
where
    Point<S, D>: Copy,
    DefaultAllocator: Allocator<S, D>,
{
    /// The width of the N-cube
    pub fn width(&self) -> S {
        self.width
    }
}

impl<D: DimName, S: Debug + PartialEq + Copy + 'static> Ncube<D, S>
where
    DefaultAllocator: Allocator<S, D>,
{
    /// The center of the N-cube
    pub fn center(&self) -> Point<S, D> {
        self.center.clone()
    }
}

impl<D, S> Subdivide for Ncube<D, S>
where
    D: DimName,
    S: Scalar + PartialOrd + NumCast + Copy + Float,
    Point<S, D>: Copy,
    DefaultAllocator: Allocator<S, D>,
{
    fn subdivide(&self) -> Vec<Ncube<D, S>> {
        let _2 = NumCast::from(2.0f64).unwrap();
        let dimension = D::dim();
        let new_width = self.width / _2;
        (0..2.pow(dimension as u32))
            .map(|n: i32| {
                let mut new_center = self.center;
                let dx = new_width / _2;
                for i in 0..dimension {
                    new_center[i] = new_center[i]
                        + match n / 2.pow(i as u32) % 2 {
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

impl<D, S> Partition<Point<S, D>> for Ncube<D, S>
where
    D: DimName,
    S: Scalar + PartialOrd + NumCast + Copy + Float,
    Point<S, D>: Copy,
    DefaultAllocator: Allocator<S, D>,
{
    fn contains(&self, elem: &Point<S, D>) -> bool {
        let _2 = NumCast::from(2.0f64).unwrap();
        (0..D::dim()).all(|i| {
            let off = (self.center[i] - elem[i]) * _2;
            (-self.width <= off) && (off < self.width)
        })
    }

    fn dispatch(&self, elem: &Point<S, D>) -> usize {
        (0..D::dim())
            .map(|k| if elem[k] < self.center[k] { 0 } else { 1 << k })
            .fold(0, |a, b| a + b)
    }
}

#[cfg(any(test, feature = "arbitrary"))]
impl<D: DimName, S: PartialOrd + Float + Arbitrary + Debug + Copy> Arbitrary for Ncube<D, S>
where
    Point<S, D>: Copy,
    DefaultAllocator: Allocator<S, D>,
    <DefaultAllocator as Allocator<S, D>>::Buffer: Send,
{
    fn arbitrary<G: Gen>(g: &mut G) -> Ncube<D, S> {
        use std::iter::repeat;
        Ncube::new(
            Arbitrary::arbitrary(g),
            repeat(())
                .map(|_| Arbitrary::arbitrary(g))
                .filter(|w: &S| w > &zero())
                .next()
                .unwrap(),
        )
    }
}

#[cfg(test)]
mod test {
    pub use super::*;
    pub use nalgebra::{Point, base::dimension::U2};

    partition_quickcheck!(
        ncube_pnt2_f32_partition,
        Ncube<U2, f32>,
        Point<f32, U2>
    );
}
