//! N-cube or hypercube partitioning scheme

use std::ops::{Index, IndexMut};
use std::num::{Int, cast};
use std::cmp::PartialOrd;
use std::iter::AdditiveIterator;
use nalgebra::{Dim, BaseFloat, Zero, zero};
use partition::Partition;


/// An N-cube based partitioning scheme
#[derive(Copy, Clone)]
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

impl<P: Copy, S: Copy> Ncube<P, S> {
    /// The width of the N-cube
    pub fn width(&self) -> S { self.width }

    /// The center of the N-cube
    pub fn center(&self) -> P { self.center }
}

impl<P, S> Partition<P> for Ncube<P, S>
    where P: Dim + Index<uint, Output=S> + IndexMut<uint, Output=S> + Copy,
          S: BaseFloat + PartialOrd,
{
    fn subdivide(&self) -> Vec<Ncube<P, S>> {
        let _2 = cast(2.0f64).unwrap();
        let dim = Dim::dim(None::<P>);
        let new_width = self.width / _2;
        range(0u, 2.pow(dim))
            .map(|n| {
                let mut new_center = self.center;
                let dx = new_width / _2;
                for i in range(0, dim) {
                    new_center[i] = new_center[i] + match n / 2.pow(i) % 2 {
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

    fn contains(&self, elem: &P) -> bool {
        let _2 = cast(2.0f64).unwrap();
        range(0u, Dim::dim(None::<P>))
            .all(|i| {
                let off = (self.center[i] - elem[i]) * _2;
                (-self.width <= off) && (off < self.width)
            })
    }

    fn dispatch(&self, elem: &P) -> uint {
        range(0, Dim::dim(None::<P>))
            .map(|k| if elem[k] < self.center[k] {0} else {1 << k})
            .sum()
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use partition::Partition;
    use nalgebra::Pnt2;
    use quickcheck::{quickcheck, TestResult};

    #[test]
    fn ncube_total() {
        fn ncube_total(center: (f64, f64), width: f64, elem: (f64, f64)) -> TestResult {
            if width > 0.0 {
                let center = Pnt2::new(center.0, center.1);
                let elem = Pnt2::new(elem.0, elem.1);
                TestResult::from_bool(Ncube::new(center, width).prop_is_total(&elem))
            }
            else {
                TestResult::discard()
            }
        }
        quickcheck(ncube_total as fn((f64, f64), f64, (f64, f64)) -> TestResult);
    }
}
