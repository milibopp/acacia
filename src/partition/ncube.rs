//! N-cube or hypercube partitioning scheme

use std::ops::{Index, IndexMut};
use std::num::{Int, cast};
use nalgebra::{Dim, BaseFloat};
use partition::Partition;


/// An N-cube based partitioning scheme
#[derive(Copy, Clone)]
pub struct Ncube<P, S> {
    /// The center of the N-cube
    pub center: P,

    /// The width of the N-cube
    pub width: S,
}

impl<P, S> Ncube<P, S> {
    /// Create a new N-cube given its center and width
    pub fn new(center: P, width: S) -> Ncube<P, S> {
        Ncube { center: center, width: width }
    }
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
}


#[cfg(test)]
mod test {
    use super::*;
    use partition::Partition;
    use nalgebra::Pnt2;
    use quickcheck::quickcheck;

    #[test]
    fn ncube_total() {
        fn ncube_total(center: (f64, f64), width: f64, elem: (f64, f64)) -> bool {
            let center = Pnt2::new(center.0, center.1);
            let elem = Pnt2::new(elem.0, elem.1);
            Ncube::new(center, width).prop_is_total(&elem)
        }
        quickcheck(ncube_total as fn((f64, f64), f64, (f64, f64)) -> bool);
    }
}
