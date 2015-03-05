use std::num::{Int, Float, cast};
#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::{Arbitrary, Gen};
use nalgebra::{Vec2, BaseFloat};
use partition::{Partition, Subdivide};


/// A partition of the unit quad [0, 1) × [0, 1)
#[derive(Copy, Clone, Debug, PartialEq, Hash, Eq)]
pub struct UnitQuad {
    scale: u8,
    offset: (u32, u32),
}

impl UnitQuad {
    /// Create a new `UnitQuad`
    ///
    /// This asserts that the offset is valid given the scale level.
    pub fn new(scale: u8, offset: (u32, u32)) -> UnitQuad {
        assert!(scale < 32); // Otherwise exponentiation will overflow
        let max_offset = 2.pow(scale as u32);
        assert!(offset.0 < max_offset && offset.1 < max_offset);
        UnitQuad { scale: scale, offset: offset }
    }

    /// Integer scale
    pub fn scale(&self) -> u8 { self.scale }

    /// Integer offset
    pub fn offset(&self) -> (u32, u32) { self.offset }

    /// Get coordinate within the partition from (u, v) coordinates
    pub fn coordinate<T: BaseFloat>(&self, coord: (T, T)) -> Vec2<T> {
        let (u, v) = coord;
        let width: T = self.width();
        Vec2::new(
            width * (cast::<_, T>(self.offset.0).unwrap() + u),
            width * (cast::<_, T>(self.offset.1).unwrap() + v),
        )
    }

    /// Center of the partitioned region
    pub fn center<T: BaseFloat>(&self) -> Vec2<T> {
        let half = cast(0.5).unwrap();
        self.coordinate((half, half))
    }

    /// Width of the partitioned region
    pub fn width<T: BaseFloat>(&self) -> T {
        cast::<_, T>(0.5).unwrap().powi(self.scale as i32)
    }
}

impl Subdivide for UnitQuad {
    fn subdivide(&self) -> Vec<UnitQuad> {
        [(0, 0), (0, 1), (1, 0), (1, 1)]
            .iter()
            .map(|&(di, dj)| {
                let (i, j) = self.offset;
                UnitQuad::new(self.scale + 1, (i * 2 + di, j * 2 + dj))
            })
            .collect()
    }
}

impl<T: BaseFloat> Partition<Vec2<T>> for UnitQuad
{
    fn contains(&self, elem: &Vec2<T>) -> bool {
        let width: T = self.width();
        let offset = Vec2::new(
            width * cast(self.offset.0).unwrap(),
            width * cast(self.offset.1).unwrap(),
        );
        (offset.x < elem.x) && (elem.x < offset.x + width) &&
        (offset.y < elem.y) && (elem.y < offset.y + width)
    }
}

#[cfg(any(test, feature = "arbitrary"))]
impl Arbitrary for UnitQuad {
    fn arbitrary<G: Gen>(g: &mut G) -> UnitQuad {
        use std::cmp;
        // FIXME: somehow this explicit usage of the `Rng` trait is required for
        // the `.gen_range` calls, even though `Rng: Gen`. The compiler
        // complains that the "source trait is private". Curiously, adding this
        // import here fixes the same situation in the `cubemap` module as well.
        use rand::Rng;
        let scale: u8 = {
            // scale >= 32 is invalid (overflow)
            // At scale >= 31 subdivision fails
            let max_scale = cmp::min(31, g.size()) as u8;
            g.gen_range(0, max_scale)
        };
        let max_offset = 2.pow(scale as u32);
        UnitQuad::new(scale, (
            g.gen_range(0, max_offset),
            g.gen_range(0, max_offset),
        ))
    }
}


#[cfg(test)]
mod test {
    pub use nalgebra::Vec2;
    pub use super::*;
    use quickcheck::{quickcheck, TestResult};
    use partition::Partition;

    partition_quickcheck!(unitquad_vec2_f32, UnitQuad, Vec2<f32>);
    partition_quickcheck!(unitquad_vec2_f64, UnitQuad, Vec2<f64>);

    #[test]
    fn unitquad_base_contains_region() {
        fn check(v: Vec2<f64>) -> TestResult {
            if v.x < 0.0 || v.x >= 1.0 || v.y < 0.0 || v.y >= 1.0 {
                TestResult::discard()
            } else {
                TestResult::from_bool(UnitQuad::new(0, (0, 0)).contains(&v))
            }
        }
        quickcheck(check as fn(Vec2<f64>) -> TestResult);
    }
}
