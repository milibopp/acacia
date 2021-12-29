use num_traits::{PrimInt, Float, NumCast};
#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::{Arbitrary, Gen};
use nalgebra::{Vector2, Scalar};
use super::{Partition, Subdivide};


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
    pub fn coordinate<T: Scalar + NumCast + Float>(&self, coord: (T, T)) -> Vector2<T> {
        let (u, v) = coord;
        let width: T = self.width();
        Vector2::new(
            width * (<T as NumCast>::from(self.offset.0).unwrap() + u),
            width * (<T as NumCast>::from(self.offset.1).unwrap() + v),
        )
    }

    /// Center of the partitioned region
    pub fn center<T: Scalar + NumCast + Float>(&self) -> Vector2<T> {
        let half = NumCast::from(0.5).unwrap();
        self.coordinate((half, half))
    }

    /// Width of the partitioned region
    pub fn width<T: Scalar + NumCast + Float>(&self) -> T {
        Float::powi(<T as NumCast>::from(0.5).unwrap(), self.scale as i32)
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

impl<T: Scalar + NumCast + Float> Partition<Vector2<T>> for UnitQuad
{
    fn contains(&self, elem: &Vector2<T>) -> bool {
        let width: T = self.width();
        let offset = Vector2::new(
            width * NumCast::from(self.offset.0).unwrap(),
            width * NumCast::from(self.offset.1).unwrap(),
        );
        (offset.x < elem.x) && (elem.x < offset.x + width) &&
        (offset.y < elem.y) && (elem.y < offset.y + width)
    }
}

#[cfg(any(test, feature = "arbitrary"))]
impl Arbitrary for UnitQuad {
    fn arbitrary(g: &mut Gen) -> UnitQuad {
        let scale: u8 = {
            // scale >= 32 is invalid (overflow)
            // At scale >= 31 subdivision fails
            u8::arbitrary(g) % 31
        };
        let mask = 2.pow(scale as u32) - 1;
        UnitQuad::new(scale, (
            u32::arbitrary(g) & mask,
            u32::arbitrary(g) & mask,
        ))
    }
}


#[cfg(test)]
mod test {
    use nalgebra::Vector2;
    use super::*;
    use quickcheck::{quickcheck, TestResult};

    partition_quickcheck!(unitquad_vec2_f32, UnitQuad, Vector2<f32>);
    partition_quickcheck!(unitquad_vec2_f64, UnitQuad, Vector2<f64>);

    #[test]
    fn unitquad_base_contains_region() {
        fn check(v: Vector2<f64>) -> TestResult {
            if v.x < 0.0 || v.x >= 1.0 || v.y < 0.0 || v.y >= 1.0 {
                TestResult::discard()
            } else {
                TestResult::from_bool(UnitQuad::new(0, (0, 0)).contains(&v))
            }
        }
        quickcheck(check as fn(Vector2<f64>) -> TestResult);
    }
}
