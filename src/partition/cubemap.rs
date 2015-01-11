use std::num::cast;
use nalgebra::{BaseFloat, Vec3, Vec2, zero};
#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::{Arbitrary, Gen};
use partition::{Partition, Subdivide};
use partition::UnitQuad;

/// An axis direction
#[derive(Copy, Clone, Show)]
pub enum Direction { Positive, Negative }

#[cfg(any(test, feature = "arbitrary"))]
impl Arbitrary for Direction {
    fn arbitrary<G: Gen>(g: &mut G) -> Direction {
        match g.gen_range(0, 2) {
            0 => Direction::Positive,
            1 => Direction::Negative,
            _ => unreachable!(),
        }
    }
}


/// A coordinate axis
#[derive(Copy, Clone, Show)]
pub enum Axis { X, Y, Z }

#[cfg(any(test, feature = "arbitrary"))]
impl Arbitrary for Axis {
    fn arbitrary<G: Gen>(g: &mut G) -> Axis {
        match g.gen_range(0, 3) {
            0 => Axis::X,
            1 => Axis::Y,
            2 => Axis::Z,
            _ => unreachable!(),
        }
    }
}


#[derive(Copy, Clone, Show)]
pub struct Quad {
    axis: Axis,
    direction: Direction,
    flat_quad: UnitQuad,
}

impl Subdivide for Quad {
    fn subdivide(&self) -> Vec<Quad> {
        self.flat_quad.subdivide()
            .into_iter()
            .map(|q| Quad {
                axis: self.axis,
                direction: self.direction,
                flat_quad: q,
            })
            .collect()
    }
}

impl<T: BaseFloat + PartialOrd> Partition<Vec3<T>> for Quad {
    fn contains(&self, elem: &Vec3<T>) -> bool {
        let _1: T = cast(1.0).unwrap();
        let _2: T = cast(2.0).unwrap();
        let (i, j, k) = match self.axis {
            Axis::X => (0, 1, 2),
            Axis::Y => (1, 2, 0),
            Axis::Z => (2, 0, 1),
        };
        match (elem[i] > zero(), self.direction) {
            (true, Direction::Positive) | (false, Direction::Negative) =>
                self.flat_quad.contains(&Vec2::new(
                    (elem[j] / elem[i] + _1) / _2,
                    (elem[k] / elem[i] + _1) / _2,
                )),
            _ => false,
        }
    }
}

#[cfg(any(test, feature = "arbitrary"))]
impl Arbitrary for Quad {
    fn arbitrary<G: Gen>(g: &mut G) -> Quad {
        Quad {
            axis: Arbitrary::arbitrary(g),
            direction: Arbitrary::arbitrary(g),
            flat_quad: Arbitrary::arbitrary(g),
        }
    }
}


#[derive(Copy, Clone, Show)]
pub enum CubeMap {
    Sphere,
    Quad(Quad),
}

impl Subdivide for CubeMap {
    fn subdivide(&self) -> Vec<CubeMap> {
        println!("{:?}.subdivide()", *self);
        match *self {
            CubeMap::Sphere => 
                vec![
                    (Direction::Positive, Axis::X),
                    (Direction::Positive, Axis::Y),
                    (Direction::Positive, Axis::Z),
                    (Direction::Negative, Axis::X),
                    (Direction::Negative, Axis::Y),
                    (Direction::Negative, Axis::Z),
                ]
                .into_iter()
                .map(|(dir, ax)| CubeMap::Quad(Quad {
                    axis: ax,
                    direction: dir,
                    flat_quad: UnitQuad::new(0, [0, 0]),
                }))
                .collect(),
            CubeMap::Quad(ref quad) =>
                quad.subdivide().into_iter().map(|q| CubeMap::Quad(q)).collect(),
        }
    }
}

impl<T: BaseFloat + PartialOrd> Partition<Vec3<T>> for CubeMap {
    fn contains(&self, elem: &Vec3<T>) -> bool {
        match *self {
            CubeMap::Sphere => true,
            CubeMap::Quad(q) => q.contains(elem),
        }
    }
}

#[cfg(any(test, feature = "arbitrary"))]
impl Arbitrary for CubeMap {
    fn arbitrary<G: Gen>(g: &mut G) -> CubeMap {
        match { let s = g.size(); g.gen_range(0, s) } {
            0 => CubeMap::Sphere,
            _ => CubeMap::Quad(Arbitrary::arbitrary(g)),
        }
    }
}


#[cfg(test)]
mod test {
    pub use nalgebra::Vec3;
    pub use super::*;
    use quickcheck::quickcheck;
    use partition::Partition;

    partition_quickcheck!(quad_vec3_f32, Quad, Vec3<f32>);
    partition_quickcheck!(quad_vec3_f64, Quad, Vec3<f64>);
    partition_quickcheck!(cubemap_vec3_f32, CubeMap, Vec3<f32>);
    partition_quickcheck!(cubemap_vec3_f64, CubeMap, Vec3<f64>);

    #[test]
    fn cubemap_covers_vec3() {
        fn check(v: Vec3<f64>) -> bool {
            CubeMap::Sphere.contains(&v)
        }
        quickcheck(check as fn(Vec3<f64>) -> bool);
    }
}
