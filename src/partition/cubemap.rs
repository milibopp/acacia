//! Cubemapping module

use nalgebra::{BaseFloat, Vec3, Vec2, Norm, zero, one};
#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::{Arbitrary, Gen};
use partition::{Partition, Subdivide, UnitQuad};


/// An axis direction
///
/// This effectively distinguishes whether we are moving in positive or negative
/// direction along some axis, i.e. +X vs -X, +Y vs. -Y etc.
#[derive(Copy, Clone, Debug, PartialEq, Hash, Eq)]
pub enum Direction {
    /// Positive direction
    Positive,

    /// Negative direction
    Negative,
}

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
#[derive(Copy, Clone, Debug, PartialEq, Hash, Eq)]
pub enum Axis {
    /// X-axis
    X,

    /// Y-axis
    Y,

    /// Z-axis
    Z,
}

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


/// Get the triple of axis vectors
///
/// The first vector is the normal vector n, the remaining are tangents t_1 and
/// t_2. They form a basis that is right-handed, i.e. n × t_1 = t_2.
pub fn axis_vector_triple<T: BaseFloat>(axis: Axis, direction: Direction) -> [Vec3<T>; 3] {
    let _p: T = one();
    let _n: T = -_p;
    let _0: T = zero();
    let sgn = match direction {
        Direction::Positive => _p,
        Direction::Negative => _n,
    };
    match axis {
        Axis::X => [
            Vec3::new(sgn, _0, _0),
            Vec3::new(_0, sgn, _0),
            Vec3::new(_0, _0, _p),
        ],
        Axis::Y => [
            Vec3::new(_0, sgn, _0),
            Vec3::new(_0, _0, sgn),
            Vec3::new(_p, _0, _0),
        ],
        Axis::Z => [
            Vec3::new(_0, _0, sgn),
            Vec3::new(sgn, _0, _0),
            Vec3::new(_0, _p, _0),
        ],
    }
}


/// A quad-shaped partition of the side of a cubemap
#[derive(Copy, Clone, Debug, PartialEq, Hash, Eq)]
pub struct Quad {
    /// Normal axis of the quad normal
    pub axis: Axis,

    /// Direction of the quad normal along the axis
    pub direction: Direction,

    /// Embedded flat unit quad
    pub flat_quad: UnitQuad,
}

impl Quad {
    /// The center of this quad on the cube
    pub fn center_on_cube<T: BaseFloat>(&self) -> Vec3<T> {
        let _1: T = one();
        let _2: T = _1 + _1;
        let c: Vec2<T> = self.flat_quad.center();
        let [n, t1, t2]: [Vec3<T>; 3] =
            axis_vector_triple(self.axis, self.direction);
        n + t1 * (c.x * _2 - _1) + t2 * (c.y * _2 - _1)
    }

    /// The center of this quad on the unit sphere
    pub fn center_on_sphere<T: BaseFloat>(&self) -> Vec3<T> {
        self.center_on_cube().normalize()
    }
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
        let _1: T = one();
        let _2: T = _1 + _1;
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


/// A cubemap partition of a 3-vector space
///
/// This has no radial partitioning, as it is intended mainly for the surface of
/// a 2-sphere (which is a subset of the full ℝ³). It is either the full
/// spherical dome or some subdivision stage on one of the six quad-shape sides
/// obtained by projecting the sphere onto a cube.
#[derive(Copy, Clone, Debug, PartialEq, Hash, Eq)]
pub enum CubeMap {
    /// The full sphere
    Sphere,

    /// A quad-based subdivision
    Quad(Quad),
}

impl Subdivide for CubeMap {
    fn subdivide(&self) -> Vec<CubeMap> {
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
                    flat_quad: UnitQuad::new(0, (0, 0)),
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
    pub use nalgebra::{Vec3, Cross};
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

    #[test]
    fn axis_vector_triples_are_right_handed() {
        fn check(axis: Axis, direction: Direction) -> bool {
            let [n, t1, t2]: [Vec3<f64>; 3] = axis_vector_triple(axis, direction);
            n.cross(&t1) == t2
        }
        quickcheck(check as fn(Axis, Direction) -> bool);
    }

    #[test]
    fn axis_vector_triples_concrete() {
        assert_eq!(
            axis_vector_triple::<f64>(Axis::X, Direction::Negative),
            [-Vec3::x(), -Vec3::y(), Vec3::z()]
        );
        assert_eq!(
            axis_vector_triple::<f64>(Axis::X, Direction::Positive),
            [Vec3::x(), Vec3::y(), Vec3::z()]
        );
    }
}
