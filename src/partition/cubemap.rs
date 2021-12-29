//! Cubemapping module

use num_traits::{NumCast, Float};
use nalgebra::{Scalar, Vector3, Vector2, zero, one, ComplexField, ClosedMul, ClosedAdd};
#[cfg(any(test, feature = "arbitrary"))]
use quickcheck::{Arbitrary, Gen};
use super::{Partition, Subdivide, UnitQuad};


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
    fn arbitrary(g: &mut Gen) -> Direction {
        *g.choose(&[Direction::Positive, Direction::Negative]).unwrap()
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
    fn arbitrary(g: &mut Gen) -> Axis {
        *g.choose(&[Axis::X, Axis::Y, Axis::Z]).unwrap()
    }
}


/// Get the triple of axis vectors
///
/// The first vector is the normal vector n, the remaining are tangents t_1 and
/// t_2. They form a basis that is right-handed, i.e. n × t_1 = t_2.
pub fn axis_vector_triple<T: Scalar + Float>(axis: Axis, direction: Direction) -> [Vector3<T>; 3] {
    let _p: T = one();
    let _n: T = -_p;
    let _0: T = zero();
    let sgn = match direction {
        Direction::Positive => _p,
        Direction::Negative => _n,
    };
    match axis {
        Axis::X => [
            Vector3::new(sgn, _0, _0),
            Vector3::new(_0, sgn, _0),
            Vector3::new(_0, _0, _p),
        ],
        Axis::Y => [
            Vector3::new(_0, sgn, _0),
            Vector3::new(_0, _0, sgn),
            Vector3::new(_p, _0, _0),
        ],
        Axis::Z => [
            Vector3::new(_0, _0, sgn),
            Vector3::new(sgn, _0, _0),
            Vector3::new(_0, _p, _0),
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
    pub fn center_on_cube<T: Scalar + NumCast + Float + ClosedMul + ClosedAdd>(&self) -> Vector3<T> {
        let _1: T = one();
        let _2: T = _1 + _1;
        let c: Vector2<T> = self.flat_quad.center();
        let triple: [Vector3<T>; 3] =
            axis_vector_triple(self.axis, self.direction);
        let n: Vector3<T> = triple[0];
        let t1 = triple[1];
        let t2 = triple[2];
        n + t1 * (c.x * _2 - _1) + t2 * (c.y * _2 - _1)
    }

    /// The center of this quad on the unit sphere
    pub fn center_on_sphere<T: Scalar + NumCast + Float + ComplexField>(&self) -> Vector3<T> {
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

impl<T: Scalar + PartialOrd + NumCast + Float> Partition<Vector3<T>> for Quad {
    fn contains(&self, elem: &Vector3<T>) -> bool {
        let _1: T = one();
        let _2: T = _1 + _1;
        let (i, j, k) = match self.axis {
            Axis::X => (0, 1, 2),
            Axis::Y => (1, 2, 0),
            Axis::Z => (2, 0, 1),
        };
        match (elem[i] > zero(), self.direction) {
            (true, Direction::Positive) | (false, Direction::Negative) =>
                self.flat_quad.contains(&Vector2::new(
                    (elem[j] / elem[i] + _1) / _2,
                    (elem[k] / elem[i] + _1) / _2,
                )),
            _ => false,
        }
    }
}

#[cfg(any(test, feature = "arbitrary"))]
impl Arbitrary for Quad {
    fn arbitrary(g: &mut Gen) -> Quad {
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

impl<T: Scalar + PartialOrd + NumCast + Float> Partition<Vector3<T>> for CubeMap {
    fn contains(&self, elem: &Vector3<T>) -> bool {
        match *self {
            CubeMap::Sphere => true,
            CubeMap::Quad(q) => q.contains(elem),
        }
    }
}

#[cfg(any(test, feature = "arbitrary"))]
impl Arbitrary for CubeMap {
    fn arbitrary(g: &mut Gen) -> CubeMap {
        if bool::arbitrary(g) {
            CubeMap::Sphere
        } else {
            CubeMap::Quad(Arbitrary::arbitrary(g))
        }
    }
}


#[cfg(test)]
mod test {
    use nalgebra::Vector3;
    use super::*;
    use quickcheck::quickcheck;

    partition_quickcheck!(quad_vec3_f32, Quad, Vector3<f32>);
    partition_quickcheck!(quad_vec3_f64, Quad, Vector3<f64>);
    partition_quickcheck!(cubemap_vec3_f32, CubeMap, Vector3<f32>);
    partition_quickcheck!(cubemap_vec3_f64, CubeMap, Vector3<f64>);

    #[test]
    fn cubemap_covers_vec3() {
        fn check(v: Vector3<f64>) -> bool {
            CubeMap::Sphere.contains(&v)
        }
        quickcheck(check as fn(Vector3<f64>) -> bool);
    }

    #[test]
    fn axis_vector_triples_are_right_handed() {
        fn check(axis: Axis, direction: Direction) -> bool {
            let triple: [Vector3<f64>; 3] =
                axis_vector_triple(axis, direction);
            let n = triple[0];
            let t1 = triple[1];
            let t2 = triple[2];
            n.cross(&t1) == t2
        }
        quickcheck(check as fn(Axis, Direction) -> bool);
    }

    #[test]
    fn axis_vector_triples_concrete() {
        assert_eq!(
            axis_vector_triple::<f64>(Axis::X, Direction::Negative),
            [-Vector3::x(), -Vector3::y(), Vector3::z()]
        );
        assert_eq!(
            axis_vector_triple::<f64>(Axis::X, Direction::Positive),
            [Vector3::x(), Vector3::y(), Vector3::z()]
        );
    }
}
