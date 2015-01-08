//! Simple integration tests oriented towards gravity computations

extern crate acacia;
extern crate nalgebra;
extern crate quickcheck;

use std::num::Float;
use nalgebra::{ApproxEq, Pnt2, Pnt3, FloatPnt, Vec2, Vec3, zero, Norm, Orig};
use quickcheck::{TestResult, quickcheck};
use acacia::tree::{Node, AssociatedData, DataQuery, Positioned};
use acacia::ntree::NTree;


#[test]
fn ntree_center_of_mass() {
    fn ntree_center_of_mass(data: Vec<(f64, (f64, f64))>) -> TestResult {
        // Only test non-empty lists with positive masses
        if data.is_empty() || data.iter().any(|&(m, _)| m <= 0.0) {
            return TestResult::discard();
        }
        // No two points should be in the same place
        for i in range(0, data.len()) {
            for j in range(0, i) {
                let (_, pi) = data[i];
                let (_, pj) = data[j];
                if pi == pj {
                    return TestResult::discard();
                }
            }
        }
        // Compute center of mass in the traditional way
        let (mps, ms) = data.iter()
            .map(|&(m, (x, y))| (Vec2::new(x, y) * m, m))
            .fold((zero::<Vec2<f64>>(), 0.0f64), |(mps, ms), (mp, m)| (mps + mp, ms + m));
        let com = mps / ms;
        // Now use the tree
        let tree = NTree::from_iter(
            data.iter().map(|&(m, (x, y))|
                Positioned { object: m, position: Pnt2::new(x, y) }
            ),
            (Vec2::new(0.0f64, 0.0), 0.0f64),
            &|obj| (obj.position.to_vec() * obj.object, obj.object),
            &|&(mps, ms), &(mp, m)| (mps + mp, ms + m)
        );
        let (tree_mps, tree_ms) = *tree.data();
        // …and compare
        TestResult::from_bool(ApproxEq::approx_eq(&(tree_mps / tree_ms), &com))
    }
    quickcheck(ntree_center_of_mass as fn(Vec<(f64, (f64, f64))>) -> TestResult);
}

#[test]
fn ntree_gravity_approx() {
    fn ntree_gravity_approx(
            starfield: Vec<(f64, (f64, f64, f64))>,
            test_point_: (f64, f64, f64)
        ) -> TestResult
    {
        // We want to have at least one star
        if starfield.is_empty() {
            return TestResult::discard();
        }
        // Only test positive masses
        if starfield.iter().any(|&(m, _)| m <= 0.0) {
            return TestResult::discard();
        }
        // The test point should not be in the same place as any star
        if starfield.iter().any(|&(_, p)| p == test_point_) {
            return TestResult::discard();
        }
        // No two stars should be in the same place
        for i in range(0, starfield.len()) {
            for j in range(0, i) {
                let (_, pi) = starfield[i];
                let (_, pj) = starfield[j];
                if pi == pj {
                    return TestResult::discard();
                }
            }
        }
        // (T, T, T) -> Pnt3<T>
        let pnt = |&:p| {
            let (x, y, z) = p;
            Pnt3::new(x, y, z)
        };
        let test_point = pnt(test_point_);
        // Newton's law of gravity for two point masses (with G = 1)
        let newton = |&:(m, p1): (f64, Pnt3<f64>), p2| {
            let diff: Vec3<f64> = p1 - p2;
            let r = diff.norm();
            diff * (m / r.powi(3))
        };
        // Calculate gravity exactly
        let simple_gravity = starfield.iter()
            .map(|&(m, p)| newton((m, pnt(p)), test_point))
            .fold(zero(), |a: Vec3<_>, b| a + b);
        // Calculate gravity using a tree
        let orig: Pnt3<f64> = Orig::orig();
        let tree = NTree::from_iter_with_geometry(
            starfield.iter().map(|&(m, (x, y, z))|
                Positioned { object: m, position: Pnt3::new(x, y, z) }
            ),
            orig,
            test_point.as_vec().norm() * 2.0,
            (orig, zero()),
            &|obj| (obj.position, obj.object),
            &|&(com1, m1), &(com2, m2)|
                if m1 + m2 > zero() {(
                    orig + (com1.to_vec() * m1 + com2.to_vec() * m2) / (m1 + m2),
                    m1 + m2,
                )}
                else {
                    (orig, zero())
                }
        );
        let theta = 0.5; // A bit arbitrary but this appears to work
        let mut tree_gravity: Vec3<_> = zero();
        tree.query_data(
            &|node| {
                let &(ref center_of_mass, _) = node.data();
                let d = FloatPnt::dist(&test_point, center_of_mass);
                let delta = FloatPnt::dist(node.center(), center_of_mass);
                d < *node.width() / theta + delta
            },
            &mut |&(com, m)| {
                tree_gravity = tree_gravity + newton((m, com), test_point);
            },
        );
        // Now the tree gravity should approximate the exact one, within 5 %
        TestResult::from_bool(simple_gravity.approx_eq_eps(&tree_gravity, &(0.05 * simple_gravity.norm())))
    }
    quickcheck(ntree_gravity_approx as fn(Vec<(f64, (f64, f64, f64))>, (f64, f64, f64)) -> TestResult)
}
