//! Simple integration tests oriented towards gravity computations

extern crate acacia;
extern crate approx;
extern crate nalgebra;
extern crate quickcheck;

use acacia::partition::Ncube;
use acacia::{AssociatedData, DataQuery, Node, Positioned, Tree};
use approx::Relative;
use nalgebra::{distance, zero, Point2, Point3, Vector2, Vector3};
use quickcheck::{quickcheck, TestResult};
use std::fmt::Debug;

#[test]
fn tree_center_of_mass() {
    fn tree_center_of_mass(data: Vec<(f64, (f64, f64))>) -> TestResult {
        // Only test non-empty lists with positive masses
        if data.is_empty() || data.iter().any(|&(m, _)| m <= 0.0) {
            return TestResult::discard();
        }
        // No two points should be in the same place
        for &(_, pi) in &data {
            for &(_, pj) in &data {
                if pi == pj {
                    return TestResult::discard();
                }
            }
        }
        // Compute center of mass in the traditional way
        let (mps, ms) = data
            .iter()
            .map(|&(m, (x, y))| (Vector2::new(x, y) * m, m))
            .fold((zero::<Vector2<f64>>(), 0.0f64), |(mps, ms), (mp, m)| {
                (mps + mp, ms + m)
            });
        let com = mps / ms;
        // Now use the tree
        let tree = Tree::new(
            data.iter().map(|&(m, (x, y))| Positioned {
                object: m,
                position: Point2::new(x, y),
            }),
            Ncube::new(Point2::origin(), 200.0f64),
            (zero(), 0.0),
            &|obj| ((obj.position - Point2::origin()) * obj.object, obj.object),
            &|&(mps, ms), &(mp, m)| (mps + mp, ms + m),
        )
        .expect("Couldn't construct tree");
        let (tree_mps, tree_ms) = *tree.data();
        // …and compare
        TestResult::from_bool(Relative::default().eq(&(tree_mps / tree_ms), &com))
    }
    quickcheck(tree_center_of_mass as fn(Vec<(f64, (f64, f64))>) -> TestResult);
}

#[test]
fn tree_gravity_approx() {
    fn tree_gravity_approx(
        starfield: Vec<(f64, (f64, f64, f64))>,
        test_point: (f64, f64, f64),
    ) -> TestResult {
        // We want to have at least one star
        if starfield.is_empty() {
            return TestResult::discard();
        }
        // Only test positive masses
        if starfield.iter().any(|&(m, _)| m <= 0.0) {
            return TestResult::discard();
        }
        // The test point should not be in the same place as any star
        if starfield.iter().any(|&(_, p)| p == test_point) {
            return TestResult::discard();
        }
        // No two stars should be in the same place
        for &(_, pi) in &starfield {
            for &(_, pj) in &starfield {
                if pi == pj {
                    return TestResult::discard();
                }
            }
        }
        // (T, T, T) -> Point3<T>
        fn pnt<T: Debug + PartialEq + Copy + 'static>(p: (T, T, T)) -> Point3<T> {
            let (x, y, z) = p;
            Point3::new(x, y, z)
        }
        let test_point = pnt(test_point);
        // Newton's law of gravity for two point masses (with G = 1)
        let newton = |(m, p1): (f64, Point3<f64>), p2| {
            let diff: Vector3<f64> = p1 - p2;
            let r = diff.norm();
            diff * (m / r.powi(3))
        };
        // Calculate gravity exactly
        let simple_gravity = starfield
            .iter()
            .map(|&(m, p)| newton((m, pnt(p)), test_point))
            .fold(zero(), |a: Vector3<_>, b| a + b);
        // Calculate gravity using a tree
        let orig: Point3<f64> = Point3::origin();
        let data_width = (test_point - orig).norm() * 2.0;
        let width = if data_width < 200.0 {
            200.0
        } else {
            data_width
        };
        let tree = Tree::new(
            starfield.iter().map(|&(m, (x, y, z))| Positioned {
                object: m,
                position: Point3::new(x, y, z),
            }),
            Ncube::new(orig, width),
            (orig, zero()),
            &|obj| (obj.position, obj.object),
            &|&(com1, m1), &(com2, m2)| {
                if m1 + m2 > zero() {
                    (com1 + (com2 - com1) * (m2 / (m1 + m2)), m1 + m2)
                } else {
                    (orig, zero())
                }
            },
        )
        .expect("Couldn't construct tree");
        let theta = 0.5; // A bit arbitrary but this appears to work
        let tree_gravity = tree
            .query_data(|node| {
                let &(ref center_of_mass, _) = node.data();
                let d = distance(&test_point, center_of_mass);
                let delta = distance(&node.partition().center(), center_of_mass);
                d < node.partition().width() / theta + delta
            })
            .map(|&(com, m)| newton((m, com), test_point))
            .fold(zero::<Vector3<f64>>(), |a, b| a + b);

        // Now the tree gravity should approximate the exact one, within 10 %
        TestResult::from_bool(
            Relative::default()
                .epsilon(0.1 * simple_gravity.norm())
                .eq(&simple_gravity, &tree_gravity),
        )
    }
    quickcheck(
        tree_gravity_approx as fn(Vec<(f64, (f64, f64, f64))>, (f64, f64, f64)) -> TestResult,
    )
}
