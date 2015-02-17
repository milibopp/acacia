extern crate rand;
extern crate nalgebra;
extern crate acacia;

use rand::thread_rng;
use rand::distributions::{Range, IndependentSample};
use std::num::Float;
use nalgebra::{Pnt3, Vec3, FloatPnt, Norm, Orig, zero};
use acacia::{Tree, DataQuery, Node, AssociatedData, Position};
use acacia::partition::Ncube;

/// Point mass
struct PointMass {
    mass: f64,
    position: Pnt3<f64>,
}

impl Position for PointMass {
    type Point = Pnt3<f64>;
    fn position(&self) -> Pnt3<f64> {
        self.position
    }
}


/// Newton's law of gravity for two point masses (with G = 1)
fn newton(m: f64, p1: Pnt3<f64>, p2: Pnt3<f64>) -> Vec3<f64> {
    let diff: Vec3<f64> = p1 - p2;
    let r = diff.norm();
    diff * m / r.powi(3)
}


fn main() {
    // Shortcut
    let origin: Pnt3<f64> = Orig::orig();

    // Generate a number of particles
    let mut rng = thread_rng();
    let coord_range = Range::new(-5.0, 5.0);
    let particles: Vec<_> = (0..1000).map(|_|
        PointMass { mass: 1.0, position: Pnt3::new(
            coord_range.ind_sample(&mut rng),
            coord_range.ind_sample(&mut rng),
            coord_range.ind_sample(&mut rng)) })
        .collect();

    // Construct the tree
    let tree = Tree::new(

        // Pass in an iterator over the objects to store in the tree
        // In this case we pass in &PointMass, which implements Positionable.
        particles.iter(),

        // Shape of the root node
        Ncube::new(origin, 11.0),

        // The value for the associated data of empty nodes. Here, we associate
        // a center of mass and a total mass to each node.
        (origin, 0.0),

        // This closure associates data to a leaf node from its object.
        &|obj| (obj.position, obj.mass),

        // This combines two pieces of associated data and thus prescribes how
        // branch nodes on higher levels get their associated data.
        &|&(com1, m1), &(com2, m2)|
            if m1 + m2 > 0.0 {(
                origin + (com1.to_vec() * m1 + com2.to_vec() * m2) / (m1 + m2),
                m1 + m2,
            )}
            else {
                (origin, 0.0)
            }
    );

    // This is the point, at which we want to know the gravitational
    // acceleration.
    let test_point = Pnt3::new(2.3, -4.1, 1.6);

    let theta = 0.5; // A bit arbitrary but this appears to work
    let tree_gravity: Vec3<f64> =
        // This is the recursion criterion. If a branch node passes this, the
        // query continues on its children.
        tree.query_data(|node| {
            let &(ref center_of_mass, _) = node.data();
            let d = FloatPnt::dist(&test_point, center_of_mass);
            let delta = FloatPnt::dist(&node.partition().center(), center_of_mass);
            d < node.partition().width() / theta + delta
        })
        // This collects our force term from each piece of associated data the
        // tree encounters during recursion.
        .map(|&(center_of_mass, mass)| newton(mass, center_of_mass, test_point))
        .fold(zero(), |a, b| a + b);

    // Calculate gravity exactly for comparison
    let exact_gravity: Vec3<f64> = particles.iter()
        .map(|particle| newton(particle.mass, particle.position, test_point))
        .fold(zero(), |a, b| a + b);

    // Print result of calculation
    println!("Tree-computed gravity at {:?} is {:?} (exact: {:?})",
        test_point, tree_gravity, exact_gravity
    );
}
