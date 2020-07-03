extern crate acacia;
extern crate nalgebra;
extern crate rand;

use acacia::partition::Ncube;
use acacia::{AssociatedData, DataQuery, Node, Position, Tree};
use nalgebra::{distance, zero, Point3, Vector3};
use rand::distributions::{Distribution, Uniform};
use rand::thread_rng;

/// Point mass
struct PointMass {
    mass: f64,
    position: Point3<f64>,
}

impl Position for PointMass {
    type Point = Point3<f64>;
    fn position(&self) -> Point3<f64> {
        self.position
    }
}

/// Newton's law of gravity for two point masses (with G = 1)
fn newton(m: f64, p1: Point3<f64>, p2: Point3<f64>) -> Vector3<f64> {
    let diff: Vector3<f64> = p1 - p2;
    let r = diff.norm();
    diff * m / r.powi(3)
}

fn main() {
    // Shortcut
    let origin = Point3::origin();

    // Generate a number of particles
    let mut rng = thread_rng();
    let coord_range = Uniform::from(-5.0..5.0);
    let particles: Vec<_> = (0..1000)
        .map(|_| PointMass {
            mass: 1.0,
            position: Point3::new(
                coord_range.sample(&mut rng),
                coord_range.sample(&mut rng),
                coord_range.sample(&mut rng),
            ),
        })
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
        &|&(com1, m1), &(com2, m2)| {
            if m1 + m2 > 0.0 {
                (com1 + (com2 - com1) * (m2 / (m1 + m2)), m1 + m2)
            } else {
                (origin, 0.0)
            }
        },
    )
    .expect("Couldn't construct tree");

    // This is the point, at which we want to know the gravitational
    // acceleration.
    let test_point = Point3::new(2.3, -4.1, 1.6);

    let theta = 0.5; // A bit arbitrary but this appears to work
    let tree_gravity: Vector3<f64> =
        // This is the recursion criterion. If a branch node passes this, the
        // query continues on its children.
        tree.query_data(|node| {
            let &(ref center_of_mass, _) = node.data();
            let d = distance(&test_point, center_of_mass);
            let delta = distance(&node.partition().center(), center_of_mass);
            d < node.partition().width() / theta + delta
        })
        // This collects our force term from each piece of associated data the
        // tree encounters during recursion.
        .map(|&(center_of_mass, mass)| newton(mass, center_of_mass, test_point))
        .fold(zero(), |a, b| a + b);

    // Calculate gravity exactly for comparison
    let exact_gravity: Vector3<f64> = particles
        .iter()
        .map(|particle| newton(particle.mass, particle.position, test_point))
        .fold(zero(), |a, b| a + b);

    // Print result of calculation
    println!(
        "Tree-computed gravity at {:?} is {:?} (exact: {:?})",
        test_point, tree_gravity, exact_gravity
    );
}
