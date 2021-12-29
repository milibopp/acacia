//! Pure tree implementation.

use std::{
    mem,
    iter::IntoIterator,
};
use crate::{
    traits::{NodeState, Node, Position},
    partition::Partition,
    iter::Iter,
    error::ConstructionError,
};


/// A pure N-dimensional tree
pub struct PureTree<P, O> {
    state: NodeState<O, Vec<PureTree<P, O>>>,
    partition: P,
}

impl<P, O> PureTree<P, O>
    where P: Partition<<O as Position>::Point>,
          O: Position,
{
    fn empty(partition: P) -> PureTree<P, O> {
        PureTree {
            state: NodeState::Empty,
            partition: partition,
        }
    }
}

impl<P, O> PureTree<P, O>
    where O: Position,
          P: Partition<<O as Position>::Point>,
{
    /// Construct a tree without checking the geometry of the input data
    pub fn new<I: Iterator<Item=O>>(iter: I, partition: P) -> Result<PureTree<P, O>, ConstructionError> {
        let mut tree = PureTree::empty(partition);
        for object in iter {
            if tree.partition.contains(&object.position()) {
                tree.insert(object)
            } else {
                return Err(ConstructionError::ObjectOutsidePartition);
            }
        }
        Ok(tree)
    }

    fn dispatch(&self, nodes: &mut Vec<PureTree<P, O>>, object: O) {
        nodes[self.partition.dispatch(&object.position())].insert(object)
    }

    fn insert(&mut self, object: O) {
        let mut tmp = NodeState::Empty;
        mem::swap(&mut tmp, &mut self.state);
        self.state = match tmp {
            NodeState::Empty => NodeState::Leaf(object),
            NodeState::Leaf(other) => {
                let mut nodes: Vec<_> = self.partition.subdivide()
                    .into_iter()
                    .map(|p| PureTree::empty(p))
                    .collect();
                self.dispatch(&mut nodes, object);
                self.dispatch(&mut nodes, other);
                NodeState::Branch(nodes)
            },
            NodeState::Branch(mut nodes) => {
                self.dispatch(&mut nodes, object);
                NodeState::Branch(nodes)
            }
        };
    }
}

impl<P: Clone, O> Node for PureTree<P, O> {
    type Partition = P;
    type Object = O;
    type Container = Vec<PureTree<P, O>>;

    fn state(&self) -> NodeState<&O, &Vec<PureTree<P, O>>> {
        use crate::traits::NodeState::*;
        match self.state {
            Empty => Empty,
            Leaf(ref obj) => Leaf(obj),
            Branch(ref vec) => Branch(vec),
        }
    }

    fn partition(&self) -> P {
        self.partition.clone()
    }
}

impl<'a, P: Clone + 'a, O: 'a> IntoIterator for &'a PureTree<P, O> {
    type Item = &'a O;
    type IntoIter = Iter<'a, PureTree<P, O>>;
    fn into_iter(self) -> Iter<'a, PureTree<P, O>> { Iter::new(self) }
}


#[cfg(test)]
mod test {
    // use test::Bencher;
    use nalgebra::Point2;
    use quickcheck::{TestResult, quickcheck};

    use crate::{
        error::ConstructionError,
        traits::Positioned,
        partition::Ncube,
    };
    use super::*;

    #[test]
    fn pure_tree_construction_error_object_outside_partition() {
        fn pure_tree_construction_error_object_outside_partition(input: (Vec<(f64, f64)>, f64)) -> TestResult {
            let (data, domain_size) = input;

            // Have atleast one object, the size of the domain should be positive, and atleast
            // one object outside the domain should exist
            if data.is_empty() ||
               domain_size <= 0.0 ||
               data.iter().all(|&(x, y)| x.abs() < domain_size && y.abs() < domain_size)
            {
                return TestResult::discard();
            }

            TestResult::from_bool(match PureTree::new(
                data.iter().map(|&(x, y)| Positioned { object: (), position: Point2::new(x, y) }),
                Ncube::new(Point2::origin(), domain_size)
            ) {
                Err(ConstructionError::ObjectOutsidePartition) => true,
                _ => false
            })
        }
        quickcheck(pure_tree_construction_error_object_outside_partition as fn(input: (Vec<(f64, f64)>, f64)) -> TestResult)
    }

    // #[bench]
    // fn pure_tree_quad_new_1000(b: &mut Bencher) {
    //     let coord_distance = Range::new(-1.0f64, 1.0);
    //     let mut rng = thread_rng();
    //     let vec: Vec<_> = (0..1000).map(|_| Positioned {
    //         object: (),
    //         position: Point2::new(
    //             coord_distance.ind_sample(&mut rng),
    //             coord_distance.ind_sample(&mut rng)
    //         ),
    //     }).collect();
    //     b.iter(|| {
    //         PureTree::new(
    //             vec.iter().map(|a| a.clone()),
    //             Ncube::new(Point2::origin(), 2.0),
    //         )
    //     })
    // }

    // #[bench]
    // fn pure_tree_query_objects(b: &mut Bencher) {
    //     let coord_distance = Range::new(-1.0f64, 1.0);
    //     let mut rng = thread_rng();
    //     let search_radius = 0.3;
    //     let tree = PureTree::new(
    //         (0..1000).map(|_| Positioned {
    //             object: (),
    //             position: Point2::new(
    //                 coord_distance.ind_sample(&mut rng),
    //                 coord_distance.ind_sample(&mut rng)
    //             ),
    //         }),
    //         Ncube::new(Point2::origin(), 200.0),
    //     ).expect("Couldn't construct tree");
    //     b.iter(|| {
    //         // Count the number of objects within the search radius 10000 times
    //         (0..10000)
    //             .map(|_|
    //                 tree.query_objects(|node|
    //                     node.partition()
    //                     .center().distance(&Point2::origin())
    //                     < search_radius + node.partition().width() / 2.0,
    //                 )
    //                 .filter(|other| other.position.distance(&Point2::origin()) < search_radius)
    //                 .count()
    //             )
    //             .fold(0, |a, b| a + b)
    //     })
    // }
}
