//! Pure tree implementation.

use std::mem;
use std::iter::IntoIterator;
use traits::{NodeState, Node, Position};
use partition::Partition;
use iter::Iter;


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
    pub fn new<I: Iterator<Item=O>>(iter: I, partition: P) -> PureTree<P, O> {
        let mut tree = PureTree::empty(partition);
        for object in iter {
            tree.insert(object)
        }
        tree
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
        use traits::NodeState::*;
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
    use rand::distributions::{IndependentSample, Range};
    use rand::thread_rng;
    use test::Bencher;
    use nalgebra::{Point2, FloatPoint, Origin};

    use traits::{Node, ObjectQuery, Positioned};
    use partition::Ncube;
    use super::*;

    #[bench]
    fn pure_tree_quad_new_1000(b: &mut Bencher) {
        let coord_distance = Range::new(-1.0f64, 1.0);
        let mut rng = thread_rng();
        let vec: Vec<_> = (0..1000).map(|_| Positioned {
            object: (),
            position: Point2::new(
                coord_distance.ind_sample(&mut rng),
                coord_distance.ind_sample(&mut rng)
            ),
        }).collect();
        b.iter(|| {
            PureTree::new(
                vec.iter().map(|a| a.clone()),
                Ncube::new(Origin::origin(), 2.0),
            )
        })
    }

    #[bench]
    fn pure_tree_query_objects(b: &mut Bencher) {
        let coord_distance = Range::new(-1.0f64, 1.0);
        let mut rng = thread_rng();
        let search_radius = 0.3;
        let tree = PureTree::new(
            (0..1000).map(|_| Positioned {
                object: (),
                position: Point2::new(
                    coord_distance.ind_sample(&mut rng),
                    coord_distance.ind_sample(&mut rng)
                ),
            }),
            Ncube::new(Origin::origin(), 200.0),
        );
        b.iter(|| {
            // Count the number of objects within the search radius 10000 times
            (0..10000)
                .map(|_|
                    tree.query_objects(|node|
                        node.partition()
                        .center().distance(&Origin::origin())
                        < search_radius + node.partition().width() / 2.0,
                    )
                    .filter(|other| other.position.distance(&Origin::origin()) < search_radius)
                    .count()
                )
                .fold(0, |a, b| a + b)
        })
    }
}
