//! Implementation of tree with associated data.

use std::mem;
use std::iter::IntoIterator;
use traits::{NodeState, DataQuery, ObjectQuery, AssociatedData, Node, Position};
use partition::Partition;
use iter::Iter;


/// An N-dimensional tree
///
/// This tree does not know the dimension of its point at compile time, as it is
/// not hard-coded and genericity over constants is unsupported in Rust.
pub struct Tree<P, O, D> {
    state: NodeState<O, Vec<Tree<P, O, D>>>,
    partition: P,
    data: D,
}

impl<P, O, D> Tree<P, O, D> {
    /// Construct an empty tree
    fn empty(partition: P, data: D) -> Tree<P, O, D> {
        Tree {
            state: NodeState::Empty,
            partition: partition,
            data: data,
        }
    }
}

impl<P, O, D: Clone> Tree<P, O, D> {
    /// Recompute the associated data
    fn recompute_data<S, C>(&mut self, default: D, single: &S, combine: &C)
        where S: Fn(&O) -> D,
              C: Fn(&D, &D) -> D,
    {
        self.data = match self.state {
            NodeState::Empty => default,
            NodeState::Leaf(ref obj) => single(obj),
            NodeState::Branch(ref mut nodes) => {
                for node in nodes.iter_mut() {
                    node.recompute_data(default.clone(), single, combine);
                }
                nodes.iter().fold(default.clone(), |current, node| combine(&current, &node.data))
            },
        };
    }
}

impl<P, O, D> Tree<P,  O, D>
    where O: Position,
          P: Partition<<O as Position>::Point>,
          D: Clone,
{
    fn dispatch(&self, nodes: &mut Vec<Tree<P, O, D>>, object: O, default: D) {
        nodes[self.partition.dispatch(&object.position())].insert(object, default)
    }

    /// Inserts a new object into the tree
    ///
    /// NOTE: this does not update the data correctly but merely places a
    /// default value in there.
    fn insert(&mut self, object: O, default: D) {
        let mut tmp = NodeState::Empty;
        mem::swap(&mut tmp, &mut self.state);
        self.state = match tmp {
            NodeState::Empty => NodeState::Leaf(object),
            NodeState::Leaf(other) => {
                let mut nodes: Vec<_> = self.partition.subdivide()
                    .into_iter()
                    .map(|p| Tree::empty(p, default.clone()))
                    .collect();
                self.dispatch(&mut nodes, object, default.clone());
                self.dispatch(&mut nodes, other, default.clone());
                NodeState::Branch(nodes)
            },
            NodeState::Branch(mut nodes) => {
                self.dispatch(&mut nodes, object, default.clone());
                NodeState::Branch(nodes)
            },
        };
    }

    /// Construct the tree from an iterator
    pub fn new<I, S, C>(objects: I, partition: P, default: D, single: &S, combine: &C) -> Tree<P, O, D>
        where I: Iterator<Item=O>,
              S: Fn(&O) -> D,
              C: Fn(&D, &D) -> D,
    {
        let mut tree = Tree::empty(partition, default.clone());
        for object in objects {
            tree.insert(object, default.clone());
        }
        tree.recompute_data(default.clone(), single, combine);
        tree
    }
}


impl<P: Clone, O, D> Node for Tree<P, O, D> {
    type Partition = P;
    type Object = O;
    type Container = Vec<Tree<P, O, D>>;

    fn state(&self) -> NodeState<&O, &Vec<Tree<P, O, D>>> {
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

impl<P, O, D> AssociatedData for Tree<P, O, D> {
    type Data = D;

    fn data(&self) -> &D {
        &self.data
    }
}

impl<P, O, D> DataQuery for Tree<P, O, D> {
    fn query_data<R, F>(&self, recurse: &R, f: &mut F)
        where R: Fn(&Tree<P, O, D>) -> bool,
              F: FnMut(&D),
    {
        match self.state {
            NodeState::Branch(ref nodes) if recurse(self) =>
                for node in nodes.iter() {
                    node.query_data(recurse, f)
                },
            _ => f(&self.data),
        }
    }
}

impl<P: Clone, O, D> ObjectQuery for Tree<P, O, D> {
    fn query_objects<R, F>(&self, recurse: &R, f: &mut F)
        where R: Fn(&Tree<P, O, D>) -> bool,
              F: FnMut(&O),
    {
        match self.state {
            NodeState::Branch(ref nodes) if recurse(self) =>
                for node in nodes.iter() {
                    node.query_objects(recurse, f)
                },
            NodeState::Leaf(ref obj) => f(obj),
            _ => (),
        }
    }
}

impl<'a, P: Clone + 'a, O: 'a, D: 'a> IntoIterator for &'a Tree<P, O, D> {
    type Iter = Iter<'a, Tree<P, O, D>>;
    fn into_iter(self) -> Iter<'a, Tree<P, O, D>> { Iter::new(self) }
}


#[cfg(test)]
mod test {
    use std::rand::distributions::{IndependentSample, Range};
    use std::rand::thread_rng;
    use test::Bencher;
    use nalgebra::{Pnt2, Vec2, Orig};
    use quickcheck::quickcheck;

    use partition::Ncube;
    use traits::{NodeState, Node, Positioned};
    use super::*;

    #[test]
    fn tree_insert_into_empty() {
        let mut n = Tree::empty(Ncube::new(Pnt2::new(0.0f32, 0.0), 10.0f32), ());
        n.insert(Positioned { object: (), position: Pnt2::new(1.0f32, 0.0) }, ());
        match n.state {
            NodeState::Leaf(_) => (),
            _ => panic!("node is no leaf")
        }
    }

    #[test]
    fn tree_branch_on_second_insert() {
        let mut n = Tree::empty(Ncube::new(Pnt2::new(0.0, 0.0), 10.0), ());
        n.insert(Positioned { object: 1i32, position: Pnt2::new(1.0, -2.0) }, ());
        n.insert(Positioned { object: 2, position: Pnt2::new(2.0, 1.0) }, ());
        match n.state {
            NodeState::Branch(nodes) => {
                for k in 1..3 {
                    assert!(nodes.iter().any(|node| match node.state {
                        NodeState::Leaf(ref entry) => entry.object == k,
                        _ => false,
                    }));
                }
            },
            _ => panic!("node is no branch"),
        }
    }

    #[test]
    fn tree_from_empty_vec() {
        let tree: Tree<Ncube<Pnt2<f64>, f64>, Positioned<u8, Pnt2<f64>>, ()> =
            Tree::new(
                vec![].into_iter(),
                Ncube::new(Pnt2::new(0.0, 0.0), 1.0),
                (), &|_| (), &|_, _| ()
            );
        match tree.state {
            NodeState::Empty => (),
            _ => panic!(),
        }
    }

    #[test]
    fn tree_from_iter_more_than_two_branches() {
        fn tree_from_iter_more_than_two_branches(data: Vec<(f64, f64)>) -> bool {
            let tree = Tree::new(
                data.iter()
                .map(|&(x, y)| Positioned {
                    object: (),
                    position: Pnt2::new(x, y),
                }),
                Ncube::new(Orig::orig(), 200.0),
                (), &|_| (), &|_, _| ()
            );
            (data.len() >= 2) == (
                match tree.state {
                    NodeState::Branch(_) => true,
                    _ => false,
                }
            )
        }
        quickcheck(tree_from_iter_more_than_two_branches as fn(data: Vec<(f64, f64)>) -> bool)
    }

    #[test]
    fn tree_from_iter_one_is_a_leaf() {
        fn tree_from_iter_one_is_a_leaf(data: Vec<(f64, f64)>) -> bool {
            let tree = Tree::new(
                data.iter()
                .map(|&(x, y)| Positioned { object: (), position: Pnt2::new(x, y) }),
                Ncube::new(Orig::orig(), 200.0),
                (), &|_| (), &|_, _| ()
            );
            (data.len() == 1) == (
                match tree.state {
                    NodeState::Leaf(_) => true,
                    _ => false,
                }
            )
        }
        quickcheck(tree_from_iter_one_is_a_leaf as fn(data: Vec<(f64, f64)>) -> bool)
    }

    #[bench]
    fn tree_quad_with_center_of_mass_new_1000(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = thread_rng();
        let vec: Vec<_> = (0..1000).map(|_| Positioned {
            object: 1.0,
            position: Pnt2::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        }).collect();
        b.iter(|| {
            Tree::new(
                vec.iter().map(|a| a.clone()),
                Ncube::new(Orig::orig(), 2.0),
                (Vec2::new(0.0f64, 0.0), 0.0f64),
                &|obj| (obj.position.to_vec() * obj.object, obj.object),
                &|&(mps, ms), &(mp, m)| (mps + mp, ms + m)
            )
        })
    }

}
