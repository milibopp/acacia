//! Dimension-unspecific tree implementations

use std::mem;
use tree::{NodeState, DataQuery, ObjectQuery, AssociatedData, Node, Position};
use partition::Partition;


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
        let mut iter = iter;
        for object in iter {
            tree.insert(object)
        }
        tree
    }

    fn dispatch(nodes: &mut Vec<PureTree<P, O>>, object: O) {
        for node in nodes.iter_mut() {
            if node.partition.contains(&object.position()) {
                node.insert(object);
                return;
            }
        }
        panic!("object could not be inserted into the tree");
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
                PureTree::dispatch(&mut nodes, object);
                PureTree::dispatch(&mut nodes, other);
                NodeState::Branch(nodes)
            },
            NodeState::Branch(mut nodes) => {
                PureTree::dispatch(&mut nodes, object);
                NodeState::Branch(nodes)
            }
        };
    }
}

impl<P, O> ObjectQuery for PureTree<P, O> {
    fn query_objects<R, F>(&self, recurse: &R, f: &mut F)
        where R: Fn(&PureTree<P, O>) -> bool,
              F: FnMut(&O),
    {
        match self.state {
            NodeState::Branch(ref nodes) => 
                if recurse(self) {
                    for node in nodes.iter() {
                        node.query_objects(recurse, f)
                    }
                },
            NodeState::Leaf(ref obj) => f(obj),
            _ => (),
        }
    }
}

impl<P, O> Node for PureTree<P, O> {
    type Partition = P;
    type Object = O;
    type Container = Vec<PureTree<P, O>>;

    fn state(&self) -> &NodeState<O, Vec<PureTree<P, O>>> {
        &self.state
    }

    fn partition(&self) -> &P {
        &self.partition
    }
}


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
    fn dispatch(nodes: &mut Vec<Tree<P, O, D>>, object: O, default: D) {
        for node in nodes.iter_mut() {
            if node.partition.contains(&object.position()) {
                node.insert(object, default);
                return;
            }
        }
        panic!("object could not be inserted into the tree");
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
                Tree::dispatch(&mut nodes, object, default.clone());
                Tree::dispatch(&mut nodes, other, default.clone());
                NodeState::Branch(nodes)
            },
            NodeState::Branch(mut nodes) => {
                Tree::dispatch(&mut nodes, object, default.clone());
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
        let mut objects = objects;
        let mut tree = Tree::empty(partition, default.clone());
        for object in objects {
            tree.insert(object, default.clone());
        }
        tree.recompute_data(default.clone(), single, combine);
        tree
    }
}


impl<P, O, D> Node for Tree<P, O, D> {
    type Partition = P;
    type Object = O;
    type Container = Vec<Tree<P, O, D>>;

    fn state(&self) -> &NodeState<O, Vec<Tree<P, O, D>>> {
        &self.state
    }

    fn partition(&self) -> &P {
        &self.partition
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

impl<P, O, D> ObjectQuery for Tree<P, O, D> {
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


#[cfg(test)]
mod test {
    use super::{Tree, PureTree};
    use tree::{NodeState, Node, ObjectQuery, Positioned};
    use std::num::Float;
    use std::rand::distributions::{IndependentSample, Range};
    use std::rand::thread_rng;
    use std::iter::AdditiveIterator;
    use test::Bencher;
    use nalgebra::{ApproxEq, Pnt2, Pnt3, FloatPnt, Vec2, Orig};
    use quickcheck::{TestResult, quickcheck};

    #[test]
    fn tree_insert_into_empty() {
        let mut n = Tree::empty(Pnt2::new(0.0f32, 0.0), 10.0f32, ());
        n.insert(Positioned { position: Pnt2::new(1.0f32, 0.0), object: 1i }, ());
        match n.state {
            NodeState::Leaf(entry) => assert_eq!(entry.object, 1),
            _ => panic!("node is no leaf")
        }
    }

    #[test]
    fn ntree_insert_into_leaf() {
        let mut n = NTree::empty(Pnt2::new(0.0f64, 0.0), 10.0f64, ());
        n.insert(Positioned { object: 1i, position: Pnt2::new(1.0f64, -2.0) }, ());
        n.insert(Positioned { object: 2i, position: Pnt2::new(2.0, 1.0) }, ());
        match n.state {
            NodeState::Branch(nodes) => {
                for &k in [1, 2].iter() {
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
    fn ntree_branch_on_second_insert() {
        let mut n = NTree::empty(Pnt2::new(0.0f64, 0.0), 8.0f64, ());
        n.insert(Positioned { object: 1u, position: Pnt2::new(1.0, 2.0) }, ());
        n.insert(Positioned { object: 1u, position: Pnt2::new(2.0, -3.0) }, ());
        match n.state {
            NodeState::Branch(_) => (),
            _ => panic!("node is no branch"),
        }
    }

    #[test]
    fn ntree_from_empty_vec() {
        let tree: NTree<Pnt2<f64>, f64, Positioned<uint, Pnt2<f64>>, ()> =
            NTree::from_iter(vec![].into_iter(), (), &|_| (), &|_, _| ());
        match tree.state {
            NodeState::Empty => (),
            _ => panic!(),
        }
    }

    #[test]
    fn ntree_from_iter_more_than_two_branches() {
        fn ntree_from_iter_more_than_two_branches(data: Vec<(uint, f64, f64)>) -> bool {
            let tree = NTree::from_iter(
                data.iter()
                .map(|&(i, x, y)| Positioned { object: i, position: Pnt2::new(x, y) }),
                (), &|_| (), &|_, _| ()
            );
            (data.len() >= 2) == (
                match tree.state {
                    NodeState::Branch(_) => true,
                    _ => false,
                }
            )
        }
        quickcheck(ntree_from_iter_more_than_two_branches as fn(data: Vec<(uint, f64, f64)>) -> bool)
    }

    #[test]
    fn ntree_from_iter_one_is_a_leaf() {
        fn ntree_from_iter_one_is_a_leaf(data: Vec<(uint, f64, f64)>) -> bool {
            let tree = NTree::from_iter(
                data.iter()
                .map(|&(i, x, y)| Positioned { object: i, position: Pnt2::new(x, y) }),
                (), &|_| (), &|_, _| ()
            );
            (data.len() == 1) == (
                match tree.state {
                    NodeState::Leaf(_) => true,
                    _ => false,
                }
            )
        }
        quickcheck(ntree_from_iter_one_is_a_leaf as fn(data: Vec<(uint, f64, f64)>) -> bool)
    }

    #[bench]
    fn pure_ntree_quad_from_iter_raw_1000(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = thread_rng();
        let vec: Vec<_> = range(0u, 1000).map(|_| Positioned {
            object: 1i,
            position: Pnt2::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        }).collect();
        b.iter(|| {
            PureNTree::from_iter_raw(
                vec.iter().map(|a| a.clone()),
                Orig::orig(), 2.0,
            )
        })
    }

    #[bench]
    fn pure_ntree_quad_from_iter_1000(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = thread_rng();
        let vec: Vec<_> = range(0u, 1000).map(|_| Positioned {
            object: 1i,
            position: Pnt2::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        }).collect();
        b.iter(||
            PureNTree::from_iter(vec.iter().map(|a| a.clone()))
        )
    }

    #[bench]
    fn pure_ntree_oc_from_iter_1000(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = thread_rng();
        let vec: Vec<_> = range(0u, 1000).map(|_| Positioned {
            object: 1i,
            position: Pnt3::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        }).collect();
        b.iter(|| {
            PureNTree::from_iter(vec.iter().map(|a| a.clone()))
        })
    }

    #[bench]
    fn ntree_quad_with_center_of_mass_from_iter_1000(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = thread_rng();
        let vec: Vec<_> = range(0u, 1000).map(|_| Positioned {
            object: 1.0,
            position: Pnt2::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        }).collect();
        b.iter(|| {
            NTree::from_iter(
                vec.iter().map(|a| a.clone()),
                (Vec2::new(0.0f64, 0.0), 0.0f64),
                &|obj| (obj.position.to_vec() * obj.object, obj.object),
                &|&(mps, ms), &(mp, m)| (mps + mp, ms + m)
            )
        })
    }

    #[bench]
    fn ntree_quad_with_center_of_mass_from_iter_raw_1000(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = thread_rng();
        let vec: Vec<_> = range(0u, 1000).map(|_| Positioned {
            object: 1.0,
            position: Pnt2::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        }).collect();
        b.iter(|| {
            NTree::from_iter_raw(
                vec.iter().map(|a| a.clone()),
                Orig::orig(), 2.0,
                (Vec2::new(0.0f64, 0.0), 0.0f64),
                &|obj| (obj.position.to_vec() * obj.object, obj.object),
                &|&(mps, ms), &(mp, m)| (mps + mp, ms + m)
            )
        })
    }

    #[bench]
    fn pure_ntree_query_objects(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = thread_rng();
        let objects: Vec<_> = range(0u, 1000).map(|_| Positioned {
            object: (),
            position: Pnt2::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        }).collect();
        let search_radius = 0.3;
        let tree = PureNTree::from_iter(objects.clone().into_iter());
        b.iter(|| {
            // Count the number of objects within the search radius 10000 times
            range(0u, 10000)
                .map(|_| {
                    let mut i = 0u;
                    tree.query_objects(
                        &|node| node.center().dist(&Orig::orig()) < search_radius + *node.width() / 2.0,
                        &mut |other| if other.position.dist(&Orig::orig()) < search_radius {i += 1},
                    );
                    i
                })
                .sum()
        })
    }
}
