use nalgebra::{Dim, BaseFloat, FloatPnt, FloatVec, zero, POrd};
use std::num::{Float, Int, cast};
use std::cmp::partial_max;
use std::iter::AdditiveIterator;
use std::mem;
use util::limits;


pub trait Positionable<P> {
    fn position(&self) -> P;
}

#[deriving(Clone)]
pub struct Entry<O, P> {
    pub object: O,
    pub position: P,
}

impl<O, P> Positionable<P> for Entry<O, P>
    where P: Copy
{
    fn position(&self) -> P {
        self.position
    }
}


/// A generic spatial tree
///
/// This trait wraps up the properties of a tree. Most prominently a tree must
/// support recursive computational queries on its structure. Closures are used
/// to determine the recursion behavior and what is to be computed.
///
/// # Type parameters
///
/// - `P` is the kind of point used to position objects and nodes spatially.
/// - `N` is the scalar of the vector space of points.
/// - The tree stores objects of type `O`. These objects need to have some
///   notion of a position.
/// - `D` is the kind of data associated with each node. This is computed
///   recursively during tree construction.
pub trait Tree<P, N, O, D> {

    /// Compute a query on the associated data using a mutable accumulator
    ///
    /// This method walks recursively through the tree, as deep as `recurse`
    /// prescribes, and `combine` data subsequently modifying an `accumulator`.
    ///
    /// If an empty or leaf node is encountered, `combine` is called on the
    /// accumulator and its associated data. For a branching node `recurse` is
    /// called on its center, width and associated data, to determine whether
    /// its subnodes should be inspected more closely. If so, the function
    /// recurses on each subnode, otherwise it acts on it as if it were not a
    /// branch.
    ///
    /// # Parameters
    ///
    /// - The `accumulator` is a mutable reference to some data, that is
    ///   modified during the query to collect.
    /// - At each node the tree is only recursed further, if
    ///   `recurse(&node.center, &node.width, &node.data)`.
    /// - `combine` is called for every node whose data is to be considered.
    fn query_data_mut<T>(&self, accumulator: &mut T, recurse: |&P, &N, &D| -> bool, combine: |&mut T, &D|);

    /// Compute a query on the associated data using an initial state
    ///
    /// This is very similar to the accumulator variant and indeed has a default
    /// implementation using it. The difference is, that an initial value is
    /// moved into the function and used to initialize the accumulator. Its
    /// final state is the method's return value.
    fn query_data<T>(&self, initial: T, recurse: |&P, &N, &D| -> bool, combine: |&T, &D| -> T) -> T {
        let mut acc = initial;
        self.query_data_mut(
            &mut acc,
            |ctr, width, data| recurse(ctr, width, data),
            |a, d| {*a = combine(a, d);}
        );
        acc
    }

    /// Compute a query on the objects using an accumulator
    ///
    /// This method walks through the tree similar to `query_data_mut`. However,
    /// the `combine` closure is only invoked, when a leaf node is encountered.
    /// It recurses on branch nodes if `recurse(...)`. Empty nodes are
    /// ignored.
    fn query_objects_mut<T>(&self, accumulator: &mut T, recurse: |&P, &N, &D| -> bool, combine: |&mut T, &O|);

    /// Compute a query on the objects using an initial state
    ///
    /// This relates to `query_objects_mut` in the same way `query_data` relates
    /// to `query_data_mut`.
    fn query_objects<T>(&self, initial: T, recurse: |&P, &N, &D| -> bool, combine: |&T, &O| -> T) -> T {
        let mut acc = initial;
        self.query_objects_mut(
            &mut acc,
            |ctr, width, data| recurse(ctr, width, data),
            |a, o| {*a = combine(a, o);}
        );
        acc
    }
}


/// A node in a tree
///
/// This is part of the essential features of a tree. Note that both a whole
/// tree and its constituents implement this.
pub trait Node<P, N, O> {
    fn state(&self) -> &NodeState<O, Self>;
    fn center(&self) -> &P;
    fn width(&self) -> &N;
}

/// A node which has data associated with it
pub trait AssociatedData<D> {
    fn data(&self) -> &D;
}


/// Subdivision helper function
///
/// Given the center of a node and its width, this function generates a vector
/// of empty nodes equally partitioning the given node geometry. This is generic
/// over the dimension of the point type.
///
/// # Params
/// - The `old_center` is the center point of the (old) node.
/// - The `old_width` is its width.
fn subdivide<P, N>(old_center: &P, old_width: &N) -> Vec<(P, N)>
    where P: Dim + Index<uint, N> + IndexMut<uint, N> + Copy,
          N: BaseFloat,
{
    let _2 = cast(2.0f64).unwrap();
    let dim = Dim::dim(None::<P>);
    let new_width = *old_width / _2;
    range(0u, 2.pow(dim))
        .map(|n| {
            let mut new_center = *old_center;
            let dx = new_width / _2;
            for i in range(0, dim) {
                new_center[i] = new_center[i] + match n / 2.pow(i) % 2 {
                    0 => -dx,
                    1 => dx,
                    _ => unreachable!(),
                };
            }
            (new_center, new_width)
        })
        .collect()
}


fn branch_dispatch<P, N>(center: &P, point: &P) -> uint
    where P: Dim + Index<uint, N> + IndexMut<uint, N>,
          N: BaseFloat,
{
    let dim = Dim::dim(None::<P>);
    range(0, dim)
        .map(|k| if point[k] < center[k] {0} else {1 << k})
        .sum()
}


pub enum NodeState<O, N> {
    Empty,
    Leaf(O),
    Branch(Vec<N>),
}


#[unstable="must be refactored to implement the new traits"]
pub struct PureNTree<O, P, N> {
    state: Option<NodeState<O, PureNTree<O, P, N>>>,
    center: P,
    width: N,
}

impl<O, P, N> PureNTree<O, P, N> {
    pub fn empty(center: P, width: N) -> PureNTree<O, P, N> {
        PureNTree {
            state: Some(NodeState::Empty),
            center: center,
            width: width,
        }
    }
}

impl<O, P, N, V> PureNTree<O, P, N>
    where O: Positionable<P>,
          P: FloatPnt<N, V> + Index<uint, N> + IndexMut<uint, N> + Copy + POrd,
          V: FloatVec<N>,
          N: BaseFloat,
{
    /// Construct a tree without checking the geometry of the input data
    ///
    /// Note: this is prone to stack overflows! By calling this you effectively
    /// assert that all positions are within the tree bounds.
    pub fn from_iter_raw<I: Iterator<O>>(iter: I, center: P, width: N) -> PureNTree<O, P, N> {
        let mut tree = PureNTree::empty(center, width);
        let mut iter = iter;
        for object in iter {
            tree.insert(object)
        }
        tree
    }

    pub fn from_iter<I: Iterator<O>>(iter: I) -> PureNTree<O, P, N> {
        let _2: N = cast(2.0f64).unwrap();
        let vec: Vec<O> = iter.collect();
        let (inf, sup) = limits(vec.iter().map(|obj| obj.position()));
        let center = (inf + sup.to_vec()) / cast(2.0f64).unwrap();
        let width = _2 * range(0, Dim::dim(None::<P>))
            .fold(zero(), |max, n| partial_max(max, sup[n] - inf[n]).unwrap());
        PureNTree::from_iter_raw(vec.into_iter(), center, width)
    }

    pub fn from_iter_with_geometry<I: Iterator<O>>(iter: I, center: P, minimal_width: N) -> PureNTree<O, P, N> {
        let _2: N = cast(2.0f64).unwrap();
        let vec: Vec<O> = iter.collect();
        let (inf, sup) = limits(vec.iter().map(|obj| obj.position()));
        let width = _2 * range(0, Dim::dim(None::<P>))
            .fold(minimal_width, |max, n|
                partial_max(
                    max, partial_max(
                        (center[n] - sup[n]).abs(),
                        (center[n] - inf[n]).abs()
                    ).unwrap()
                ).unwrap()
            );
        PureNTree::from_iter_raw(vec.into_iter(), center, width)
    }
}

impl<O, P, N, V> PureNTree<O, P, N>
    where O: Positionable<P>,
          P: FloatPnt<N, V> + Index<uint, N> + IndexMut<uint, N> + Copy,
          N: BaseFloat,
{
    pub fn insert(&mut self, object: O) {
        let tmp = self.state.take().unwrap();
        self.state = Some(match tmp {
            NodeState::Empty => NodeState::Leaf(object),
            NodeState::Leaf(other) => {
                let mut nodes: Vec<PureNTree<O, P, N>> = subdivide(&self.center, &self.width)
                    .into_iter()
                    .map(|(p, n)| PureNTree::empty(p, n))
                    .collect();
                nodes[branch_dispatch(&self.center, &other.position())].insert(other);
                nodes[branch_dispatch(&self.center, &object.position())].insert(object);
                NodeState::Branch(nodes)
            },
            NodeState::Branch(mut nodes) => {
                nodes[branch_dispatch(&self.center, &object.position())].insert(object);
                NodeState::Branch(nodes)
            },
        });
    }
}


/// An N-dimensional tree
///
/// This tree does not know the dimension of its point at compile time, as it is
/// not hard-coded and genericity over constants is unsupported in Rust.
///
pub struct NTree<P, N, O, D> {
    state: NodeState<O, NTree<P, N, O, D>>,
    center: P,
    width: N,
    data: D,
}

impl<P, N, O, D> NTree<P, N, O, D> {

    /// Construct an empty tree
    fn empty(center: P, width: N, data: D) -> NTree<P, N, O, D> {
        NTree {
            state: NodeState::Empty,
            center: center,
            width: width,
            data: data,
        }
    }
}

impl<P, N, O, D: Clone> NTree<P, N, O, D> {

    /// Recompute the associated data
    fn recompute_data(&mut self, default: D, single: |&O| -> D, combine: |&D, &D| -> D) {
        self.data = match self.state {
            NodeState::Empty => default,
            NodeState::Leaf(ref obj) => single(obj),
            NodeState::Branch(ref mut nodes) => {
                for node in nodes.iter_mut() {
                    node.recompute_data(default.clone(), |o| single(o), |d1, d2| combine(d1, d2));
                }
                nodes.iter().fold(default.clone(), |current, node| combine(&current, &node.data))
            },
        };
    }
}

impl<P, N, O, D> NTree<P, N, O, D>
    where O: Positionable<P>,
          P: Dim + Index<uint, N> + IndexMut<uint, N> + Copy,
          N: BaseFloat,
          D: Clone,
{

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
                let mut nodes: Vec<NTree<P, N, O, D>> = subdivide(&self.center, &self.width)
                    .into_iter()
                    .map(|(p, n)| NTree::empty(p, n, default.clone()))
                    .collect();
                nodes[branch_dispatch(&self.center, &other.position())].insert(other, default.clone());
                nodes[branch_dispatch(&self.center, &object.position())].insert(object, default.clone());
                NodeState::Branch(nodes)
            },
            NodeState::Branch(mut nodes) => {
                nodes[branch_dispatch(&self.center, &object.position())].insert(object, default.clone());
                NodeState::Branch(nodes)
            },
        };
    }

    /// Construct the tree from an iterator with fixed center and width
    ///
    /// NOTE: this is prone to stack overflows! By calling this you effectively
    /// assert that all positions are within the tree bounds.
    fn from_iter_raw<I: Iterator<O>>(objects: I, center: P, width: N, default: D, single: |&O| -> D, combine: |&D, &D| -> D) -> NTree<P, N, O, D> {
        let mut objects = objects;
        let mut tree = NTree::empty(center, width, default.clone());
        for object in objects {
            tree.insert(object, default.clone());
        }
        tree.recompute_data(default.clone(), |o| single(o), |d1, d2| combine(d1, d2));
        tree
    }
}

impl<P, N, O, D, V> NTree<P, N, O, D>
    where O: Positionable<P>,
          P: FloatPnt<N, V> + Index<uint, N> + IndexMut<uint, N> + Copy + POrd,
          V: FloatVec<N>,
          N: BaseFloat,
          D: Clone,
{

    /// Construct a tree from an iterator
    ///
    /// The center and width of the root node should be determined automatically
    /// from the structure of the objects provided by the iterator.
    ///
    /// # Parameters
    ///
    /// - `objects` is an iterator over the objects that should comprise the
    ///   tree.
    /// - The `default` value is the associated data of an empty node.
    /// - A terminal leaf will be associated with data using the `single`
    ///   closure.
    /// - More complex structures will combine their constituent data by
    ///   subsequent invokations of the `combine` function as an operator.
    pub fn from_iter<I: Iterator<O>>(objects: I, default: D, single: |&O| -> D, combine: |&D, &D| -> D) -> NTree<P, N, O, D> {
        let _2: N = cast(2.0f64).unwrap();
        let vec: Vec<O> = objects.collect();
        let (inf, sup) = limits(vec.iter().map(|obj| obj.position()));
        let center = (inf + sup.to_vec()) / _2;
        let width = _2 * range(0, Dim::dim(None::<P>))
            .fold(zero(), |max, n| partial_max(max, sup[n] - inf[n]).unwrap());
        NTree::from_iter_raw(vec.into_iter(), center, width, default, |o| single(o), |d1, d2| combine(d1, d2))
    }

    /// Same as `from_iter` but with geometrical constraints
    ///
    /// # Parameters
    ///
    /// - `objects`, `default`, `single`, `combine` are the same as in
    ///   `from_iter`.
    /// - `center` will be the center of the tree.
    /// - `minimal_width` is the minimal width of the root node of the tree.
    ///   This may be increased by the geometrical requirements of the objects,
    ///   thus only a lower bound on the width can be imposed on the tree.
    pub fn from_iter_with_geometry<I: Iterator<O>>(objects: I, center: P, minimal_width: N, default: D, single: |&O| -> D, combine: |&D, &D| -> D) -> NTree<P, N, O, D> {
        let _2: N = cast(2.0f64).unwrap();
        let vec: Vec<O> = objects.collect();
        let (inf, sup) = limits(vec.iter().map(|obj| obj.position()));
        let width = _2 * range(0, Dim::dim(None::<P>))
            .fold(minimal_width, |max, n|
                partial_max(
                    max, partial_max(
                        (center[n] - sup[n]).abs(),
                        (center[n] - inf[n]).abs()
                    ).unwrap()
                ).unwrap()
            );
        NTree::from_iter_raw(vec.into_iter(), center, width, default, |o| single(o), |d1, d2| combine(d1, d2))
    }
}

impl<P, N, O, D> Node<P, N, O> for NTree<P, N, O, D> {
    fn state(&self) -> &NodeState<O, NTree<P, N, O, D>> {
        &self.state
    }

    fn center(&self) -> &P {
        &self.center
    }

    fn width(&self) -> &N {
        &self.width
    }
}

impl<P, N, O, D> AssociatedData<D> for NTree<P, N, O, D> {
    fn data(&self) -> &D {
        &self.data
    }
}

impl<P, N, O, D> Tree<P, N, O, D> for NTree<P, N, O, D> {
    fn query_data_mut<T>(&self, acc: &mut T, recurse: |&P, &N, &D| -> bool, combine: |&mut T, &D|) {
        match self.state {
            NodeState::Branch(ref nodes) if recurse(&self.center, &self.width, &self.data) =>
                for node in nodes.iter() {
                    node.query_data_mut(acc, |p, n, d| recurse(p, n, d), |t, d| combine(t, d))
                },
            _ => combine(acc, &self.data),
        }
    }

    fn query_objects_mut<T>(&self, acc: &mut T, recurse: |&P, &N, &D| -> bool, combine: |&mut T, &O|) {
        match self.state {
            NodeState::Branch(ref nodes) if recurse(&self.center, &self.width, &self.data) =>
                for node in nodes.iter() {
                    node.query_objects_mut(acc, |p, n, d| recurse(p, n, d), |t, o| combine(t, o))
                },
            NodeState::Leaf(ref obj) => combine(acc, obj),
            _ => (),
        }
    }
}


#[cfg(test)]
mod test {
    use super::{Entry, NodeState, NTree, Tree, subdivide, branch_dispatch};
    use std::num::Float;
    use std::rand::distributions::{IndependentSample, Range};
    use std::rand::task_rng;
    use test::Bencher;
    use nalgebra::{ApproxEq, Pnt2, Pnt3, FloatPnt, Vec2, Vec3, zero, Norm, Orig};
    use quickcheck::TestResult;

    #[test]
    fn subdivide_2d_4nodes() {
        let nodes = subdivide(&Pnt2::new(0.0f64, 0.0), &1.0f64);
        assert_eq!(nodes.len(), 4);
    }

    #[test]
    fn subdivide_3d_8nodes() {
        let nodes = subdivide(&Pnt3::new(0.0f64, 0.0, 0.0), &1.0f64);
        assert_eq!(nodes.len(), 8);
    }

    #[quickcheck]
    fn subdivide_new_width_half((x, y, width): (f64, f64, f64)) -> bool {
        let new_coords = subdivide(&Pnt2::new(x, y), &width);
        new_coords.iter().all(|&(_, new_width)|
            new_width == width / 2.0
        )
    }

    #[quickcheck]
    fn subdivide_new_centers_dist((x, y, width): (f64, f64, f64)) -> TestResult {
        if width > 0.0 {
            let center = Pnt2::new(x, y);
            let new_coords = subdivide(&center, &width);
            TestResult::from_bool(new_coords.iter().all(|&(new_center, _)| {
                ApproxEq::approx_eq(
                    &FloatPnt::dist(&new_center, &center),
                    &(width * 2.0.powf(-1.5))
                )
            }))
        } else {
            TestResult::discard()
        }
    }

    #[quickcheck]
    fn branch_dispatch_index_range_2d((px, py): (f64, f64)) -> bool {
        // TODO: make the center variable
        branch_dispatch(&Pnt2::new(0.0, 0.0), &Pnt2::new(px, py)) < 4
    }

    #[quickcheck]
    fn branch_dispatch_index_range_3d((px, py, pz): (f64, f64, f64)) -> bool {
        // TODO: make the center variable
        branch_dispatch(&Pnt3::new(0.0, 0.0, 0.0), &Pnt3::new(px, py, pz)) < 8
    }

    #[test]
    fn branch_dispatch_cases_2d() {
        let center = Pnt2::new(0.0f64, 0.0);
        let width = 10.0f64;
        let subs = subdivide(&center, &width);
        assert_eq!(
            subs[branch_dispatch(&center, &Pnt2::new(2.0, 4.0))].0,
            Pnt2::new(2.5, 2.5)
        );
        assert_eq!(
            subs[branch_dispatch(&center, &Pnt2::new(-1.0, 2.0))].0,
            Pnt2::new(-2.5, 2.5)
        );
    }

    #[test]
    fn branch_dispatch_cases_3d() {
        let center = Pnt3::new(0.0f64, 0.0, 0.0);
        let width = 8.0f64;
        let subs = subdivide(&center, &width);
        assert_eq!(
            subs[branch_dispatch(&center, &Pnt3::new(3.0, 1.0, -1.2))].0,
            Pnt3::new(2.0, 2.0, -2.0)
        );
        assert_eq!(
            subs[branch_dispatch(&center, &Pnt3::new(-2.0, 2.0, -0.1))].0,
            Pnt3::new(-2.0, 2.0, -2.0)
        );
    }

    #[test]
    fn ntree_insert_into_empty() {
        let mut n = NTree::empty(Pnt2::new(0.0f32, 0.0), 10.0f32, ());
        n.insert(Entry { position: Pnt2::new(1.0f32, 0.0), object: 1i }, ());
        match n.state {
            NodeState::Leaf(entry) => assert_eq!(entry.object, 1),
            _ => panic!("node is no leaf")
        }
    }

    #[test]
    fn ntree_insert_into_leaf() {
        let mut n = NTree::empty(Pnt2::new(0.0f64, 0.0), 10.0f64, ());
        n.insert(Entry { object: 1i, position: Pnt2::new(1.0f64, -2.0) }, ());
        n.insert(Entry { object: 2i, position: Pnt2::new(2.0, 1.0) }, ());
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
        n.insert(Entry { object: 1u, position: Pnt2::new(1.0, 2.0) }, ());
        n.insert(Entry { object: 1u, position: Pnt2::new(2.0, -3.0) }, ());
        match n.state {
            NodeState::Branch(_) => (),
            _ => panic!("node is no branch"),
        }
    }

    #[test]
    fn ntree_from_empty_vec() {
        let tree: NTree<Pnt2<f64>, f64, Entry<uint, Pnt2<f64>>, ()> =
            NTree::from_iter(vec![].into_iter(), (), |_| (), |_, _| ());
        match tree.state {
            NodeState::Empty => (),
            _ => panic!(),
        }
    }

    #[quickcheck]
    fn ntree_from_iter_more_than_two_branches(data: Vec<(uint, f64, f64)>) -> bool {
        let tree = NTree::from_iter(
            data.iter()
            .map(|&(i, x, y)| Entry { object: i, position: Pnt2::new(x, y) }),
            (), |_| (), |_, _| ()
        );
        (data.len() >= 2) == (
            match tree.state {
                NodeState::Branch(_) => true,
                _ => false,
            }
        )
    }

    #[quickcheck]
    fn ntree_from_iter_one_is_a_leaf(data: Vec<(uint, f64, f64)>) -> bool {
        let tree = NTree::from_iter(
            data.iter()
            .map(|&(i, x, y)| Entry { object: i, position: Pnt2::new(x, y) }),
            (), |_| (), |_, _| ()
        );
        (data.len() == 1) == (
            match tree.state {
                NodeState::Leaf(_) => true,
                _ => false,
            }
        )
    }

    #[bench]
    fn quadtree_raw_tree_construction_uniform(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = task_rng();
        let vec = Vec::from_fn(1000, |_| Entry {
            object: 1i,
            position: Pnt2::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        });
        b.iter(|| {
            NTree::from_iter_raw(
                vec.iter().map(|&a| a.clone()),
                Orig::orig(), 2.0,
                (), |_| (), |_, _| ()
            )
        })
    }

    #[bench]
    fn quadtree_construction_uniform(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = task_rng();
        let vec = Vec::from_fn(1000, |_| Entry {
            object: 1i,
            position: Pnt2::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        });
        b.iter(|| {
            let ntree = NTree::from_iter(
                vec.iter().map(|&a| a.clone()),
                (), |_| (), |_, _| (),
            );
            ntree
        })
    }

    #[bench]
    fn octree_construction_uniform(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = task_rng();
        let vec = Vec::from_fn(1000, |_| Entry {
            object: 1i,
            position: Pnt3::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        });
        b.iter(|| {
            let ntree = NTree::from_iter(
                vec.iter().map(|&a| a.clone()),
                (), |_| (), |_, _| (),
            );
            ntree
        })
    }

    #[bench]
    fn quadtree_associate_data(b: &mut Bencher) {
        let coord_dist = Range::new(-1.0f64, 1.0);
        let mut rng = task_rng();
        let vec = Vec::from_fn(1000, |_| Entry {
            object: 1.0,
            position: Pnt2::new(
                coord_dist.ind_sample(&mut rng),
                coord_dist.ind_sample(&mut rng)
            ),
        });
        b.iter(|| {
            let ntree = NTree::from_iter(
                vec.iter().map(|&a| a.clone()),
                (Vec2::new(0.0f64, 0.0), 0.0f64),
                |obj| (obj.position.to_vec() * obj.object, obj.object),
                |&(mps, ms), &(mp, m)| (mps + mp, ms + m)
            );
            ntree
        })
    }

    #[quickcheck]
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
                Entry { object: m, position: Pnt2::new(x, y) }
            ),
            (Vec2::new(0.0f64, 0.0), 0.0f64),
            |obj| (obj.position.to_vec() * obj.object, obj.object),
            |&(mps, ms), &(mp, m)| (mps + mp, ms + m)
        );
        let (tree_mps, tree_ms) = tree.data;
        // …and compare
        TestResult::from_bool(ApproxEq::approx_eq(&(tree_mps / tree_ms), &com))
    }

    #[quickcheck]
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
        let pnt = |p| {
            let (x, y, z) = p;
            Pnt3::new(x, y, z)
        };
        let test_point = pnt(test_point_);
        // Newton's law of gravity for two point masses (with G = 1)
        let newton = |(m, p1): (f64, Pnt3<f64>), p2| {
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
                Entry { object: m, position: Pnt3::new(x, y, z) }
            ),
            orig,
            test_point.as_vec().norm() * 2.0,
            (orig, zero()),
            |obj| (obj.position, obj.object),
            |&(com1, m1), &(com2, m2)|
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
        tree.query_data_mut(
            &mut tree_gravity,
            |&node_center, &node_size, &(center_of_mass, _)| {
                let d = FloatPnt::dist(&test_point, &center_of_mass);
                let delta = FloatPnt::dist(&node_center, &center_of_mass);
                d < node_size / theta + delta
            },
            |g, &(com, m)| {
                *g = *g + newton((m, com), test_point);
            },
        );
        // Now the tree gravity should approximate the exact one, within 5 %
        TestResult::from_bool(simple_gravity.approx_eq_eps(&tree_gravity, &(0.05 * simple_gravity.norm())))
    }
}
