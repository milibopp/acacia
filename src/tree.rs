use nalgebra::{Dim, BaseFloat, FloatPnt, FloatVec, zero, POrd};
use std::num::{Int, cast};
use std::cmp::partial_max;
use std::iter::AdditiveIterator;
use util::limits;


pub trait Positionable<P> {
    fn position(&self) -> P;
}

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


/// Subdivision helper function
///
/// Given the center of a node and its extent, this function generates a vector
/// of empty nodes equally partitioning the given node geometry. This is generic
/// over the dimension of the point type.
///
/// # Params
/// - The `old_center` is the center point of the (old) node.
/// - The `old_extent` is its extent.
fn subdivide<P, N>(old_center: &P, old_extent: &N) -> Vec<(P, N)>
    where P: Dim + Index<uint, N> + IndexMut<uint, N> + Copy,
          N: BaseFloat,
{
    let dim = Dim::dim(None::<P>);
    let new_extent = *old_extent / cast(2.0f64).unwrap();
    range(0u, 2.pow(dim))
        .map(|n| {
            let mut new_center = *old_center;
            for i in range(0, dim) {
                new_center[i] = new_center[i] + match n / 2.pow(i) % 2 {
                    0 => -new_extent,
                    1 => new_extent,
                    _ => unreachable!(),
                };
            }
            (new_center, new_extent)
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


enum NodeState<O, P, N> {
    Empty,
    Leaf(O),
    Branch(Vec<Node<O, P, N>>),
}


pub struct Node<O, P, N> {
    state: Option<NodeState<O, P, N>>,
    center: P,
    extent: N,
}

impl<O, P, N> Node<O, P, N> {
    pub fn empty(center: P, extent: N) -> Node<O, P, N> {
        Node {
            state: Some(NodeState::Empty),
            center: center,
            extent: extent,
        }
    }
}

impl<O, P, N, V> Node<O, P, N>
    where O: Positionable<P>,
          P: FloatPnt<N, V> + Index<uint, N> + IndexMut<uint, N> + Copy + POrd,
          V: FloatVec<N>,
          N: BaseFloat,
{
    pub fn from_vec(vec: Vec<O>) -> Node<O, P, N> {
        let (inf, sup) = limits(vec.iter().map(|obj| obj.position()));
        let center = (inf + sup.to_vec()) / cast(2.0f64).unwrap();
        let extent = range(0, Dim::dim(None::<P>))
            .fold(zero(), |max, n| partial_max(max, sup[n] - inf[n]).unwrap());
        let mut tree = Node::empty(center, extent);
        for object in vec.into_iter() {
            tree.insert(object)
        }
        tree
    }
}

impl<O, P, N, V> Node<O, P, N>
    where O: Positionable<P>,
          P: FloatPnt<N, V> + Index<uint, N> + IndexMut<uint, N> + Copy,
          N: BaseFloat,
{
    pub fn insert(&mut self, object: O) {
        let tmp = self.state.take().unwrap();
        self.state = Some(match tmp {
            NodeState::Empty => NodeState::Leaf(object),
            NodeState::Leaf(other) => {
                let mut nodes: Vec<Node<O, P, N>> = subdivide(&self.center, &self.extent)
                    .into_iter()
                    .map(|(p, n)| Node::empty(p, n))
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


#[cfg(test)]
mod test {
    use super::{Node, Entry, NodeState, subdivide, branch_dispatch};
    use nalgebra::{ApproxEq, Pnt2, Pnt3, FloatPnt};
    use std::num::Float;
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
    fn subdivide_new_extents_half_width((x, y, ext): (f64, f64, f64)) -> bool {
        let new_coords = subdivide(&Pnt2::new(x, y), &ext);
        new_coords.iter().all(|&(_, new_ext)|
            new_ext == ext / 2.0
        )
    }

    #[quickcheck]
    fn subdivide_new_centers_dist((x, y, ext): (f64, f64, f64)) -> TestResult {
        if ext > 0.0 {
            let center = Pnt2::new(x, y);
            let new_coords = subdivide(&center, &ext);
            TestResult::from_bool(new_coords.iter().all(|&(new_center, _)| {
                ApproxEq::approx_eq(
                    &FloatPnt::dist(&new_center, &center),
                    &(ext * 2.0.sqrt() / 2.0)
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
        let ext = 10.0f64;
        let subs = subdivide(&center, &ext);
        assert_eq!(
            subs[branch_dispatch(&center, &Pnt2::new(2.0, 2.0))].0,
            Pnt2::new(5.0, 5.0)
        );
        assert_eq!(
            subs[branch_dispatch(&center, &Pnt2::new(-2.0, 2.0))].0,
            Pnt2::new(-5.0, 5.0)
        );
    }

    #[test]
    fn branch_dispatch_cases_3d() {
        let center = Pnt3::new(0.0f64, 0.0, 0.0);
        let ext = 8.0f64;
        let subs = subdivide(&center, &ext);
        assert_eq!(
            subs[branch_dispatch(&center, &Pnt3::new(3.0, 1.0, -1.2))].0,
            Pnt3::new(4.0, 4.0, -4.0)
        );
        assert_eq!(
            subs[branch_dispatch(&center, &Pnt3::new(-2.0, 2.0, -0.1))].0,
            Pnt3::new(-4.0, 4.0, -4.0)
        );
    }

    #[test]
    fn node_insert_into_empty() {
        let mut n = Node::empty(Pnt2::new(0.0f32, 0.0), 10.0f32);
        n.insert(Entry { position: Pnt2::new(1.0f32, 0.0), object: 1i });
        match n.state {
            Some(NodeState::Leaf(entry)) => assert_eq!(entry.object, 1),
            _ => panic!("node is no leaf")
        }
    }

    #[test]
    fn node_insert_into_leaf() {
        let mut n = Node::empty(Pnt2::new(0.0f64, 0.0), 10.0f64);
        n.insert(Entry { object: 1i, position: Pnt2::new(1.0f64, -2.0) });
        n.insert(Entry { object: 2i, position: Pnt2::new(2.0, 1.0) });
        match n.state {
            Some(NodeState::Branch(nodes)) => {
                for &k in [1, 2].iter() {
                    assert!(nodes.iter().any(|node| match node.state {
                        Some(NodeState::Leaf(ref entry)) => entry.object == k,
                        _ => false,
                    }));
                }
            },
            _ => panic!("node is no branch"),
        }
    }

    #[test]
    fn node_branch_on_second_insert() {
        let mut n = Node::empty(Pnt2::new(0.0f64, 0.0), 8.0f64);
        n.insert(Entry { object: 1u, position: Pnt2::new(1.0, 2.0) });
        n.insert(Entry { object: 1u, position: Pnt2::new(2.0, -3.0) });
        match n.state.unwrap() {
            NodeState::Branch(_) => (),
            _ => panic!("node is no branch"),
        }
    }

    #[test]
    fn node_from_empty_vec() {
        let tree: Node<Entry<uint, Pnt2<f64>>, Pnt2<f64>, f64> =
            Node::from_vec(vec![]);
        match tree.state {
            Some(NodeState::Empty) => (),
            _ => panic!(),
        }
    }

    #[quickcheck]
    fn node_from_vec_more_than_two_branches(data: Vec<(uint, f64, f64)>) -> bool {
        let tree = Node::from_vec(
            data.iter()
            .map(|&(i, x, y)| Entry { object: i, position: Pnt2::new(x, y) })
            .collect()
        );
        (data.len() >= 2) == (
            match tree.state {
                Some(NodeState::Branch(_)) => true,
                _ => false,
            }
        )
    }

    #[quickcheck]
    fn node_from_vec_one_is_a_leaf(data: Vec<(uint, f64, f64)>) -> bool {
        let tree = Node::from_vec(
            data.iter()
            .map(|&(i, x, y)| Entry { object: i, position: Pnt2::new(x, y) })
            .collect()
        );
        (data.len() == 1) == (
            match tree.state {
                Some(NodeState::Leaf(_)) => true,
                _ => false,
            }
        )
    }
}
