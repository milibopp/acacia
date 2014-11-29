use nalgebra::{Dim, BaseFloat};
use std::num::{Int, cast};
use std::iter::AdditiveIterator;

pub struct Entry<O, P> {
    pub object: O,
    pub position: P,
}

pub enum NodeState<O, P, N> {
    Empty,
    Leaf(Entry<O, P>),
    Branch(Vec<Node<O, P, N>>),
}

pub struct Node<O, P, N> {
    pub state: Option<NodeState<O, P, N>>,
    pub center: P,
    pub extent: N,
}

impl<O, P, N> Node<O, P, N> {
    pub fn new(state: NodeState<O, P, N>, center: P, extent: N) -> Node<O, P, N> {
        Node {
            state: Some(state),
            center: center,
            extent: extent,
        }
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
pub fn subdivide<P, N>(old_center: &P, old_extent: &N) -> Vec<(P, N)>
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

pub fn branch_dispatch<P, N>(center: &P, point: &P) -> uint
    where P: Dim + Index<uint, N> + IndexMut<uint, N>,
          N: BaseFloat,
{
    let dim = Dim::dim(None::<P>);
    range(0, dim)
        .map(|k| if point[k] < center[k] {0} else {1 << k})
        .sum()
}

impl<O, P, N> Node<O, P, N>
    where P: Dim + Index<uint, N> + IndexMut<uint, N> + Copy,
          N: BaseFloat,
{
    pub fn insert(&mut self, entry: Entry<O, P>) {
        let tmp = self.state.take().unwrap();
        self.state = Some(match tmp {
            NodeState::Empty => NodeState::Leaf(entry),
            NodeState::Leaf(other) => {
                let mut nodes: Vec<Node<O, P, N>> = subdivide(&self.center, &self.extent)
                    .into_iter()
                    .map(|(p, n)| Node::new(NodeState::Empty, p, n))
                    .collect();
                nodes[branch_dispatch(&self.center, &other.position)].insert(other);
                nodes[branch_dispatch(&self.center, &entry.position)].insert(entry);
                NodeState::Branch(nodes)
            },
            NodeState::Branch(mut nodes) => {
                nodes[branch_dispatch(&self.center, &entry.position)].insert(entry);
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
        let mut n = Node::new(
            NodeState::Empty,
            Pnt2::new(0.0f32, 0.0),
            10.0f32,
        );
        n.insert(Entry { position: Pnt2::new(1.0f32, 0.0), object: 1i });
        match n.state {
            Some(NodeState::Leaf(entry)) => assert_eq!(entry.object, 1),
            _ => panic!("node is no leaf")
        }
    }

    #[test]
    fn node_insert_into_leaf() {
        let mut n = Node::new(
            NodeState::Leaf(Entry {object: 1i, position: Pnt2::new(1.0f64, -2.0)}),
            Pnt2::new(0.0f64, 0.0),
            10.0f64,
        );
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
        let mut n = Node::new(NodeState::Empty, Pnt2::new(0.0f64, 0.0), 8.0f64);
        n.insert(Entry { object: 1u, position: Pnt2::new(1.0, 2.0) });
        n.insert(Entry { object: 1u, position: Pnt2::new(2.0, -3.0) });
        match n.state.unwrap() {
            NodeState::Branch(_) => (),
            _ => panic!("node is no branch"),
        }
    }
}
