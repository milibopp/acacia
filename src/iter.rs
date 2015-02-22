//! Generic tree iterators.

use traits::{Node, AssociatedData};


/// An iterator over the objects in a tree.
pub struct Iter<'a, T: 'a> {
    nodes: Vec<(usize, &'a T)>,
}

impl<'a, T> Iter<'a, T> {
    /// Create a new iterator.
    pub fn new(tree: &'a T) -> Iter<'a, T> {
        Iter { nodes: vec![(0, tree)] }
    }
}

impl<'a, T> Iterator for Iter<'a, T>
    where T: Node<Container = Vec<T>>,
{
    type Item = &'a <T as Node>::Object;

    fn next(&mut self) -> Option<&'a <T as Node>::Object> {
        use traits::NodeState::*;
        match self.nodes.pop() {
            None => None,
            Some((count, node)) => match (count, node.state()) {
                (_, Empty) => self.next(),
                (_, Leaf(obj)) => Some(obj),
                (n, Branch(vec)) => {
                    if let Some(child) = vec.get(n) {
                        self.nodes.push((n + 1, node));
                        self.nodes.push((0, child));
                    }
                    self.next()
                },
            }
        }
    }
}


/// An iterator over the objects in a tree, that only recurses as deep as
/// specified by some predicate.
pub struct RecurseObjects<'a, T: 'a, R> {
    nodes: Vec<(usize, &'a T)>,
    recurse: R,
}

impl<'a, T, R> RecurseObjects<'a, T, R> {
    /// Create a new iterator.
    pub fn new(tree: &'a T, recurse: R) -> RecurseObjects<'a, T, R> {
        RecurseObjects { nodes: vec![(0, tree)], recurse: recurse }
    }
}

impl<'a, T, R> Iterator for RecurseObjects<'a, T, R>
    where T: Node<Container = Vec<T>>,
          R: FnMut(&T) -> bool,
{
    type Item = &'a <T as Node>::Object;

    fn next(&mut self) -> Option<&'a <T as Node>::Object> {
        use traits::NodeState::*;
        match self.nodes.pop() {
            None => None,
            Some((count, node)) => match (count, node.state()) {
                (_, Empty) => self.next(),
                (_, Leaf(obj)) => Some(obj),
                (n, Branch(vec)) => {
                    if let Some(child) = vec.get(n) {
                        if (self.recurse)(node) {
                            self.nodes.push((n + 1, node));
                            self.nodes.push((0, child));
                        }
                    }
                    self.next()
                },
            }
        }
    }
}


/// An iterator over the objects in a tree, that only recurses as deep as
/// specified by some predicate.
pub struct RecurseData<'a, T: 'a, R> {
    nodes: Vec<(Option<usize>, &'a T)>,
    recurse: R,
}

impl<'a, T, R> RecurseData<'a, T, R> {
    /// Create a new iterator.
    pub fn new(tree: &'a T, recurse: R) -> RecurseData<'a, T, R> {
        RecurseData { nodes: vec![(None, tree)], recurse: recurse }
    }
}

impl<'a, T, R> Iterator for RecurseData<'a, T, R>
    where T: Node<Container = Vec<T>> + AssociatedData,
          <T as Node>::Object: 'a,
          R: FnMut(&T) -> bool,
{
    type Item = &'a <T as AssociatedData>::Data;

    fn next(&mut self) -> Option<&'a <T as AssociatedData>::Data> {
        use traits::NodeState::*;
        match self.nodes.pop() {
            None => None,
            Some((count, node)) => match (count, node.state()) {
                (None, Branch(_)) => {
                    if (self.recurse)(node) {
                        self.nodes.push((Some(0), node));
                        self.next()
                    } else {
                        Some(node.data())
                    }
                }
                (Some(n), Branch(vec)) => {
                    if let Some(child) = vec.get(n) {
                        self.nodes.push((Some(n + 1), node));
                        self.nodes.push((None, child));
                    }
                    self.next()
                },
                _ => Some(node.data()),
            }
        }
    }
}


#[cfg(test)]
mod test {
    use nalgebra::{Pnt2, Orig};

    use partition::Ncube;
    use traits::{Positioned, Node};
    use pure_tree::PureTree;
    use super::*;

    #[test]
    fn iter_pure_tree() {
        let tree = PureTree::new(
            vec![
                Positioned { object: 1, position: Pnt2::new(-0.1, 1.0) },
                Positioned { object: 2, position: Pnt2::new(0.5, -0.3) },
            ].into_iter(),
            Ncube::new(Orig::orig(), 2.0)
        );
        let all: Vec<_> = Iter::new(&tree).collect();
        assert_eq!(all.len(), 2);
        assert!(all[0].object == 1 || all[1].object == 1);
    }

    #[test]
    fn recurse_objects_pure_tree() {
        let tree = PureTree::new(
            vec![
                Positioned { object: 1i32, position: Pnt2::new(-0.1, 0.8) },
                Positioned { object: 2, position: Pnt2::new(-0.2, 0.7) },
                Positioned { object: 3, position: Pnt2::new(0.5, -0.3) },
            ].into_iter(),
            Ncube::new(Orig::orig(), 2.0)
        );
        let all: Vec<_> =
            RecurseObjects::new(
                &tree,
                |n: &PureTree<Ncube<_, f64>, _>|
                    n.partition().width() > 1.5
            )
            .collect();
        assert_eq!(all.len(), 1);
        assert_eq!(all[0].object, 3);
    }
}
