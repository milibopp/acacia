//! Common abstractions for all trees

/// The state of a node
///
/// A node may either be empty, a leaf with exactly one object or a branch with
/// a list of other nodes. This enum encodes these states and the data
/// associated with each of them.
///
/// # Type parameters
///
/// - `O` is the type of object stored in the tree structure.
/// - `C` is a collection of nodes in a branch.
pub enum NodeState<O, C> {

    /// An empty node does not contain any object
    Empty,

    /// A leaf node contains exactly one object
    Leaf(O),

    /// A branch node contains a collection of nodes
    Branch(C),
}


/// A tree which allows recursive queries on its associated data
///
/// Closures are used to determine the recursion behavior and what is to be
/// computed.
///
/// # Type parameters
///
/// - `D` is the type of the associated data.
pub trait DataQuery<D> {

    /// Compute a query on the associated data using a mutable accumulator
    ///
    /// This method walks recursively through the tree, as deep as `recurse`
    /// prescribes, and calls a function on the associated data of each node
    /// encountered.
    ///
    /// If an empty or leaf node is encountered, the function is called on its
    /// associated data. For a branching node `recurse` checks, whether its
    /// subnodes should be inspected more closely. If so, this method recurses
    /// on each subnode, otherwise it simply calls the function on its
    /// associated data.
    ///
    /// # Parameters
    ///
    /// - At each branching node the tree is only recursed further, iff
    ///   `recurse(&node)`.
    /// - `f` is called on the associated data of every node reached by the
    ///   recursion. This may mutably borrow its environment, which is currently
    ///   the only way to obtain a result from this function.
    fn query_data<R, F>(&self, recurse: &R, f: &mut F)
        where R: Fn(&Self) -> bool,
              F: FnMut(&D);
}


/// A tree that allows recursive queries on its objects
///
/// Closures are used to determine the recursion behavior and what is to be
/// computed.
///
/// # Type parameters
///
/// - `O` is the type of the objects stored in the tree.
pub trait ObjectQuery<O> {

    /// Compute a query on the objects using an accumulator
    ///
    /// This method walks recursively through the tree, as deep as `recurse`
    /// prescribes, and calls a function on each object encountered.
    ///
    /// Empty nodes and branch nodes with `!recurse(&node)` are ignored, whereas
    /// the callback is called on every object in a leaf node.
    ///
    /// # Parameters
    ///
    /// - At each branching node the tree is only recursed further, iff
    ///   `recurse(&node)`.
    /// - `f` is the callback function. This may mutably borrow its environment,
    ///   which is currently the only way to obtain a result from this function.
    fn query_objects<R, F>(&self, recurse: &R, f: &mut F)
        where R: Fn(&Self) -> bool,
              F: FnMut(&O);
}


/// A tree node
///
/// This is part of the essential features of a tree. Note that both a whole
/// tree and its constituents implement this.
pub trait Node<P, N, O, C> {

    /// The state of the node
    fn state(&self) -> &NodeState<O, C>;

    /// The center point of the node
    fn center(&self) -> &P;

    /// The width of the node
    fn width(&self) -> &N;
}


/// A tree with associated data
pub trait AssociatedData<D> {

    /// Data associated to the node
    fn data(&self) -> &D;
}


/// A pure spatial tree
///
/// This trait wraps up the properties of a tree that contains only spatial
/// information, but no associated data.
///
/// # Type parameters
///
/// - `P` is the kind of point used to position objects and nodes spatially.
/// - `N` is the scalar of the vector space of points.
/// - The tree stores objects of type `O`. These objects need to have some
///   notion of a position.
pub trait PureTree<P, N, O, C>: ObjectQuery<O> + Node<P, N, O, C> {}


/// A spatial tree with associated data
///
/// This trait wraps up the properties of a tree that contains associated data.
///
/// # Type parameters
///
/// - `P` is the kind of point used to position objects and nodes spatially.
/// - `N` is the scalar of the vector space of points.
/// - The tree stores objects of type `O`. These objects need to have some
///   notion of a position.
/// - `D` is the kind of data associated with each node. This is computed
///   recursively during tree construction.
pub trait Tree<P, N, O, C, D>: DataQuery<D> + AssociatedData<D> + PureTree<P, N, O, C> {}


/// A type that has a notion of a position
pub trait Positionable<P> {

    /// The position
    fn position(&self) -> P;
}

impl<'a, P, O> Positionable<P> for &'a O
    where O: Positionable<P>
{
    fn position(&self) -> P {
        // Explicitly dereference here to avoid infinite recursion
        (*self).position()
    }
}


/// A positioned object
///
/// This is the most simple generic implementation of Positionable and serves as
/// a wrapper for types that do not have a notion of a position themselves. It
/// equips these with an additional generic position as an attribute.
#[derive(Clone)]
pub struct Positioned<O, P> {

    /// The object wrapped in this type
    pub object: O,

    /// The position stored along with it
    pub position: P,
}

impl<O, P> Positionable<P> for Positioned<O, P>
    where P: Copy
{
    fn position(&self) -> P {
        self.position
    }
}


#[cfg(test)]
mod test {

    use super::{Positionable, Positioned};

    #[test]
    fn positioned_position() {
        assert_eq!(Positioned { object: 1u, position: 14i }.position(), 14i);
    }

    #[test]
    fn positionable_by_ref() {
        fn twice_pos<O: Positionable<int>>(obj: O) -> int {
            2 * obj.position()
        }
        let obj = Positioned { object: 1u, position: 77i };
        assert_eq!(twice_pos(&obj), 154i);
    }
}
