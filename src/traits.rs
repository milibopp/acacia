//! Common abstractions for all trees

use iter::{RecurseObjects, RecurseData};


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


/// A tree node
///
/// This is part of the essential features of a tree. Note that both a whole
/// tree and its constituents implement this.
pub trait Node {

    /// Type of spatial partitioning scheme
    type Partition;

    /// The type of object stored
    type Object;

    /// Type of container used to store subnodes
    type Container;

    /// The state of the node
    fn state(&self) -> NodeState<&Self::Object, &Self::Container>;

    /// The partitioning scheme
    fn partition(&self) -> Self::Partition;
}


/// A tree with associated data
pub trait AssociatedData {

    /// Type of the associated data
    type Data;

    /// Data associated to the node
    fn data(&self) -> &Self::Data;
}


/// A tree that allows recursive queries on its objects. A closure is used to
/// determine the recursion behavior.
///
/// This is an extension trait of `Node`.
pub trait ObjectQuery: Node {
    /// Iterate over objects through all nodes.
    ///
    /// This method yields an iterator that walks recursively through the tree,
    /// as deep as `recurse` prescribes.
    ///
    /// Empty nodes and branch nodes with `!recurse(&node)` are omitted, whereas
    /// the iterator considers every object in a leaf node.
    ///
    /// # Parameters
    ///
    /// - At each branching node the tree is only recursed further, iff
    ///   `recurse(&node)`.
    fn query_objects<'a, R>(&'a self, recurse: R) -> RecurseObjects<'a, Self, R>
        where R: Fn(&Self) -> bool;
}

impl<T> ObjectQuery for T
    where T: Node<Container = Vec<T>>,
{
    fn query_objects<'a, R>(&'a self, recurse: R) -> RecurseObjects<'a, Self, R>
        where R: Fn(&T) -> bool
    {
        RecurseObjects::new(self, recurse)
    }
}


/// A tree which allows recursive queries on its associated data
///
/// Closures are used to determine the recursion behavior and what is to be
/// computed.
///
/// # Type parameters
///
/// - `D` is the type of the associated data.
pub trait DataQuery: AssociatedData {

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
    fn query_data<'a, R>(&'a self, recurse: R) -> RecurseData<'a, Self, R>
        where R: Fn(&Self) -> bool;
}

impl<T> DataQuery for T
    where T: Node<Container=Vec<T>> + AssociatedData,
{
    fn query_data<'a, R>(&'a self, recurse: R) -> RecurseData<'a, T, R>
        where R: Fn(&Self) -> bool
    {
        RecurseData::new(self, recurse)
    }
}


/// A type that has a notion of a position
pub trait Position {
    /// The underlying point type
    type Point;

    /// The position
    fn position(&self) -> Self::Point;
}

impl<'a, O> Position for &'a O
    where O: Position
{
    type Point = <O as Position>::Point;

    fn position(&self) -> <O as Position>::Point {
        // Explicitly dereference here to avoid infinite recursion
        (*self).position()
    }
}


/// A positioned object
///
/// This is the most simple generic implementation of Position and serves as a
/// wrapper for types that do not have a notion of a position themselves. It
/// equips these with an additional generic position as an attribute.
#[derive(Clone)]
pub struct Positioned<O, P> {
    /// The object wrapped in this type
    pub object: O,

    /// The position stored along with it
    pub position: P,
}

impl<O, P> Position for Positioned<O, P>
    where P: Copy
{
    type Point = P;

    fn position(&self) -> P {
        self.position
    }
}


#[cfg(test)]
mod test {

    use super::{Position, Positioned};

    #[test]
    fn positioned_position() {
        assert_eq!(Positioned { object: (), position: 14 }.position(), 14);
    }

    #[test]
    fn positionable_by_ref() {
        fn twice_pos<O: Position<Point=i32>>(obj: O) -> i32 {
            2 * obj.position()
        }
        let obj = Positioned { object: (), position: 77 };
        assert_eq!(twice_pos(&obj), 154);
    }
}
