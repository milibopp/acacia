//! Common abstractions for all trees

/// The state of a node
///
/// A node may either be empty, a leaf with exactly one object or a branch with
/// a list of other nodes.
///
/// TODO: abstract over this, so that it works different structural
/// representations
pub enum NodeState<O, N> {
    Empty,
    Leaf(O),
    Branch(Vec<N>),
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
    ///   recursion.
    fn query_data(&self, recurse: |&Self| -> bool, f: |&D|);
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
    /// - `f` is the callback function.
    fn query_objects(&self, recurse: |&Self| -> bool, f: |&O|);
}


/// A tree node
///
/// This is part of the essential features of a tree. Note that both a whole
/// tree and its constituents implement this.
pub trait Node<P, N, O> {
    fn state(&self) -> &NodeState<O, Self>;
    fn center(&self) -> &P;
    fn width(&self) -> &N;
}


/// A tree with associated data
pub trait AssociatedData<D> {
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
pub trait PureTree<P, N, O>: ObjectQuery<O> + Node<P, N, O> {}


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
pub trait Tree<P, N, O, D>: DataQuery<D> + AssociatedData<D> + PureTree<P, N, O> {}
