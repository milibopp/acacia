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
    /// prescribes, and `combine` data subsequently modifying an `accumulator`.
    ///
    /// If an empty or leaf node is encountered, `combine` is called on the
    /// accumulator and its associated data. For a branching node `recurse` is
    /// called on it, to determine whether its subnodes should be inspected more
    /// closely. If so, the function recurses on each subnode, otherwise it acts
    /// on it as if it were not a branch.
    ///
    /// # Parameters
    ///
    /// - The `accumulator` is a mutable reference to some data, that is
    ///   modified during the query to collect.
    /// - At each node the tree is only recursed further, iff `recurse(&node)`.
    /// - `combine` is called for every node whose data is to be considered.
    fn query_data_mut<T>(&self, accumulator: &mut T, recurse: |&Self| -> bool, combine: |&mut T, &D|);

    /// Compute a query on the associated data using an initial state
    ///
    /// This is very similar to the accumulator variant and indeed has a default
    /// implementation using it. The difference is, that an initial value is
    /// moved into the function and used to initialize the accumulator. Its
    /// final state is the method's return value.
    fn query_data<T>(&self, initial: T, recurse: |&Self| -> bool, combine: |&T, &D| -> T) -> T {
        let mut acc = initial;
        self.query_data_mut(
            &mut acc,
            |node| recurse(node),
            |a, d| {*a = combine(a, d);}
        );
        acc
    }
}


/// A tree that allows recursive queries on its objects
///
/// Closures are used to determine the recursion behavior and what is to be
/// computed.
///
/// TODO: do not mention DataQuery methods in documentation, to make this
/// self-contained
///
/// # Type parameters
///
/// - `O` is the type of the objects stored in the tree.
pub trait ObjectQuery<O> {

    /// Compute a query on the objects using an accumulator
    ///
    /// This method walks through the tree similar to `query_data_mut`. However,
    /// the `combine` closure is only invoked, when a leaf node is encountered.
    /// `recurse` determines if it recurses into branch nodes more deeply. Empty
    /// nodes are ignored.
    fn query_objects_mut<T>(&self, accumulator: &mut T, recurse: |&Self| -> bool, combine: |&mut T, &O|);

    /// Compute a query on the objects using an initial state
    ///
    /// This relates to `query_objects_mut` in the same way `query_data` relates
    /// to `query_data_mut`.
    fn query_objects<T>(&self, initial: T, recurse: |&Self| -> bool, combine: |&T, &O| -> T) -> T {
        let mut acc = initial;
        self.query_objects_mut(
            &mut acc,
            |node| recurse(node),
            |a, o| {*a = combine(a, o);}
        );
        acc
    }
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
