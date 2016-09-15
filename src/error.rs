//! Error types.

/// Errors during tree construction
#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum ConstructionError {
    /// An object lays outside the tree domain
    ObjectOutsidePartition
}
