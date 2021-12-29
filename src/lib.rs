//! A spatial tree library

#![warn(missing_docs)]

#[macro_use]
extern crate itertools;

pub use self::{
	traits::*,
	pure_tree::PureTree,
	data_tree::Tree,
	error::ConstructionError,
};

pub mod partition;
pub mod pure_tree;
pub mod data_tree;
pub mod traits;
pub mod iter;
pub mod error;
