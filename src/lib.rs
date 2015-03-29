//! A spatial tree library

#![warn(missing_docs)]
#![feature(core)]
#![cfg_attr(test, feature(test))]

extern crate nalgebra;
#[macro_use]
extern crate itertools;
extern crate rand;

#[cfg(any(test, feature = "arbitrary"))]
extern crate quickcheck;
#[cfg(test)]
extern crate test;

pub use traits::*;
pub use pure_tree::PureTree;
pub use data_tree::Tree;

pub mod partition;
pub mod pure_tree;
pub mod data_tree;
pub mod traits;
pub mod iter;
