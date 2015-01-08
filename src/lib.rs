//! A spatial tree library

#![warn(missing_docs)]

#[macro_use]
extern crate nalgebra;

#[cfg(test)]
extern crate quickcheck;
#[cfg(test)]
extern crate test;

pub mod ntree;
pub mod tree;
pub mod util;
pub mod partition;
