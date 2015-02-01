//! A spatial tree library

#![warn(missing_docs)]
#![feature(core, hash)]

extern crate nalgebra;
#[macro_use]
extern crate itertools;

#[cfg(any(test, feature = "arbitrary"))]
extern crate quickcheck;
#[cfg(test)]
extern crate test;

pub mod ntree;
pub mod tree;
pub mod partition;
