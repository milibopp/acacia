//! A spatial tree library

#![warn(missing_docs)]
#![allow(unused_features)]
#![feature(core, hash, rand, test)]

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
