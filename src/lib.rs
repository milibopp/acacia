//! A spatial tree library

#![feature(phase)]
#![feature(unboxed_closures)]
#![warn(missing_docs)]

extern crate nalgebra;

#[cfg(test)]
extern crate quickcheck;
#[cfg(test)]
#[phase(plugin)]
extern crate quickcheck_macros;
#[cfg(test)]
#[phase(plugin)]
extern crate nalgebra;
#[cfg(test)]
extern crate test;

pub mod ntree;
pub mod tree;
pub mod util;
