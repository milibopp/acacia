//! A spatial tree library

#![feature(phase)]

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

pub mod tree;
pub mod util;
