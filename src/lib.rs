//! TODO: write crate documentation

#![feature(phase)]

extern crate nalgebra;

#[cfg(test)]
extern crate quickcheck;
#[cfg(test)]
#[phase(plugin)]
extern crate quickcheck_macros;

pub mod tree;
