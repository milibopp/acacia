//! TODO: write crate documentation

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

#[phase(plugin, link)] extern crate log;

pub mod tree;
pub mod util;
