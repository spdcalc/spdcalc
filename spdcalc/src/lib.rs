pub extern crate dimensioned as dim;
pub extern crate nalgebra as na;

pub mod types;
pub use types::*;
pub mod constants;
pub use constants::*;

pub mod math;

pub mod photon;
pub mod spd;

pub mod crystal;
pub mod junk;
pub mod utils;
pub use crystal::Crystal;
pub mod plotting;

#[allow(unused_imports)]
#[cfg(test)]
#[macro_use]
extern crate float_cmp;
