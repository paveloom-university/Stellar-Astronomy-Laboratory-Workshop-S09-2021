//! This module provides the models of the Galactic potential
//!
//! There are two types of models provided, depending on how they represent
//! the Galactic potential:
//! 1. A sum of potentials: bulge + disk + halo;
//! 2. A sum of potentials: bulge + thin disk + thick disk + halo.
//!
//! These are combinations of the available [potentials](crate::orbit::calc::potentials).
//!
//! # Group 1
//!
//! | # | Bulge | Disk | Halo |
//! | - | :---: | :--: | :--: |
//! | 1 | P     | MN   | NFW  |
//! | 2 | MN    | MN   | NFW  |
//!
//! # Group 2
//!
//! | # | Bulge | Thin Disk | Thick Disk | Halo |
//! | - | :---: | :-------: | :--------: | :--: |
//! | 3 | P     | ?         | ?          | NFW  |

mod m1;
mod model;

pub use m1::M1;
pub use model::Model;

/// An array of models: indices of the items are the same as the indices of the models
pub const MODELS: &[&dyn Model] = &[&M1 {}];
