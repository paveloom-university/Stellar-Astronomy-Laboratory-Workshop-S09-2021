//! This module provides the models of the Galactic potential
//!
//! There are two groups of models provided, depending on how they represent
//! the Galactic potential:
//! 1. A sum of potentials: bulge + disk + halo;
//! 2. A sum of potentials: bulge + thin disk + thick disk + halo.
//!
//! These are combinations of the available [potentials](crate::orbit::calc::potentials).
//!
//! # Group 1
//!
//! | # | Bulge | Disk | Halo | Values of the parameters    |
//! | - | :---: | :--: | :--: | :-------------------------: |
//! | 1 | P     | MN   | NFW  | Bajkova, Bobylev (2020, v1) |
//!
//! # Group 2
//!
//! | # | Bulge | Thin Disk  | Thick Disk  | Halo | Values of the parameters                         |
//! | - | :---: | :--------: | :---------: | :--: | :----------------------------------------------: |
//! | 2 | P     | MN         | MN          | NFW  | Pouliasis et al. (2017, model I) + Eilers (2018) |

mod m1;
mod m2;
mod model;

pub use m1::M1;
pub use m2::M2;
pub use model::Model;

/// An array of models: indices of the items are the same as the indices of the models
pub const MODELS: &[&dyn Model] = &[&M1 {}, &M2 {}];
