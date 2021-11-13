//! This module provides the models of the Galactic potential

mod m1;
mod model;

pub use m1::M1;
pub use model::Model;

/// An array of models: indices of the items are the same as the indices of the models
pub const MODELS: &[&dyn Model] = &[&M1 {}];
