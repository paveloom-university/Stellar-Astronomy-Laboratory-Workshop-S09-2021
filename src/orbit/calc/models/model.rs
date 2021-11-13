//! This module provides the [`Model`] trait

use crate::F;

/// This trait defines that a model needs to provide functions
/// that return the values of Galactic potential and its derivatives.
pub trait Model {
    /// Calculate the value of the Galactic potential $ \Phi(R, Z) $
    /// $\[ 100 \\, \text{km}^2 \\, \text{s}^{-2} \]$
    fn phi(&self, r: F, z: F) -> F;
    /// Calculate the value of $ \partial \Phi(R, Z) / \partial R $
    /// $\[ \text{kpc} \\, \text{Myr}^{-2} \]$
    fn phi_dr(&self, r: F, z: F) -> F;
    /// Calculate the value of $ \partial \Phi(R, Z) / \partial Z $
    /// $\[ \text{kpc} \\, \text{Myr}^{-2} \]$
    fn phi_dz(&self, r: F, z: F) -> F;
}
