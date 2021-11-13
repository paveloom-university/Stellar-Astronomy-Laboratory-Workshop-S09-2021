use crate::F;

pub trait Model {
    /// Calculate the value of the Galactic potential $ \Phi(R, Z) $
    /// $\[ 100 \\, \text{km}^2 \\, \text{s}^{-2} \]$
    fn phi(r: F, z: F) -> F;
    /// Calculate the value of $ \partial \Phi(R, Z) / \partial R $
    /// $\[ \text{kpc} \\, \text{Myr}^{-2} \]$
    fn phi_dr(r: F, z: F) -> F;
    /// Calculate the value of $ \partial \Phi(R, Z) / \partial Z $
    /// $\[ \text{kpc} \\, \text{Myr}^{-2} \]$
    fn phi_dz(r: F, z: F) -> F;
}
