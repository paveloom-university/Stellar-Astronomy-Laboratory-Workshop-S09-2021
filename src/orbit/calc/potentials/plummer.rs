//! Plummer potential

use crate::orbit::F;

/// $ M $ parameter of the Plummer potential
const M: F = 443.0;
/// $ b $ parameter of the Plummer potential
const B: F = 0.2672;

/// Calculate the value of the Plummer potential
pub fn phi(r: F, z: F) -> F {
    -M / (r.powi(2) + z.powi(2) + B.powi(2)).sqrt()
}

/// Calculate the value of the R derivative of the Plummer potential
pub fn phi_dr(r: F, z: F) -> F {
    M * r / (r.powi(2) + z.powi(2) + B.powi(2)).powf(1.5)
}

/// Calculate the value of the Z derivative of the Plummer potential
pub fn phi_dz(r: F, z: F) -> F {
    M * z / (r.powi(2) + z.powi(2) + B.powi(2)).powf(1.5)
}
