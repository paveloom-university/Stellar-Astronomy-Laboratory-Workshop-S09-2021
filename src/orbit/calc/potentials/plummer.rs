//! Plummer potential

use crate::orbit::F;

/// Calculate the value of the Plummer potential
pub fn phi(r: F, z: F, m: F, b: F) -> F {
    -m / (r.powi(2) + z.powi(2) + b.powi(2)).sqrt()
}

/// Calculate the value of the R derivative of the Plummer potential
pub fn phi_dr(r: F, z: F, m: F, b: F) -> F {
    m * r / (r.powi(2) + z.powi(2) + b.powi(2)).powf(1.5)
}

/// Calculate the value of the Z derivative of the Plummer potential
pub fn phi_dz(r: F, z: F, m: F, b: F) -> F {
    m * z / (r.powi(2) + z.powi(2) + b.powi(2)).powf(1.5)
}
