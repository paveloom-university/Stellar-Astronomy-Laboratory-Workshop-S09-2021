//! Miyamoto & Nagai potential

use crate::orbit::F;

/// Calculate the value of the Miyamoto & Nagai potential
#[must_use]
pub fn phi(r: F, z: F, m: F, a: F, b: F) -> F {
    -m / (r.powi(2) + (a + (z.powi(2) + b.powi(2)).sqrt()).powi(2)).sqrt()
}

/// Calculate the value of the R derivative of the Miyamoto & Nagai potential
#[must_use]
pub fn phi_dr(r: F, z: F, m: F, a: F, b: F) -> F {
    m * r / (r.powi(2) + (a + (z.powi(2) + b.powi(2)).sqrt()).powi(2)).powf(1.5)
}

/// Calculate the value of the Z derivative of the Miyamoto & Nagai potential
#[must_use]
pub fn phi_dz(r: F, z: F, m: F, a: F, b: F) -> F {
    let k_1 = (z.powi(2) + b.powi(2)).sqrt();
    let k_2 = a + k_1;
    m * z * k_2 / k_1 / (k_2.powi(2) + r.powi(2)).powf(1.5)
}
