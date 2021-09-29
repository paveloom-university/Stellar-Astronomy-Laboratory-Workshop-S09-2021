//! Miyamoto & Nagai potential

use crate::orbit::F;

/// $ M $ parameter of the Miyamoto & Nagai potential
const M: F = 2798.0;
/// $ a $ parameter of the Miyamoto & Nagai potential
const A: F = 4.40;
/// $ b $ parameter of the Miyamoto & Nagai potential
const B: F = 0.3084;

/// Calculate the value of the Miyamoto & Nagai potential
pub fn phi(r: F, z: F) -> F {
    -M / (r.powi(2) + (A + (z.powi(2) + B.powi(2)).sqrt()).powi(2)).sqrt()
}

/// Calculate the value of the R derivative of the Miyamoto & Nagai potential
pub fn phi_dr(r: F, z: F) -> F {
    M * r / (r.powi(2) + (A + (z.powi(2) + B.powi(2)).sqrt()).powi(2)).powf(1.5)
}

/// Calculate the value of the Z derivative of the Miyamoto & Nagai potential
pub fn phi_dz(r: F, z: F) -> F {
    let k_1 = (z.powi(2) + B.powi(2)).sqrt();
    let k_2 = A + k_1;
    M * z * k_2 / k_1 / (k_2.powi(2) + r.powi(2)).powf(1.5)
}
