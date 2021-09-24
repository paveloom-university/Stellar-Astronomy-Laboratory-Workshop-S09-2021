//! Navarro-Frenk-White potential

use crate::orbit::F;

/// Calculate the value of the Navarro-Frenk-White potential
pub fn phi(r: F, z: F, m: F, a: F) -> F {
    let dist = (r.powi(2) + z.powi(2)).sqrt();
    -m / dist * F::ln(1.0 + dist / a)
}

/// Calculate the value of the R derivative of the Navarro-Frenk-White potential
pub fn phi_dr(r: F, z: F, m: F, a: F) -> F {
    let k_1 = r.powi(2) + z.powi(2);
    let k_2 = a + k_1;
    2.0 * m * r * (k_2 * F::ln(k_2 / a) - k_1) / k_1.powi(2) / k_2
}

/// Calculate the value of the Z derivative of the Navarro-Frenk-White potential
pub fn phi_dz(r: F, z: F, m: F, a: F) -> F {
    let k_1 = r.powi(2) + z.powi(2);
    let k_2 = a + k_1;
    2.0 * m * z * (k_2 * F::ln(k_2 / a) - k_1) / k_1.powi(2) / k_2
}
