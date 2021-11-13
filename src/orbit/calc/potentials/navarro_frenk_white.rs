//! Navarro-Frenk-White potential

use crate::orbit::F;

/// Calculate the value of the Navarro-Frenk-White potential
#[must_use]
pub fn phi(r: F, z: F, m: F, a: F) -> F {
    let dist = (r.powi(2) + z.powi(2)).sqrt();
    -m / dist * F::ln(1.0 + dist / a)
}

/// Calculate the value of the R derivative of the Navarro-Frenk-White potential
#[must_use]
pub fn phi_dr(r: F, z: F, m: F, a: F) -> F {
    let sq_sum = r.powi(2) + z.powi(2);
    let dist = sq_sum.sqrt();
    let k_1 = dist / a + 1.0;
    m * r * F::ln(k_1) / sq_sum.powf(1.5) - m * r / a / sq_sum / k_1
}

/// Calculate the value of the Z derivative of the Navarro-Frenk-White potential
#[must_use]
pub fn phi_dz(r: F, z: F, m: F, a: F) -> F {
    let sq_sum = r.powi(2) + z.powi(2);
    let dist = sq_sum.sqrt();
    let k_1 = dist / a + 1.0;
    m * z * F::ln(k_1) / sq_sum.powf(1.5) - m * z / a / sq_sum / k_1
}
