//! Navarro-Frenk-White potential

use crate::orbit::F;

/// $ M $ parameter of the Navarro-Frenk-White potential
const M: F = 12474.0;
/// $ A $ parameter of the Navarro-Frenk-White potential
const A: F = 7.7;

/// Calculate the value of the Navarro-Frenk-White potential
#[must_use]
pub fn phi(r: F, z: F) -> F {
    let dist = (r.powi(2) + z.powi(2)).sqrt();
    -M / dist * F::ln(1.0 + dist / A)
}

/// Calculate the value of the R derivative of the Navarro-Frenk-White potential
#[must_use]
pub fn phi_dr(r: F, z: F) -> F {
    let sq_sum = r.powi(2) + z.powi(2);
    let dist = sq_sum.sqrt();
    let k_1 = dist / A + 1.0;
    M * r * F::ln(k_1) / sq_sum.powf(1.5) - M * r / A / sq_sum / k_1
}

/// Calculate the value of the Z derivative of the Navarro-Frenk-White potential
#[must_use]
pub fn phi_dz(r: F, z: F) -> F {
    let sq_sum = r.powi(2) + z.powi(2);
    let dist = sq_sum.sqrt();
    let k_1 = dist / A + 1.0;
    M * z * F::ln(k_1) / sq_sum.powf(1.5) - M * z / A / sq_sum / k_1
}
