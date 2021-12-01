//! This module provides the [`M1`] model

use super::Model;
use crate::orbit::calc::potentials::{miyamoto_nagai, navarro_frenk_white, plummer};
use crate::F;

// Bulge

/// $ M $ parameter of the Plummer potential
const M_B: F = 443.0;
/// $ b $ parameter of the Plummer potential
const B_B: F = 0.2672;

// Disk

/// $ M $ parameter of the Miyamoto & Nagai potential
const M_D: F = 2798.0;
/// $ a $ parameter of the Miyamoto & Nagai potential
const A_D: F = 4.40;
/// $ b $ parameter of the Miyamoto & Nagai potential
const B_D: F = 0.3084;

// Halo

/// $ M $ parameter of the Navarro-Frenk-White potential
const M_H: F = 12474.0;
/// $ a $ parameter of the Navarro-Frenk-White potential
const A_H: F = 7.7;

/// This model uses the Plummer potential for bulge, the Miyamoto & Nagai potential
/// for disk, and the Navarro-Frenk-White potential for halo. The values of the parameters
/// are taken from Bajkova, Bobylev (2020, v1).
pub struct M1 {}

impl Model for M1 {
    fn phi(&self, r: F, z: F) -> F {
        plummer::phi(r, z, M_B, B_B)
            + miyamoto_nagai::phi(r, z, M_D, A_D, B_D)
            + navarro_frenk_white::phi(r, z, M_H, A_H)
    }
    fn phi_dr(&self, r: F, z: F) -> F {
        plummer::phi_dr(r, z, M_B, B_B)
            + miyamoto_nagai::phi_dr(r, z, M_D, A_D, B_D)
            + navarro_frenk_white::phi_dr(r, z, M_H, A_H)
    }
    fn phi_dz(&self, r: F, z: F) -> F {
        plummer::phi_dz(r, z, M_B, B_B)
            + miyamoto_nagai::phi_dz(r, z, M_D, A_D, B_D)
            + navarro_frenk_white::phi_dz(r, z, M_H, A_H)
    }
}
