//! This module provides the [`M2`] model

use super::Model;
use crate::orbit::calc::potentials::{miyamoto_nagai, navarro_frenk_white, plummer};
use crate::{F, PI};

// Bulge

/// $ M $ parameter of the Plummer potential
const M_B: F = 460.0;
/// $ b $ parameter of the Plummer potential
const B_B: F = 0.3;

// Thin disk

/// $ M $ parameter of the Miyamoto & Nagai potential
const M_THIN_D: F = 1700.0;
/// $ a $ parameter of the Miyamoto & Nagai potential
const A_THIN_D: F = 5.3;
/// $ b $ parameter of the Miyamoto & Nagai potential
const B_THIN_D: F = 0.25;

// Thick disk

/// $ M $ parameter of the Miyamoto & Nagai potential
const M_THICK_D: F = 1700.0;
/// $ a $ parameter of the Miyamoto & Nagai potential
const A_THICK_D: F = 2.6;
/// $ b $ parameter of the Miyamoto & Nagai potential
const B_THICK_D: F = 0.8;

// Halo

/// $ M $ parameter of the Navarro-Frenk-White potential
const M_H: F = 4.0 * PI * (1.06 * 1e7) * (14.8 * 14.8 * 14.8) / (2.325 * 1e7);
/// $ a $ parameter of the Navarro-Frenk-White potential
const A_H: F = 14.8;

/// This model uses the Plummer potential for bulge, the Miyamoto & Nagai potential
/// for thin and thick disks, and the Navarro-Frenk-White potential for halo. The
/// values of the parameters for the bulge, the thin disk and the thick disk are
/// taken from Pouliasis et al. (2017, model I), values of the parameters for the
/// halo are taken from Eilers (2018).
pub struct M2 {}

impl Model for M2 {
    fn phi(&self, r: F, z: F) -> F {
        plummer::phi(r, z, M_B, B_B)
            + miyamoto_nagai::phi(r, z, M_THIN_D, A_THIN_D, B_THIN_D)
            + miyamoto_nagai::phi(r, z, M_THICK_D, A_THICK_D, B_THICK_D)
            + navarro_frenk_white::phi(r, z, M_H, A_H)
    }
    fn phi_dr(&self, r: F, z: F) -> F {
        plummer::phi_dr(r, z, M_B, B_B)
            + miyamoto_nagai::phi_dr(r, z, M_THIN_D, A_THIN_D, B_THIN_D)
            + miyamoto_nagai::phi_dr(r, z, M_THICK_D, A_THICK_D, B_THICK_D)
            + navarro_frenk_white::phi_dr(r, z, M_H, A_H)
    }
    fn phi_dz(&self, r: F, z: F) -> F {
        plummer::phi_dz(r, z, M_B, B_B)
            + miyamoto_nagai::phi_dz(r, z, M_THIN_D, A_THIN_D, B_THIN_D)
            + miyamoto_nagai::phi_dz(r, z, M_THICK_D, A_THICK_D, B_THICK_D)
            + navarro_frenk_white::phi_dz(r, z, M_H, A_H)
    }
}
