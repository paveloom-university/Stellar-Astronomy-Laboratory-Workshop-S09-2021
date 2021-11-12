//! Constants and methods that are used to correct the
//! coordinates for the solar motion or the Local Standard of Rest

use crate::F;

/// Galactocentric distance of the Local Standard
/// of Rest around the Galactic center $\[ \text{kpc} \]$
const LSR_R: F = 8.3;
/// Galactocentric linear velocity of the Local Standard
/// of Rest around the Galactic center $\[ \text{km} \\, \text{s}^{-1} \]$
const LSR_V: F = 244.0;
/// Height of the Sun above the Galactic plane $\[ \text{kpc} \]$
const SUN_H: F = 16.0 / 1000.0;
/// The Sun’s peculiar velocity U with respect
/// to the Local Standard of Rest $\[ \text{km} \\, \text{s}^{-1} \]$
const SUN_U: F = 11.1;
/// The Sun’s peculiar velocity V with respect
/// to the Local Standard of Rest $\[ \text{km} \\, \text{s}^{-1} \]$
const SUN_V: F = 12.2;
/// The Sun’s peculiar velocity W with respect
/// to the Local Standard of Rest $\[ \text{km} \\, \text{s}^{-1} \]$
const SUN_W: F = 7.3;

/// Adds correction methods to float types
pub trait Corrections {
    /// A float type. Should be the same as Self.
    type Output;
    /// Convert X component of the radius vector from the
    /// Heliocentric Cartesian system to Galactic Cylindrical system
    fn to_gca_x(&self) -> Self::Output;
    /// Convert Z component of the radius vector from the
    /// Heliocentric Cartesian system to Galactic Cylindrical system
    fn to_gca_z(&self) -> Self::Output;
    /// Convert U component of the velocity vector from the
    /// Heliocentric Cartesian system to Galactic Cylindrical system
    fn to_gca_u(&self) -> Self::Output;
    /// Convert V component of the velocity vector from the
    /// Heliocentric Cartesian system to Galactic Cylindrical system
    fn to_gca_v(&self) -> Self::Output;
    /// Convert W component of the velocity vector from the
    /// Heliocentric Cartesian system to Galactic Cylindrical system
    fn to_gca_w(&self) -> Self::Output;
}

impl Corrections for F {
    type Output = F;
    fn to_gca_x(&self) -> Self::Output {
        LSR_R - self
    }
    fn to_gca_z(&self) -> Self::Output {
        self + SUN_H
    }
    fn to_gca_u(&self) -> Self::Output {
        self + SUN_U
    }
    fn to_gca_v(&self) -> Self::Output {
        self + SUN_V + LSR_V
    }
    fn to_gca_w(&self) -> Self::Output {
        self + SUN_W
    }
}
