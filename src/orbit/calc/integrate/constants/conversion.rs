//! Constants and methods that are used to transfer
//! numbers from one unit of measurement to another

use crate::F;

/// A constant to convert $ \text{kpc} \\, \text{Myr}^{-1} $ to
/// $ \text{km} \\, \text{s}^{-1} $
const KPC_PER_MYR_TO_KM_PER_S: F = 977.812_299_895_122;
/// A constant to convert $ \text{km} \\, \text{s}^{-1} $ to
/// $ \text{kpc} \\, \text{Myr}^{-1} $
const KM_PER_S_TO_KPC_PER_MYR: F = 0.001_022_73;
/// A constant to convert $ 100 \\; \text{km}^2 \\, \text{s}^{-2} $ to
/// $ \text{kpc}^2 \\, \text{Myr}^{-2} $
///
/// This constant is a replacement of the gravitational constant $ G $
/// (see Appendix A in Irrgang et al.(2013)). Since $ G $ can be taken
/// out of any potential, it is also used for potential differentials.
const HUNDREDS_KM_2_PER_S_2_TO_KPC_2_PER_MYR_2: F = 1.045_897_218_694_908e-4;
/// A constant to convert $ \text{Myr} $ to seconds
const MYR_TO_S: F = 3.155_695_2e13;
/// A constant to convert $ \text{km} $ to $ \text{kpc} $
const KM_TO_KPC: F = 3.240_779_289_666_4e-17;

/// Adds conversion methods to float types
pub trait Conversions {
    /// A float type. Should be the same as Self.
    type Output;
    /// Convert $ \text{kpc} \\, \text{Myr}^{-1} $ to
    /// $ \text{km} \\, \text{s}^{-1} $
    fn to_km_per_s(&self) -> Self::Output;
    /// Convert $ \text{km} \\, \text{s}^{-1} $ to
    /// $ \text{kpc} \\, \text{Myr}^{-1} $
    fn to_kpc_per_myr(&self) -> Self::Output;
    /// Convert $ 100 \\; \text{km}^2 \\, \text{s}^{-2} $ to
    /// $ \text{kpc}^2 \\, \text{Myr}^{-2} $
    fn to_kpc_per_myr_2(&self) -> Self::Output;
    /// Convert $ \text{Myr} $ to seconds
    fn to_seconds(&self) -> Self::Output;
    /// Convert $ \text{km} $ to $ \text{kpc} $
    fn to_kpc(&self) -> Self::Output;
}

impl Conversions for F {
    type Output = F;
    fn to_km_per_s(&self) -> Self::Output {
        self * KPC_PER_MYR_TO_KM_PER_S
    }
    fn to_kpc_per_myr(&self) -> Self::Output {
        self * KM_PER_S_TO_KPC_PER_MYR
    }
    fn to_kpc_per_myr_2(&self) -> Self::Output {
        self * HUNDREDS_KM_2_PER_S_2_TO_KPC_2_PER_MYR_2
    }
    fn to_seconds(&self) -> Self::Output {
        self * MYR_TO_S
    }
    fn to_kpc(&self) -> Self::Output {
        self * KM_TO_KPC
    }
}
