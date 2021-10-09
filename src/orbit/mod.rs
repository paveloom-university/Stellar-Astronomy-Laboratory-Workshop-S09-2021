//! This module defines the [`Orbit`] struct
//!
//! Additional modules split the implementation of the struct
//! into different files.

use serde::Deserialize;

use crate::{F, I};

pub mod calc;
mod io;

/// Create a struct with a common fields type
macro_rules! def_struct {
    (
        $(#[$attr:meta])*
        $vis:vis struct $name:ident($type:ty) {
            $(
                #[$field_doc:meta]
                $field:ident,
            )*
        }
    ) => {
        $(#[$attr])*
        $vis struct $name {
            $(
                #[$field_doc]
                $field: $type,
            )*
        }
    }
}

/// Create a struct with a common fields type and a default initializer
macro_rules! init_struct {
    (
        $(#[$attr:meta])*
        $vis:vis struct $name:ident($type:ty, $expr:expr) {
            $(
                #[$field_doc:meta]
                $field:ident,
            )*
        }
    ) => {
        $(#[$attr])*
        $vis struct $name {
            $(
                #[$field_doc]
                $field: $type,
            )*
        }

        impl $name {
            /// Initialize a new struct with default values
            fn new() -> Self {
                Self {
                    $(
                        $field: $expr,
                    )*
                }
            }
        }
    }
}

/// Create a struct with a common fields type, a default initializer,
/// and a method to get fields as static string slices
macro_rules! fields_struct {
    (
        $(#[$attr:meta])*
        $vis:vis struct $name:ident($type:ty, $expr:expr) {
            $(
                #[$field_doc:meta]
                $field:ident,
            )*
        }
    ) => {
        $(#[$attr])*
        $vis struct $name {
            $(
                #[$field_doc]
                $field: $type,
            )*
        }

        impl $name {
            /// Initialize a new struct with default values
            fn new() -> Self {
                Self {
                    $(
                        $field: $expr,
                    )*
                }
            }
            /// Get the fields of the struct as static string slices
            #[must_use]
            pub fn fields() -> &'static [&'static str] {
                static NAMES: &'static [&'static str] = &[$(stringify!($field)),*];
                NAMES
            }
        }
    }
}

def_struct! {
/// Initial coordinates and velocities of a globular
/// cluster in the Heliocentric Cartesian system
#[derive(Deserialize)]
pub struct HCInitials(F) {
    /// X component of the radius vector $\[ \text{kpc} \]$
    x,
    /// Standard deviation of the X component of the radius vector $\[ \text{kpc} \]$
    x_err,
    /// Y component of the radius vector $\[ \text{kpc} \]$
    y,
    /// Standard deviation of the Y component of the radius vector $\[ \text{kpc} \]$
    y_err,
    /// Z component of the radius vector $\[ \text{kpc} \]$
    z,
    /// Standard deviation of the Z component of the radius vector $\[ \text{kpc} \]$
    z_err,
    /// U component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    u,
    /// Standard deviation of the U component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    u_err,
    /// V component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    v,
    /// Standard deviation of the V component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    v_err,
    /// W component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    w,
    /// Standard deviation of the W component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    w_err,
}
}

init_struct! {
/// Initial coordinates and velocities of a globular
/// cluster in the Galactic Cartesian system
pub struct GCaInitials(F, 0.0) {
    /// X component of the radius vector $\[ \text{kpc} \]$
    x,
    /// Y component of the radius vector $\[ \text{kpc} \]$
    y,
    /// Z component of the radius vector $\[ \text{kpc} \]$
    z,
    /// U component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    u,
    /// V component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    v,
    /// W component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    w,
}
}

init_struct! {
/// Initial coordinates and velocities of a globular
/// cluster in the Galactic Cylindrical system
pub struct GCyInitials(F, 0.0) {
    /// Radius $\[ \text{kpc} \]$
    r,
    /// Azimuth $\[ \text{rad} \]$
    psi,
    /// Height $\[ \text{kpc} \]$
    z,
    /// Time derivative of R $\[ \text{km} \\, \text{s}^{-1} \]$
    dr,
    /// Angular velocity $\[ \text{rad} \\, \text{s}^{-1} \]$
    dpsi,
    /// Time derivative of Z $\[ \text{km} \\, \text{s}^{-1} \]$
    dz,
}
}

fields_struct! {
/// Values integrated / calculated during runtime
pub struct Results(Vec<F>, Vec::<F>::new()) {
    /// Radius in the Galactic Cylindrical system $\[ \text{kpc} \]$
    r,
    /// Azimuth in the Galactic Cylindrical system $\[ \text{kpc} \]$
    psi,
    /// Height in the Galactic Cartesian / Cylindrical system $\[ \text{kpc} \]$
    z,
    /// Momentum canonically conjugate to the radius $\[ \text{kpc} \\, \text{Myr}^{-1} \]$
    p_r,
    /// Momentum canonically conjugate to the azimuth $\[ \text{kpc}^2 \\, \text{rad} \\, \text{Myr}^{-1} \]$
    p_psi,
    /// Momentum canonically conjugate to the height $\[ \text{kpc} \\, \text{Myr}^{-1} \]$
    p_z,
    /// X component of the radius vector in the Galactic Cartesian system $\[ \text{kpc} \]$
    x,
    /// Y component of the radius vector in the Galactic Cartesian system $\[ \text{kpc} \]$
    y,
    /// Total energy $\[ \text{km}^2 \\, \text{s}^{-2} \]$
    e,
    /// Apocentric distance $\[ \text{kpc} \]$
    apo,
    /// Pericentric distance $\[ \text{kpc} \]$
    peri,
}
}

/// An orbit of a globular cluster
pub struct Orbit {
    /// ID of the object
    id: String,
    /// Initial coordinates and velocities in the Heliocentric Cartesian system
    hc_initials: HCInitials,
    /// Initial coordinates and velocities in the Galactic Cartesian system
    gca_initials: GCaInitials,
    /// Initial coordinates and velocities in the Galactic Cylindrical system
    gcy_initials: GCyInitials,
    /// Number of iterations used in the last integration
    n: I,
    /// Integrated coordinates of a globular cluster in the Galactic Cylindrical system
    results: Results,
}

impl Orbit {
    /// Initialize an orbit with an ID from the initial coordinates
    /// and velocities in the Heliocentric Cartesian system
    #[must_use]
    pub fn from(id: String, hc_initials: HCInitials) -> Self {
        Orbit {
            id,
            hc_initials,
            gca_initials: GCaInitials::new(),
            gcy_initials: GCyInitials::new(),
            n: 0,
            results: Results::new(),
        }
    }
    /// Get the ID of the object
    #[must_use]
    pub fn id(&self) -> &String {
        &self.id
    }
}
