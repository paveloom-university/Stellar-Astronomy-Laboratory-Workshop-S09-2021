//! This module defines the [`Orbit`] struct
//!
//! Additional modules split the implementation of the struct
//! into different files.

use serde::Deserialize;

use crate::{F, I};

pub mod calc;
mod io;

/// Initial coordinates and velocities of a globular
/// cluster in the Heliocentric Cartesian system
#[derive(Deserialize)]
pub struct HCInitials {
    /// X component of the radius vector $\[ \text{kpc} \]$
    x: F,
    /// Standard deviation of the X component of the radius vector $\[ \text{kpc} \]$
    x_err: F,
    /// Y component of the radius vector $\[ \text{kpc} \]$
    y: F,
    /// Standard deviation of the Y component of the radius vector $\[ \text{kpc} \]$
    y_err: F,
    /// Z component of the radius vector $\[ \text{kpc} \]$
    z: F,
    /// Standard deviation of the Z component of the radius vector $\[ \text{kpc} \]$
    z_err: F,
    /// U component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    u: F,
    /// Standard deviation of the U component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    u_err: F,
    /// V component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    v: F,
    /// Standard deviation of the V component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    v_err: F,
    /// W component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    w: F,
    /// Standard deviation of the W component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    w_err: F,
}

/// Initial coordinates and velocities of a globular
/// cluster in the Galactic Cartesian system
pub struct GCaInitials {
    /// X component of the radius vector $\[ \text{kpc} \]$
    x: F,
    /// Y component of the radius vector $\[ \text{kpc} \]$
    y: F,
    /// Z component of the radius vector $\[ \text{kpc} \]$
    z: F,
    /// U component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    u: F,
    /// V component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    v: F,
    /// W component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    w: F,
}

impl GCaInitials {
    /// Initialize a new struct with default values
    fn new() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            u: 0.0,
            v: 0.0,
            w: 0.0,
        }
    }
}

/// Initial coordinates and velocities of a globular
/// cluster in the Galactic Cylindrical system
pub struct GCyInitials {
    /// Radius $\[ \text{kpc} \]$
    r: F,
    /// Azimuth $\[ \text{rad} \]$
    psi: F,
    /// Height $\[ \text{kpc} \]$
    z: F,
    /// Time derivative of R $\[ \text{km} \\, \text{s}^{-1} \]$
    dr: F,
    /// Angular velocity $\[ \text{rad} \\, \text{s}^{-1} \]$
    dpsi: F,
    /// Time derivative of Z $\[ \text{km} \\, \text{s}^{-1} \]$
    dz: F,
}

impl GCyInitials {
    /// Initialize a new struct with default values
    fn new() -> Self {
        Self {
            r: 0.0,
            psi: 0.0,
            z: 0.0,
            dr: 0.0,
            dpsi: 0.0,
            dz: 0.0,
        }
    }
}

/// Fields in the [`Results`] struct. This constant is used in the CLI
pub const RESULTS_FIELDS: &[&str] = &[
    "r", "psi", "z", "p_r", "p_psi", "p_z", "x", "y", "e", "apo", "peri",
];

/// Values integrated / calculated during runtime
pub struct Results {
    /// Radius in the Galactic Cylindrical system $\[ \text{kpc} \]$
    r: Vec<F>,
    /// Azimuth in the Galactic Cylindrical system $\[ \text{kpc} \]$
    psi: Vec<F>,
    /// Height in the Galactic Cartesian / Cylindrical system $\[ \text{kpc} \]$
    z: Vec<F>,
    /// Momentum canonically conjugate to the radius $\[ \text{kpc} \\, \text{Myr}^{-1} \]$
    p_r: Vec<F>,
    /// Momentum canonically conjugate to the azimuth $\[ \text{kpc}^2 \\, \text{rad} \\, \text{Myr}^{-1} \]$
    p_psi: Vec<F>,
    /// Momentum canonically conjugate to the height $\[ \text{kpc} \\, \text{Myr}^{-1} \]$
    p_z: Vec<F>,
    /// X component of the radius vector in the Galactic Cartesian system $\[ \text{kpc} \]$
    x: Vec<F>,
    /// Y component of the radius vector in the Galactic Cartesian system $\[ \text{kpc} \]$
    y: Vec<F>,
    /// Total energy $\[ \text{km}^2 \\, \text{s}^{-2} \]$
    e: Vec<F>,
    /// Apocentric distance $\[ \text{kpc} \]$
    apo: Vec<F>,
    /// Pericentric distance $\[ \text{kpc} \]$
    peri: Vec<F>,
}

impl Results {
    /// Initialize a new struct with default values
    fn new() -> Self {
        Self {
            r: Vec::<F>::new(),
            psi: Vec::<F>::new(),
            z: Vec::<F>::new(),
            p_r: Vec::<F>::new(),
            p_psi: Vec::<F>::new(),
            p_z: Vec::<F>::new(),
            x: Vec::<F>::new(),
            y: Vec::<F>::new(),
            e: Vec::<F>::new(),
            apo: Vec::<F>::new(),
            peri: Vec::<F>::new(),
        }
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
