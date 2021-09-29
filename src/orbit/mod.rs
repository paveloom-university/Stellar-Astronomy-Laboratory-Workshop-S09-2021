//! This module defines the [`Orbit`] struct
//!
//! Additional modules split the implementation of the struct
//! into different files.

use serde::Deserialize;

use crate::{F, I};

mod calc;
mod io;

/// Initial coordinates and velocities of a globular
/// cluster in the Heliocentric Cartesian system
#[derive(Deserialize)]
pub struct HCInitials {
    /// X component of the radius vector (kpc)
    x: F,
    /// Y component of the radius vector (kpc)
    y: F,
    /// Z component of the radius vector (kpc)
    z: F,
    /// U component of the velocity vector (kpc)
    u: F,
    /// V component of the velocity vector (km/s)
    v: F,
    /// W component of the velocity vector (km/s)
    w: F,
}

/// Initial coordinates and velocities of a globular
/// cluster in the Galactic Cartesian system
pub struct GCaInitials {
    /// X component of the radius vector (kpc)
    x: F,
    /// Y component of the radius vector (kpc)
    y: F,
    /// Z component of the radius vector (kpc)
    z: F,
    /// U component of the velocity vector (km/s)
    u: F,
    /// V component of the velocity vector (km/s)
    v: F,
    /// W component of the velocity vector (km/s)
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
    /// Radius (kpc)
    r: F,
    /// Angle (rad)
    psi: F,
    /// Height (kpc)
    z: F,
    /// Time derivative of R (km/s)
    dr: F,
    /// Angular velocity (rad/s)
    dpsi: F,
    /// Time derivative of Z (km/s)
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

/// Some of the values integrated during runtime
///
/// These include:
/// - R and Z coordinates of a globular cluster in the Galactic Cylindrical system
/// - X and Y coordinates of a globular cluster in the Galactic Cartesian system
/// - Total energy
pub struct Integrated {
    /// R component of the radius vector in the Galactic Cylindrical system (kpc)
    r: Vec<F>,
    /// Z component of the radius vector in the Galactic Cylindrical system (kpc)
    z: Vec<F>,
    /// X component of the radius vector in the Galactic Cartesian system (kpc)
    x: Vec<F>,
    /// Y component of the radius vector in the Galactic Cartesian system (kpc)
    y: Vec<F>,
    /// Total energy (km^2/s^2)
    e: Vec<F>,
    vel: Vec<F>,
    phi: Vec<F>,
}

impl Integrated {
    /// Initialize a new struct with default values
    fn new() -> Self {
        Self {
            r: Vec::<F>::new(),
            z: Vec::<F>::new(),
            x: Vec::<F>::new(),
            y: Vec::<F>::new(),
            e: Vec::<F>::new(),
            vel: Vec::<F>::new(),
            phi: Vec::<F>::new(),
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
    integrated: Integrated,
}

impl Orbit {
    /// Initialize an orbit with an ID from the initial coordinates
    /// and velocities in the Heliocentric Cartesian system
    pub fn from(id: String, hc_initials: HCInitials) -> Self {
        Orbit {
            id,
            hc_initials,
            gca_initials: GCaInitials::new(),
            gcy_initials: GCyInitials::new(),
            n: 0,
            integrated: Integrated::new(),
        }
    }
    /// Get the ID of the object
    pub fn id(&self) -> &String {
        &self.id
    }
}
