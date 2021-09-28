//! This module defines the [`Orbit`] struct
//!
//! Additional modules split the implementation of the struct
//! into different files.

use serde::Deserialize;

mod calc;
mod io;

/// The floating point type used for calculations
type F = f64;

/// Initial coordinates and velocities of a globular
/// cluster in the Heliocentric Cartesian system
#[derive(Deserialize)]
pub struct HeliocentricCartesianInitials {
    /// X coordinate (kpc)
    x: F,
    /// Error of the X coordinate (kpc)
    x_err: F,
    /// Y coordinate
    y: F,
    /// Error of the Y coordinate (kpc)
    y_err: F,
    /// Z coordinate (kpc)
    z: F,
    /// Error of the Z coordinate (kpc)
    z_err: F,
    /// U velocity (kpc)
    u: F,
    /// Error of the U velocity (km s^{-1})
    u_err: F,
    /// V velocity (km s^{-1})
    v: F,
    /// Error of the V velocity (km s^{-1})
    v_err: F,
    /// W velocity (km s^{-1})
    w: F,
    /// Error of the W velocity (km s^{-1})
    w_err: F,
}

/// Initial coordinates and velocities of a globular
/// cluster in the Galactic Cartesian system
pub struct GalacticCartesianInitials {
    /// X coordinate (kpc)
    x: F,
    /// Error of the X coordinate (kpc)
    x_err: F,
    /// Y coordinate (kpc)
    y: F,
    /// Error of the Y coordinate (kpc)
    y_err: F,
    /// Z coordinate (kpc)
    z: F,
    /// Error of the Z coordinate (kpc)
    z_err: F,
    /// U velocity (km s^{-1})
    u: F,
    /// Error of the U velocity (km s^{-1})
    u_err: F,
    /// V velocity (km s^{-1})
    v: F,
    /// Error of the V velocity (km s^{-1})
    v_err: F,
    /// W velocity (km s^{-1})
    w: F,
    /// Error of the W velocity (km s^{-1})
    w_err: F,
}

impl GalacticCartesianInitials {
    /// Initialize a new struct with default values
    fn new() -> Self {
        GalacticCartesianInitials {
            x: 0.0,
            x_err: 0.0,
            y: 0.0,
            y_err: 0.0,
            z: 0.0,
            z_err: 0.0,
            u: 0.0,
            u_err: 0.0,
            v: 0.0,
            v_err: 0.0,
            w: 0.0,
            w_err: 0.0,
        }
    }
}

/// Initial coordinates and velocities of a globular
/// cluster in the Galactic Cylindrical system
pub struct GalacticCylindricalInitials {
    /// R coordinate (kpc)
    r: F,
    /// Error of the R coordinate (kpc)
    r_err: F,
    /// $ \psi $ coordinate (radians)
    psi: F,
    /// Error of the $ \psi $ coordinate (radians)
    psi_err: F,
    /// Z coordinate (kpc)
    z: F,
    /// Error of the Z coordinate (kpc)
    z_err: F,
    /// R velocity (km s^{-1})
    r_vel: F,
    /// Error of the R velocity (km s^{-1})
    r_vel_err: F,
    /// $ \psi $ velocity (radians)
    psi_vel: F,
    /// Error of the $ \psi $ velocity (radians)
    psi_vel_err: F,
    /// Z velocity (km s^{-1})
    z_vel: F,
    /// Error of the Z velocity (km s^{-1})
    z_vel_err: F,
}

impl GalacticCylindricalInitials {
    /// Initialize a new struct with default values
    fn new() -> Self {
        Self {
            r: 0.0,
            r_err: 0.0,
            psi: 0.0,
            psi_err: 0.0,
            z: 0.0,
            z_err: 0.0,
            r_vel: 0.0,
            r_vel_err: 0.0,
            psi_vel: 0.0,
            psi_vel_err: 0.0,
            z_vel: 0.0,
            z_vel_err: 0.0,
        }
    }
}

/// Values integrated during runtime:
/// - Coordinates of a globular cluster in the Galactic Cylindrical system
/// - Total energy
pub struct Integrated {
    /// R coordinate in the Galactic Cylindrical system (kpc)
    r: Vec<F>,
    /// $ \psi $ coordinate in the Galactic Cylindrical system (radians)
    psi: Vec<F>,
    /// Z coordinate in the Galactic Cylindrical system (kpc)
    z: Vec<F>,
    /// X coordinate in the Galactic Cartesian system (kpc)
    x: Vec<F>,
    /// Y coordinate in the Galactic Cartesian system (kpc)
    y: Vec<F>,
    /// Total energy (km^2 / s^2)
    e: Vec<F>,
    vel: Vec<F>,
    phi: Vec<F>,
}

impl Integrated {
    /// Initialize a new struct with default values
    fn new() -> Self {
        Self {
            r: Vec::<F>::new(),
            psi: Vec::<F>::new(),
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
    hc_initials: HeliocentricCartesianInitials,
    /// Initial coordinates and velocities in the Galactic Cartesian system
    gca_initials: GalacticCartesianInitials,
    /// Initial coordinates and velocities in the Galactic Cylindrical system
    gcy_initials: GalacticCylindricalInitials,
    /// Number of iterations used in the last integration
    n: usize,
    /// Integrated coordinates of a globular cluster in the Galactic Cylindrical system
    integrated: Integrated,
}

impl Orbit {
    /// Initialize an orbit with an ID from the initial coordinates
    /// and velocities in the Heliocentric Cartesian system
    pub fn from(id: String, hc_initials: HeliocentricCartesianInitials) -> Self {
        Orbit {
            id,
            hc_initials,
            gca_initials: GalacticCartesianInitials::new(),
            gcy_initials: GalacticCylindricalInitials::new(),
            n: 0,
            integrated: Integrated::new(),
        }
    }
    /// Get the ID of the object
    pub fn id(&self) -> &String {
        &self.id
    }
}
