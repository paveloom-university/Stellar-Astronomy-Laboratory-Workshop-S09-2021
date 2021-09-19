//! This module defines the [`Orbit`] struct
//!
//! Additional modules split the implementation of the struct
//! into different files.

use serde::Deserialize;

mod io;

/// The floating point type used for calculations
type F = f64;

/// A representation of the initial coordinates and velocities of a
/// globular cluster in the heliocentric system
#[derive(Deserialize)]
pub struct Initials {
    /// Initial X coordinate in the heliocentric system
    x_0: F,
    /// Error of the initial X coordinate in the heliocentric system
    x_0_err: F,
    /// Initial Y coordinate in the heliocentric system
    y_0: F,
    /// Error of the initial Y coordinate in the heliocentric system
    y_0_err: F,
    /// Initial Z coordinate in the heliocentric system
    z_0: F,
    /// Error of the initial Z coordinate in the heliocentric system
    z_0_err: F,
    /// Initial U velocity in the heliocentric system
    u_0: F,
    /// Error of the initial U velocity in the heliocentric system
    u_0_err: F,
    /// Initial V velocity in the heliocentric system
    v_0: F,
    /// Error of the initial V velocity in the heliocentric system
    v_0_err: F,
    /// Initial W velocity in the heliocentric system
    w_0: F,
    /// Error of the initial W velocity in the heliocentric system
    w_0_err: F,
}

/// A representation of an orbit of a globular cluster
pub struct Orbit {
    /// Initial coordinates and velocities in the heliocentric system
    initials: Initials,
}

impl Orbit {
    /// Initialize an orbit from the initial coordinates and velocities
    /// in the heliocentric system
    pub fn from(initials: Initials) -> Self {
        Orbit { initials }
    }
}
