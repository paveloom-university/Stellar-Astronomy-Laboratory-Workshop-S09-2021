//! Integrate the orbit

use super::super::{Orbit, F};

/// Galactocentric distance of the Local Standard
/// of Rest around the Galactic center (kpc)
const G_R_SUN: F = 8.3;

/// Galactocentric linear velocity of the Local Standard
/// of Rest around the Galactic center (km s^{-1})
const G_V_SUN: F = 244.0;

/// Height of the Sun above the Galactic plane (pc)
const G_H_SUN: F = 16.0;

/// The Sun’s peculiar velocity U with respect
/// to the Local Standard of Rest (km s^{-1})
const P_U_SUN: F = 11.1;

/// The Sun’s peculiar velocity V with respect
/// to the Local Standard of Rest (km s^{-1})
const P_V_SUN: F = 12.2;

/// The Sun’s peculiar velocity W with respect
/// to the Local Standard of Rest (km s^{-1})
const P_W_SUN: F = 7.3;

impl Orbit {
    /// Integrate the orbit
    pub fn integrate(&mut self) {
        // Convert the initial coordinates and velocities
        // from the Heliocentric Cartesian system to the
        // Galactic Cartesian system
        self.gca_initials.x = G_R_SUN - self.hc_initials.x_0;
        self.gca_initials.y = self.hc_initials.y_0;
        self.gca_initials.z = self.hc_initials.z_0 + G_H_SUN / 1000.0;
        self.gca_initials.u = self.hc_initials.u_0 + P_U_SUN;
        self.gca_initials.v = self.hc_initials.v_0 + P_V_SUN + G_V_SUN;
        self.gca_initials.w = self.hc_initials.w_0 + P_W_SUN;
        // Convert the initial coordinates and velocities
        // from the Galactic Cartesian system to the
        // Galactic Cylindrical system
        self.gcy_initials.r =
            (self.gca_initials.x.powf(2.0) + self.gca_initials.y.powf(2.0)).sqrt();
        self.gcy_initials.psi = F::atan(self.gca_initials.x / self.gca_initials.y);
        self.gcy_initials.z = self.gca_initials.z;
        self.gcy_initials.r_vel =
            (self.gca_initials.u.powf(2.0) + self.gca_initials.v.powf(2.0)).sqrt();
        self.gcy_initials.psi_vel = F::atan(self.gca_initials.u / self.gca_initials.v);
        self.gcy_initials.r_vel = self.gca_initials.w;
    }
}
