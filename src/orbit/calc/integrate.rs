//! Integrate the orbit

use super::potentials::{miyamoto_nagai, navarro_frenk_white, plummer};
use crate::orbit::{Orbit, F};

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
    /// Integrate the orbit with the specified
    /// number of iteration `n` and the time step `h`
    /// using the fourth-order Runge-Kutta algorithm
    pub fn integrate(&mut self, n: usize, h: F) {
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
        self.gcy_initials.r = (self.gca_initials.x.powi(2) + self.gca_initials.y.powi(2)).sqrt();
        self.gcy_initials.psi = F::atan(self.gca_initials.x / self.gca_initials.y);
        self.gcy_initials.z = self.gca_initials.z;
        self.gcy_initials.r_vel =
            (self.gca_initials.u.powi(2) + self.gca_initials.v.powi(2)).sqrt();
        self.gcy_initials.psi_vel = F::atan(self.gca_initials.u / self.gca_initials.v);
        self.gcy_initials.z_vel = self.gca_initials.w;

        // Prepare the vectors for the solution
        let mut r = Vec::<F>::with_capacity(n + 1);
        let mut psi = Vec::<F>::with_capacity(n + 1);
        let mut z = Vec::<F>::with_capacity(n + 1);
        let mut p_r = Vec::<F>::with_capacity(n + 1);
        let mut p_psi = Vec::<F>::with_capacity(n + 1);
        let mut p_z = Vec::<F>::with_capacity(n + 1);
        r.push(self.gcy_initials.r);
        psi.push(self.gcy_initials.psi);
        z.push(self.gcy_initials.z);
        p_r.push(self.gcy_initials.r_vel * 3.24078e-17);
        p_psi.push(self.gcy_initials.r.powi(2) * self.gcy_initials.psi_vel);
        p_z.push(self.gcy_initials.z_vel * 3.24078e-17);

        // Integrate the orbit using the fourth-order Runge-Kutta algorithm
        for i in 1..=n {
            r.push(r[i - 1] + runge_kutta_step(p_r[i - 1], h));
            psi.push(psi[i - 1] + runge_kutta_step(p_psi[i - 1] / r[i - 1].powi(2), h));
            z.push(z[i - 1] + runge_kutta_step(p_z[i - 1], h));
            p_r.push(
                p_r[i - 1]
                    + runge_kutta_step(
                        -phi_dr(r[i - 1], z[i - 1]) + p_psi[i - 1].powi(2) / r[i - 1].powi(3),
                        h,
                    ),
            );
            p_psi.push(p_psi[i - 1]);
            p_z.push(p_z[i - 1] + runge_kutta_step(-phi_dz(r[i - 1], z[i - 1]), h));
        }

        // Save the results
        self.gcy_integrated.r = r;
        self.gcy_integrated.psi = psi;
        self.gcy_integrated.z = z;
    }
}

/// Calculate the value step with the fourth-order Runge-Kutta method
fn runge_kutta_step(v: F, h: F) -> F {
    let k_1 = v;
    let k_2 = v + h * k_1 / 2.0;
    let k_3 = v + h * k_2 / 2.0;
    let k_4 = v + h * k_3;
    1.0 / 6.0 * h * (k_1 + 2.0 * k_2 + 2.0 * k_3 + k_4)
}

/// Calculate the value of the R derivative of the Galactic potential
fn phi_dr(r: F, z: F) -> F {
    plummer::phi_dr(r, z, 443.0, 0.2672)
        + miyamoto_nagai::phi_dr(r, z, 2798.0, 4.40, 0.3084)
        + navarro_frenk_white::phi_dr(r, z, 12474.0, 7.7)
}

/// Calculate the value of the Z derivative of the Galactic potential
fn phi_dz(r: F, z: F) -> F {
    plummer::phi_dz(r, z, 443.0, 0.2672)
        + miyamoto_nagai::phi_dz(r, z, 2798.0, 4.40, 0.3084)
        + navarro_frenk_white::phi_dz(r, z, 12474.0, 7.7)
}
