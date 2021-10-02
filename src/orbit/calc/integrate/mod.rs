//! This module provides routines for orbit integration
//!
//! Here is a quick breakdown of what is going on (see
//! Bajkova & Bobylev ([2020, v1](https://arxiv.org/abs/2008.13624v1))
//! for more details).
//!
//! Let the initial positions and space velocities of a
//! test particle in the heliocentric coordinate system be
//! $ (x_0, y_0 , z_0, u_0, v_0 , w_0) $. The initial positions
//! $ (X, Y, Z) $ and velocities $ (U, V, W) $ of the test
//! particle in Galactic Cartesian coordinates are then given
//! by the formulas:
//!
//! $$
//! X \\, [\text{kpc}] = R_\odot - x_0; \\\\
//! Y \\, [\text{kpc}] = y_0; \\\\
//! Z \\, [\text{kpc}] = z_0 + h_\odot; \\\\
//! U \\, [\text{km} \\, \text{s}^{-1}] = u_0 + u_\odot; \\\\
//! V \\, [\text{km} \\, \text{s}^{-1}] = v_0 + v_\odot + V_\odot; \\\\
//! W \\, [\text{km} \\, \text{s}^{-1}] = w_0 + w_\odot,
//! $$
//!
//! where $ R_\odot = 8.3 \\, \text{kpc} $ and $ V_\odot = 244 \\,
//! \text{km} \\, \text{s}^{-1} $ are the Galactocentric distance
//! and the linear velocity of the Local Standard of Rest around
//! the Galactic center, $ h_\odot = 16 \\, \text{pc} $ (Bobylev &
//! Bajkova (2016)) is the height of the Sun above the Galactic plane.
//!
//! The initial positions $ (R, \psi, Z) $ and their time derivatives
//! $ (\dot{R}, \dot{\psi}, \dot{Z}) $ are then obtained by the
//! formulas:
//!
//! $$
//! R \\, [\text{kpc}] = \sqrt{X^2 + Y^2}; \\\\
//! \psi \\, [\text{rad}] = \text{atan2}(Y, X); \\\\
//! \dot{R} \\, [\text{km} \\, \text{s}^{-1}] = -U \cos{\psi} + V \sin{\psi}; \\\\
//! \dot{\psi} \\, [\text{rad} \\, \text{s}^{-1}] = (U \sin{\psi} + V \cos{\psi})
//! / R \\, [\text{kpc} \rightarrow \text{km}]; \\\\
//! \dot{Z} \\, [\text{km} \\, \text{s}^{-1}] = W.
//! $$
//!
//! From here we obtain the initial values of the canonical moments
//! $ (p_R, p_\psi, p_Z) $:
//!
//! $$
//! p_R \\, [\text{kpc} \\, \text{Myr}^{-1}] = \dot{R} \\, [\text{km} \\, \text{s}^{-1}
//! \rightarrow \text{kpc} \\, \text{Myr}^{-1}]; \\\\
//! p_{\psi} \\, [\text{kpc}^2 \\, \text{rad} \\, \text{Myr}^{-1}] = R^2 \dot{\psi} \\, [\text{s}^{-1}
//! \rightarrow \text{Myr}^{-1}]; \\\\
//! p_Z \\, [\text{kpc} \\, \text{Myr}^{-1}] = \dot{Z} \\, [\text{km} \\, \text{s}^{-1}
//! \rightarrow \text{kpc} \\, \text{Myr}^{-1}].
//! $$
//!
//! Now we're all set to integrate the Lagrangian equations
//!
//! $$
//! \dot{R} \\, [\text{kpc} \\, \text{Myr}^{-1}] = p_R; \\\\
//! \dot{\psi} \\, [\text{rad} \\, \text{Myr}^{-1}] = p_\psi / R^2; \\\\
//! \dot{Z} \\, [\text{kpc} \\, \text{Myr}^{-1}] = p_Z; \\\\
//! \dot{p_R} \\, [\text{kpc} \\, \text{Myr}^{-2}] = - \partial \Phi(R, Z) / \partial R
//! \\, [100 \\, \text{km} \\, \text{s}^{-2} \rightarrow \text{kpc} \\, \text{Myr}^{-2}]
//! + p_\psi^2 / R^3; \\\\
//! \dot{p_\psi} \\, [\text{kpc}^2 \\, \text{Myr}^{-2}] = 0; \\\\
//! \dot{p_Z} \\, [\text{kpc} \\, \text{Myr}^{-2}] = - \partial \Phi(R, Z) / \partial Z
//! \\, [100 \\, \text{km} \\, \text{s}^{-2} \rightarrow \text{kpc} \\, \text{Myr}^{-2}].
//! $$
//!
//! To calculate the total energy, we first obtain radial velocity
//!
//! $$
//! \Pi \\, [\text{km} \\, \text{s}^{-1}] = -U \frac{X}{R} + V \frac{Y}{R} =
//! -(p_r \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\, \text{s}^{-1}]
//! \cos{\psi} - (p_\psi / r) \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\,
//! \text{s}^{-1}] \sin{\psi}) \cos{\psi} \\\\
//! + (p_r \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\, \text{s}^{-1}]
//! \sin{\psi} - (p_\psi / r) \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\,
//! \text{s}^{-1}] \cos{\psi}) \sin{\psi},
//! $$
//!
//! tangential velocity
//!
//! $$
//! \Theta \\, [\text{km} \\, \text{s}^{-1}] = U \frac{Y}{R} + V \frac{X}{R} =
//! (p_r \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\, \text{s}^{-1}]
//! \cos{\psi} - (p_\psi / r) \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\,
//! \text{s}^{-1}] \sin{\psi}) \sin{\psi} \\\\
//! + (p_r \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\, \text{s}^{-1}]
//! \sin{\psi} - (p_\psi / r) \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\,
//! \text{s}^{-1}] \cos{\psi}) \cos{\psi},
//! $$
//!
//! and total 3D velocity $ V_\text{tot} \\, [\text{km}^2 \\, \text{s}^{-2}] =
//! \sqrt{\Pi^2 + \Theta^2 + W^2} = \sqrt{\Pi^2 + \Theta^2 + (p_z \\, [\text{kpc} \\,
//! \text{Myr}^{-1} \rightarrow \text{km} \\, \text{s}^{-1}])^2} $.
//!
//! Then we can calculate the total energy $ E \\, [\text{km}^2 \\, \text{s}^{-2}] =
//! \Phi(R, Z) \\, [\text{km}^2 \\, \text{s}^{-2} \rightarrow \text{100} \\,
//! \text{km}^2 \\, \text{s}^{-2}] + V_\text{tot}^2 / 2 $.

use super::potentials::{miyamoto_nagai, navarro_frenk_white, plummer};
use crate::{orbit::Orbit, F, I};

mod constants;

use constants::{conversion::Conversions, correction::Corrections};

/// Calculate the value of the Galactic potential $ \Phi(R, Z) $
/// $\[ 100 \\, \text{km}^2 \\, \text{s}^{-2} \]$
fn phi(r: F, z: F) -> F {
    plummer::phi(r, z) + miyamoto_nagai::phi(r, z) + navarro_frenk_white::phi(r, z)
}

/// Calculate the value of $ \partial \Phi(R, Z) / \partial R $
/// $\[ \text{kpc} \\, \text{Myr}^{-2} \]$
fn phi_dr(r: F, z: F) -> F {
    (plummer::phi_dr(r, z) + miyamoto_nagai::phi_dr(r, z) + navarro_frenk_white::phi_dr(r, z))
        .to_kpc_per_myr_2()
}

/// Calculate the value of $ \partial \Phi(R, Z) / \partial Z $
/// $\[ \text{kpc} \\, \text{Myr}^{-2} \]$
fn phi_dz(r: F, z: F) -> F {
    (plummer::phi_dz(r, z) + miyamoto_nagai::phi_dz(r, z) + navarro_frenk_white::phi_dz(r, z))
        .to_kpc_per_myr_2()
}

/// Calculate the right-hand part of the equation for
/// $ \dot{R} $ $\[ \text{kpc} \\, \text{Myr}^{-1} \]$
fn f_r(p_r: F) -> F {
    p_r
}

/// Calculate the right-hand part of the equation for
/// $ \dot{\psi} $ $\[ \text{Myr}^{-1} \]$
fn f_psi(r: F, p_psi: F) -> F {
    p_psi / r.powi(2)
}

/// Calculate the right-hand part of the equation for
/// $ \dot{z} $ $\[ \text{kpc} \\, \text{Myr}^{-1} \]$
fn f_z(p_z: F) -> F {
    p_z
}

/// Calculate the right-hand part of the equation for
/// $ \dot{p_r} $ $\[ \text{kpc} \\, \text{Myr}^{-2} \]$
fn f_p_r(r: F, z: F, p_psi: F) -> F {
    -phi_dr(r, z) + p_psi.powi(2) / r.powi(3)
}

/// Calculate the right-hand part of the equation for
/// $ \dot{p_z} $ $\[ \text{kpc} \\, \text{Myr}^{-2} \]$
fn f_p_z(r: F, z: F) -> F {
    -phi_dz(r, z)
}

/// Calculate the value of radial velocity
/// $ \Pi $ $\[ \text{km} \\, \text{s}^{-1} \]$
fn radial_velocity(r: F, psi: F, p_r: F, p_psi: F) -> F {
    -(p_r.to_km_per_s() * psi.cos() - (p_psi / r).to_km_per_s() * psi.sin()) * psi.cos()
        + (p_r.to_km_per_s() * psi.sin() + (p_psi / r).to_km_per_s() * psi.cos()) * psi.sin()
}

/// Calculate the value of tangential velocity
/// $ \Theta $ $\[ \text{km} \\, \text{s}^{-1} \]$
fn tangential_velocity(r: F, psi: F, p_r: F, p_psi: F) -> F {
    (p_r.to_km_per_s() * psi.cos() - (p_psi / r).to_km_per_s() * psi.sin()) * psi.sin()
        + (p_r.to_km_per_s() * psi.sin() + (p_psi / r).to_km_per_s() * psi.cos()) * psi.cos()
}

/// Calculate the value of total energy E $\[ \text{km}^2 \\, \text{s}^{-2} \]$
fn total_energy(r: F, psi: F, z: F, p_r: F, p_psi: F, p_z: F) -> F {
    (radial_velocity(r, psi, p_r, p_psi).powi(2)
        + tangential_velocity(r, psi, p_r, p_psi).powi(2)
        + p_z.to_km_per_s().powi(2))
        / 2.0
        + phi(r, z) * 100.0
}

/// Integrate the orbit using the fourth-order Runge-Kutta algorithm
fn rk4(o: &mut Orbit, h: F, s_idx: I, f_idx: I) {
    // Prepare buffers for the increments
    let mut k_r_1: F;
    let mut k_r_2: F;
    let mut k_r_3: F;
    let mut k_r_4: F;
    let mut k_psi_1: F;
    let mut k_psi_2: F;
    let mut k_psi_3: F;
    let mut k_psi_4: F;
    let mut k_z_1: F;
    let mut k_z_2: F;
    let mut k_z_3: F;
    let mut k_z_4: F;
    let mut k_p_r_1: F;
    let mut k_p_r_2: F;
    let mut k_p_r_3: F;
    let mut k_p_r_4: F;
    let mut k_p_z_1: F;
    let mut k_p_z_2: F;
    let mut k_p_z_3: F;
    let mut k_p_z_4: F;

    // Indices range depends on the forward or reverse mode
    for i in s_idx..=f_idx {
        // Calculate the increments for all integrable variables
        // (except for `p_psi`, since its derivative is 0)

        k_r_1 = h * f_r(o.results.p_r[i - 1]);
        k_psi_1 = h * f_psi(o.results.r[i - 1], o.results.p_psi[i - 1]);
        k_z_1 = h * f_z(o.results.p_z[i - 1]);
        k_p_r_1 = h * f_p_r(
            o.results.r[i - 1],
            o.results.z[i - 1],
            o.results.p_psi[i - 1],
        );
        k_p_z_1 = h * f_p_z(o.results.r[i - 1], o.results.z[i - 1]);

        k_r_2 = h * f_r(o.results.p_r[i - 1] + 0.5 * k_p_r_1);
        k_psi_2 = h * f_psi(o.results.r[i - 1] + 0.5 * k_r_1, o.results.p_psi[i - 1]);
        k_z_2 = h * f_z(o.results.p_z[i - 1] + 0.5 * k_p_z_1);
        k_p_r_2 = h * f_p_r(
            o.results.r[i - 1] + 0.5 * k_r_1,
            o.results.z[i - 1] + 0.5 * k_z_1,
            o.results.p_psi[i - 1],
        );
        k_p_z_2 = h * f_p_z(
            o.results.r[i - 1] + 0.5 * k_r_1,
            o.results.z[i - 1] + 0.5 * k_z_1,
        );

        k_r_3 = h * f_r(o.results.p_r[i - 1] + 0.5 * k_p_r_2);
        k_psi_3 = h * f_psi(o.results.r[i - 1] + 0.5 * k_r_2, o.results.p_psi[i - 1]);
        k_z_3 = h * f_z(o.results.p_z[i - 1] + 0.5 * k_p_z_2);
        k_p_r_3 = h * f_p_r(
            o.results.r[i - 1] + 0.5 * k_r_2,
            o.results.z[i - 1] + 0.5 * k_z_2,
            o.results.p_psi[i - 1],
        );
        k_p_z_3 = h * f_p_z(
            o.results.r[i - 1] + 0.5 * k_r_2,
            o.results.z[i - 1] + 0.5 * k_z_2,
        );

        k_r_4 = h * f_r(o.results.p_r[i - 1] + k_p_r_3);
        k_psi_4 = h * f_psi(o.results.r[i - 1] + k_r_3, o.results.p_psi[i - 1]);
        k_z_4 = h * f_z(o.results.p_z[i - 1] + k_p_z_3);
        k_p_r_4 = h * f_p_r(
            o.results.r[i - 1] + k_r_3,
            o.results.z[i - 1] + k_z_3,
            o.results.p_psi[i - 1],
        );
        k_p_z_4 = h * f_p_z(o.results.r[i - 1] + k_r_3, o.results.z[i - 1] + k_z_3);

        // Push the next values

        // Galactic Cylindrical system
        o.results
            .r
            .push(o.results.r[i - 1] + 1.0 / 6.0 * (k_r_1 + 2.0 * k_r_2 + 2.0 * k_r_3 + k_r_4));
        o.results.psi.push(
            o.results.psi[i - 1] + 1.0 / 6.0 * (k_psi_1 + 2.0 * k_psi_2 + 2.0 * k_psi_3 + k_psi_4),
        );
        o.results
            .z
            .push(o.results.z[i - 1] + 1.0 / 6.0 * (k_z_1 + 2.0 * k_z_2 + 2.0 * k_z_3 + k_z_4));
        o.results.p_r.push(
            o.results.p_r[i - 1] + 1.0 / 6.0 * (k_p_r_1 + 2.0 * k_p_r_2 + 2.0 * k_p_r_3 + k_p_r_4),
        );
        o.results.p_psi.push(o.results.p_psi[i - 1]);
        o.results.p_z.push(
            o.results.p_z[i - 1] + 1.0 / 6.0 * (k_p_z_1 + 2.0 * k_p_z_2 + 2.0 * k_p_z_3 + k_p_z_4),
        );

        // Galactic Cartesian system
        o.results.x.push(o.results.r[i] * o.results.psi[i].cos());
        o.results.y.push(o.results.r[i] * o.results.psi[i].sin());

        // Total energy
        o.results.e.push(total_energy(
            o.results.r[i],
            o.results.psi[i],
            o.results.z[i],
            o.results.p_r[i],
            o.results.p_psi[i],
            o.results.p_z[i],
        ));
    }
}

impl Orbit {
    /// Integrate the orbit with the specified
    /// number of iterations `n` and the time step `h`
    /// using the fourth-order Runge-Kutta algorithm
    pub fn integrate(&mut self, n: I, h: F, rev: bool) {
        // Save the number of integrations
        self.n = if rev { n * 2 } else { n };

        // Convert the initial coordinates and velocities
        // from the Heliocentric Cartesian system to the
        // Galactic Cartesian system
        self.gca_initials.x = self.hc_initials.x.to_gca_x();
        self.gca_initials.y = self.hc_initials.y;
        self.gca_initials.z = self.hc_initials.z.to_gca_z();
        self.gca_initials.u = self.hc_initials.u.to_gca_u();
        self.gca_initials.v = self.hc_initials.v.to_gca_v();
        self.gca_initials.w = self.hc_initials.w.to_gca_w();

        // Convert the initial coordinates and velocities
        // from the Galactic Cartesian system to the
        // Galactic Cylindrical system
        self.gcy_initials.r = (self.gca_initials.x.powi(2) + self.gca_initials.y.powi(2)).sqrt();
        self.gcy_initials.psi = self.gca_initials.y.atan2(self.gca_initials.x);
        self.gcy_initials.z = self.gca_initials.z;
        self.gcy_initials.dr = -self.gca_initials.u * self.gcy_initials.psi.cos()
            + self.gca_initials.v * self.gcy_initials.psi.sin();
        self.gcy_initials.dpsi = ((self.gca_initials.u * self.gcy_initials.psi.sin()
            + self.gca_initials.v * self.gcy_initials.psi.cos())
            / self.gcy_initials.r)
            .to_kpc();
        self.gcy_initials.dz = self.gca_initials.w;

        // Prepare the vectors for the solutions

        // Lagrangian equations
        self.results.r = Vec::<F>::with_capacity(self.n + 1);
        self.results.psi = Vec::<F>::with_capacity(self.n + 1);
        self.results.z = Vec::<F>::with_capacity(self.n + 1);
        self.results.p_r = Vec::<F>::with_capacity(self.n + 1);
        self.results.p_psi = Vec::<F>::with_capacity(self.n + 1);
        self.results.p_z = Vec::<F>::with_capacity(self.n + 1);

        // Galactic Cartesian system
        self.results.x = Vec::<F>::with_capacity(self.n + 1);
        self.results.y = Vec::<F>::with_capacity(self.n + 1);

        // Total energy
        self.results.e = Vec::<F>::with_capacity(self.n + 1);

        // Push the initial values

        // Lagrangian equations
        self.results.r.push(self.gcy_initials.r);
        self.results.psi.push(self.gcy_initials.psi);
        self.results.z.push(self.gcy_initials.z);
        self.results.p_r.push(self.gcy_initials.dr.to_kpc_per_myr());
        self.results
            .p_psi
            .push((self.gcy_initials.r.powi(2) * self.gcy_initials.dpsi).to_seconds());
        self.results.p_z.push(self.gcy_initials.dz.to_kpc_per_myr());

        // Galactic Cartesian system
        self.results
            .x
            .push(self.results.r[0] * self.results.psi[0].cos());
        self.results
            .y
            .push(self.results.r[0] * self.results.psi[0].sin());

        // Total energy
        self.results.e.push(total_energy(
            self.results.r[0],
            self.results.psi[0],
            self.results.z[0],
            self.results.p_r[0],
            self.results.p_psi[0],
            self.results.p_z[0],
        ));

        // Integrate the orbit using the fourth-order Runge-Kutta algorithm
        rk4(self, h, 1, n);
        // Integrate in reverse if specified by the user
        if rev {
            rk4(self, -h, n + 1, self.n + 1);
        }
    }
}
