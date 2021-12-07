//! This module provides routines for the orbit integration
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
//! -(-p_r \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\, \text{s}^{-1}]
//! \cos{\psi} + (p_\psi / r) \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\,
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
//! (-p_r \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\, \text{s}^{-1}]
//! \cos{\psi} + (p_\psi / r) \\, [\text{kpc} \\, \text{Myr}^{-1} \rightarrow \text{km} \\,
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

use crate::orbit::calc::models::Model;
use crate::{orbit::Orbit, F, I};

mod constants;

use constants::{conversion::Conversions, correction::Corrections};

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
fn f_p_r(r: F, z: F, p_psi: F, model: &dyn Model) -> F {
    -model.phi_dr(r, z).to_kpc_per_myr_2() + p_psi.powi(2) / r.powi(3)
}

/// Calculate the right-hand part of the equation for
/// $ \dot{p_z} $ $\[ \text{kpc} \\, \text{Myr}^{-2} \]$
fn f_p_z(r: F, z: F, model: &dyn Model) -> F {
    -model.phi_dz(r, z).to_kpc_per_myr_2()
}

/// Calculate the value of radial velocity
/// $ \Pi $ $\[ \text{km} \\, \text{s}^{-1} \]$
fn radial_velocity(r: F, psi: F, p_r: F, p_psi: F) -> F {
    -(-p_r.to_km_per_s() * psi.cos() + (p_psi / r).to_km_per_s() * psi.sin()) * psi.cos()
        + (p_r.to_km_per_s() * psi.sin() + (p_psi / r).to_km_per_s() * psi.cos()) * psi.sin()
}

/// Calculate the value of tangential velocity
/// $ \Theta $ $\[ \text{km} \\, \text{s}^{-1} \]$
fn tangential_velocity(r: F, psi: F, p_r: F, p_psi: F) -> F {
    (-p_r.to_km_per_s() * psi.cos() + (p_psi / r).to_km_per_s() * psi.sin()) * psi.sin()
        + (p_r.to_km_per_s() * psi.sin() + (p_psi / r).to_km_per_s() * psi.cos()) * psi.cos()
}

/// Calculate the value of total energy E $\[ \text{km}^2 \\, \text{s}^{-2} \]$
fn total_energy(r: F, psi: F, z: F, p_r: F, p_psi: F, p_z: F, model: &dyn Model) -> F {
    (radial_velocity(r, psi, p_r, p_psi).powi(2)
        + tangential_velocity(r, psi, p_r, p_psi).powi(2)
        + p_z.to_km_per_s().powi(2))
        / 2.0
        + model.phi(r, z) * 100.0
}

/// Integrate the orbit using the fourth-order Runge-Kutta algorithm
fn rk4(orbit: &mut Orbit, model: &dyn Model, h: F, s_idx: I, f_idx: I) {
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

    // Prepare aliases
    let (r, psi, z, p_r, p_psi, p_z, x, y, e, _, _) = orbit.results.unpack();

    // Indices range depends on the forward or reverse mode
    for i in s_idx..=f_idx {
        // Calculate the increments for all integrable variables
        // (except for `p_psi`, since its derivative is 0)

        k_r_1 = h * f_r(p_r[i - 1]);
        k_psi_1 = h * f_psi(r[i - 1], p_psi[i - 1]);
        k_z_1 = h * f_z(p_z[i - 1]);
        k_p_r_1 = h * f_p_r(r[i - 1], z[i - 1], p_psi[i - 1], model);
        k_p_z_1 = h * f_p_z(r[i - 1], z[i - 1], model);

        k_r_2 = h * f_r(p_r[i - 1] + 0.5 * k_p_r_1);
        k_psi_2 = h * f_psi(r[i - 1] + 0.5 * k_r_1, p_psi[i - 1]);
        k_z_2 = h * f_z(p_z[i - 1] + 0.5 * k_p_z_1);
        k_p_r_2 = h * f_p_r(
            r[i - 1] + 0.5 * k_r_1,
            z[i - 1] + 0.5 * k_z_1,
            p_psi[i - 1],
            model,
        );
        k_p_z_2 = h * f_p_z(r[i - 1] + 0.5 * k_r_1, z[i - 1] + 0.5 * k_z_1, model);

        k_r_3 = h * f_r(p_r[i - 1] + 0.5 * k_p_r_2);
        k_psi_3 = h * f_psi(r[i - 1] + 0.5 * k_r_2, p_psi[i - 1]);
        k_z_3 = h * f_z(p_z[i - 1] + 0.5 * k_p_z_2);
        k_p_r_3 = h * f_p_r(
            r[i - 1] + 0.5 * k_r_2,
            z[i - 1] + 0.5 * k_z_2,
            p_psi[i - 1],
            model,
        );
        k_p_z_3 = h * f_p_z(r[i - 1] + 0.5 * k_r_2, z[i - 1] + 0.5 * k_z_2, model);

        k_r_4 = h * f_r(p_r[i - 1] + k_p_r_3);
        k_psi_4 = h * f_psi(r[i - 1] + k_r_3, p_psi[i - 1]);
        k_z_4 = h * f_z(p_z[i - 1] + k_p_z_3);
        k_p_r_4 = h * f_p_r(r[i - 1] + k_r_3, z[i - 1] + k_z_3, p_psi[i - 1], model);
        k_p_z_4 = h * f_p_z(r[i - 1] + k_r_3, z[i - 1] + k_z_3, model);

        // Push the next values

        // Galactic Cylindrical system
        r.push(r[i - 1] + 1.0 / 6.0 * (k_r_1 + 2.0 * k_r_2 + 2.0 * k_r_3 + k_r_4));
        psi.push(psi[i - 1] + 1.0 / 6.0 * (k_psi_1 + 2.0 * k_psi_2 + 2.0 * k_psi_3 + k_psi_4));
        z.push(z[i - 1] + 1.0 / 6.0 * (k_z_1 + 2.0 * k_z_2 + 2.0 * k_z_3 + k_z_4));
        p_r.push(p_r[i - 1] + 1.0 / 6.0 * (k_p_r_1 + 2.0 * k_p_r_2 + 2.0 * k_p_r_3 + k_p_r_4));
        p_psi.push(p_psi[i - 1]);
        p_z.push(p_z[i - 1] + 1.0 / 6.0 * (k_p_z_1 + 2.0 * k_p_z_2 + 2.0 * k_p_z_3 + k_p_z_4));

        // Galactic Cartesian system
        x.push(r[i] * psi[i].cos());
        y.push(r[i] * psi[i].sin());

        // Total energy
        e.push(total_energy(
            r[i], psi[i], z[i], p_r[i], p_psi[i], p_z[i], model,
        ));
    }
}

impl Orbit {
    /// Integrate the orbit with the specified
    /// number of iterations `n` and the time step `h`
    /// using the fourth-order Runge-Kutta algorithm
    pub fn integrate(&mut self, model: &dyn Model, n: I, h: F, rev: bool) {
        // Save the number of integrations
        self.n = if rev { n * 2 } else { n };

        // Prepare aliases
        let (x_hc, _, y_hc, _, z_hc, _, u_hc, _, v_hc, _, w_hc, _) = self.hc_initials.unpack();
        let (x_gca, y_gca, z_gca, u_gca, v_gca, w_gca) = self.gca_initials.unpack();
        let (r_gcy, psi_gcy, z_gcy, dr_gcy, dpsi_gcy, dz_gcy) = self.gcy_initials.unpack();
        let (r_res, psi_res, z_res, p_r_res, p_psi_res, p_z_res, x_res, y_res, e_res, _, _) =
            self.results.unpack();

        // Convert the initial coordinates and velocities
        // from the Heliocentric Cartesian system to the
        // Galactic Cartesian system
        *x_gca = x_hc.to_gca_x();
        *y_gca = *y_hc;
        *z_gca = z_hc.to_gca_z();
        *u_gca = u_hc.to_gca_u();
        *v_gca = v_hc.to_gca_v();
        *w_gca = w_hc.to_gca_w();

        // Convert the initial coordinates and velocities
        // from the Galactic Cartesian system to the
        // Galactic Cylindrical system
        *r_gcy = (x_gca.powi(2) + y_gca.powi(2)).sqrt();
        *psi_gcy = y_gca.atan2(*x_gca);
        *z_gcy = *z_gca;
        *dr_gcy = -*u_gca * psi_gcy.cos() + *v_gca * psi_gcy.sin();
        *dpsi_gcy = ((*u_gca * psi_gcy.sin() + *v_gca * psi_gcy.cos()) / *r_gcy).to_kpc();
        *dz_gcy = *w_gca;

        // Prepare the results vectors

        // Solutions of the Lagrangian equations
        *r_res = Vec::<F>::with_capacity(self.n + 1);
        *psi_res = Vec::<F>::with_capacity(self.n + 1);
        *z_res = Vec::<F>::with_capacity(self.n + 1);
        *p_r_res = Vec::<F>::with_capacity(self.n + 1);
        *p_psi_res = Vec::<F>::with_capacity(self.n + 1);
        *p_z_res = Vec::<F>::with_capacity(self.n + 1);

        // Galactic Cartesian system
        *x_res = Vec::<F>::with_capacity(self.n + 1);
        *y_res = Vec::<F>::with_capacity(self.n + 1);

        // Total energy
        *e_res = Vec::<F>::with_capacity(self.n + 1);

        // Push the initial values

        // Solutions of the Lagrangian equations
        r_res.push(*r_gcy);
        psi_res.push(*psi_gcy);
        z_res.push(*z_gcy);
        p_r_res.push(dr_gcy.to_kpc_per_myr());
        p_psi_res.push((r_gcy.powi(2) * *dpsi_gcy).to_seconds());
        p_z_res.push(dz_gcy.to_kpc_per_myr());

        // Galactic Cartesian system
        x_res.push(r_res[0] * psi_res[0].cos());
        y_res.push(r_res[0] * psi_res[0].sin());

        // Total energy
        e_res.push(total_energy(
            r_res[0],
            psi_res[0],
            z_res[0],
            p_r_res[0],
            p_psi_res[0],
            p_z_res[0],
            model,
        ));

        // Integrate the orbit using the fourth-order Runge-Kutta algorithm
        rk4(self, model, h, 1, n);
        // Integrate in reverse if specified by the user
        if rev {
            rk4(self, model, -h, n + 1, self.n + 1);
        }
    }
}
