//! Integrate the orbit

use ode_solvers::{Rk4, *};

use super::potentials::{miyamoto_nagai, navarro_frenk_white, plummer};
use crate::orbit::{Orbit, F};

type State = Vector6<F>;
type Time = F;

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

const M_B: F = 443.0;
const M_D: F = 2798.0;
const M_H: F = 12474.0;
const B_B: F = 0.2672;
const A_D: F = 4.40;
const B_D: F = 0.3084;
const A_H: F = 7.7;

/// Calculate the value of the Galactic potential
fn phi(r: F, z: F) -> F {
    plummer::phi(r, z, M_B, B_B)
        + miyamoto_nagai::phi(r, z, M_D, A_D, B_D)
        + navarro_frenk_white::phi(r, z, M_H, A_H)
}

/// Calculate the value of the R derivative of the Galactic potential
fn phi_dr(r: F, z: F) -> F {
    plummer::phi_dr(r, z, M_B, B_B)
        + miyamoto_nagai::phi_dr(r, z, M_D, A_D, B_D)
        + navarro_frenk_white::phi_dr(r, z, M_H, A_H)
}

/// Calculate the value of the Z derivative of the Galactic potential
fn phi_dz(r: F, z: F) -> F {
    plummer::phi_dz(r, z, M_B, B_B)
        + miyamoto_nagai::phi_dz(r, z, M_D, A_D, B_D)
        + navarro_frenk_white::phi_dz(r, z, M_H, A_H)
}

/// Calculate the right-hand part of the equation for r'
fn f_r(p_r: F) -> F {
    p_r * 0.00102273
}

/// Calculate the right-hand part of the equation for psi'
fn f_psi(r: F, p_psi: F) -> F {
    p_psi / r.powi(2)
}

/// Calculate the right-hand part of the equation for z'
fn f_z(p_z: F) -> F {
    p_z * 0.00102273
}

/// Calculate the right-hand part of the equation for p_r'
fn f_p_r(r: F, z: F, p_psi: F) -> F {
    -phi_dr(r, z) * 1.04598e-4 + p_psi.powi(2) / r.powi(3) / CONVERT
}

/// Calculate the right-hand part of the equation for p_z'
fn f_p_z(r: F, z: F) -> F {
    -phi_dz(r, z) * 1.04598e-4
}

const CONVERT: F = 3.240_779_289_666_4e-17;

impl Orbit {
    /// Integrate the orbit with the specified
    /// number of iteration `n` and the time step `h`
    /// using the fourth-order Runge-Kutta algorithm
    pub fn integrate(&mut self, n: usize, h: F) {
        // Save the number of integrations
        self.n = n;

        // Convert the initial coordinates and velocities
        // from the Heliocentric Cartesian system to the
        // Galactic Cartesian system
        self.gca_initials.x = G_R_SUN - self.hc_initials.x;
        self.gca_initials.y = self.hc_initials.y;
        self.gca_initials.z = self.hc_initials.z + G_H_SUN / 1000.0;
        self.gca_initials.u = self.hc_initials.u + P_U_SUN;
        self.gca_initials.v = self.hc_initials.v + P_V_SUN + G_V_SUN;
        self.gca_initials.w = self.hc_initials.w + P_W_SUN;

        // Convert the initial coordinates and velocities
        // from the Galactic Cartesian system to the
        // Galactic Cylindrical system
        self.gcy_initials.r = (self.gca_initials.x.powi(2) + self.gca_initials.y.powi(2)).sqrt();
        self.gcy_initials.psi = self.gca_initials.y.atan2(self.gca_initials.x);
        self.gcy_initials.z = self.gca_initials.z;
        self.gcy_initials.r_vel =
            (self.gca_initials.u.powi(2) + self.gca_initials.v.powi(2)).sqrt();
        self.gcy_initials.psi_vel = self.gca_initials.v.atan2(self.gca_initials.u);
        self.gcy_initials.z_vel = self.gca_initials.w;

        // println!(
        //     "{}",
        //     -self.gcy_initials.psi_vel.cos()
        //         * self.gcy_initials.psi.cos()
        //         * self.gcy_initials.r_vel
        //         + self.gcy_initials.psi_vel.sin()
        //             * self.gcy_initials.psi.sin()
        //             * self.gcy_initials.r_vel
        // );

        // println!(
        //     "{}",
        //     self.gcy_initials.psi_vel.cos() * self.gcy_initials.psi.sin() * self.gcy_initials.r_vel
        //         + self.gcy_initials.psi_vel.sin()
        //             * self.gcy_initials.psi.cos()
        //             * self.gcy_initials.r_vel
        // );

        // println!(
        //     "{}",
        //     ((-self.gcy_initials.psi_vel.cos()
        //         * self.gcy_initials.psi.cos()
        //         * self.gcy_initials.r_vel
        //         + self.gcy_initials.psi_vel.sin()
        //             * self.gcy_initials.psi.sin()
        //             * self.gcy_initials.r_vel)
        //         .powi(2)
        //         + (self.gcy_initials.psi_vel.cos()
        //             * self.gcy_initials.psi.sin()
        //             * self.gcy_initials.r_vel
        //             + self.gcy_initials.psi_vel.sin()
        //                 * self.gcy_initials.psi.cos()
        //                 * self.gcy_initials.r_vel)
        //             .powi(2)
        //         + self.gcy_initials.z_vel.powi(2))
        //     .sqrt()
        // );

        println!(
            "{}",
            ((-self.gcy_initials.psi_vel.cos()
                * self.gcy_initials.psi.cos()
                * self.gcy_initials.r_vel
                + self.gcy_initials.psi_vel.sin()
                    * self.gcy_initials.psi.sin()
                    * self.gcy_initials.r_vel)
                .powi(2)
                + (self.gcy_initials.psi_vel.cos()
                    * self.gcy_initials.psi.sin()
                    * self.gcy_initials.r_vel
                    + self.gcy_initials.psi_vel.sin()
                        * self.gcy_initials.psi.cos()
                        * self.gcy_initials.r_vel)
                    .powi(2)
                + self.gcy_initials.z_vel.powi(2))
                / 2.0
                + phi(self.gcy_initials.r, self.gcy_initials.z) * 100.0
        );

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
        p_r.push(self.gcy_initials.r_vel);
        p_psi.push(self.gcy_initials.r.powi(2) * self.gcy_initials.psi_vel);
        p_z.push(self.gcy_initials.z_vel);

        // println!(
        //     "{}",
        //     -(p_psi[0] / r[0].powi(2)).cos() * psi[0].cos() * p_r[0] / CONVERT
        //         + (p_psi[0] / r[0].powi(2)).sin() * psi[0].sin() * p_r[0] / CONVERT
        // );

        // println!(
        //     "{}",
        //     -(p_psi[0] / r[0].powi(2)).cos() * psi[0].sin() * p_r[0] / CONVERT
        //         + (p_psi[0] / r[0].powi(2)).sin() * psi[0].cos() * p_r[0] / CONVERT
        // );

        // println!(
        //     "{}",
        //     ((-(p_psi[0] / r[0].powi(2)).cos() * psi[0].cos() * p_r[0] / CONVERT
        //         + (p_psi[0] / r[0].powi(2)).sin() * psi[0].sin() * p_r[0] / CONVERT)
        //         .powi(2)
        //         + (-(p_psi[0] / r[0].powi(2)).cos() * psi[0].sin() * p_r[0] / CONVERT
        //             + (p_psi[0] / r[0].powi(2)).sin() * psi[0].cos() * p_r[0] / CONVERT)
        //             .powi(2)
        //         + (p_z[0] / CONVERT).powi(2))
        //     .sqrt()
        // );

        println!(
            "{}",
            ((-(p_psi[0] / r[0].powi(2)).cos() * psi[0].cos() * p_r[0]
                + (p_psi[0] / r[0].powi(2)).sin() * psi[0].sin() * p_r[0])
                .powi(2)
                + (-(p_psi[0] / r[0].powi(2)).cos() * psi[0].sin() * p_r[0]
                    + (p_psi[0] / r[0].powi(2)).sin() * psi[0].cos() * p_r[0])
                    .powi(2)
                + (p_z[0]).powi(2))
                / 2.0
                + phi(r[0], z[0]) * 100.0
        );

        // Prepare buffers for the coefficients
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

        // Integrate the orbit using the fourth-order Runge-Kutta algorithm
        for i in 1..=n {
            k_r_1 = h * f_r(p_r[i - 1]);
            k_psi_1 = h * f_psi(r[i - 1], p_psi[i - 1]);
            k_z_1 = h * f_z(p_z[i - 1]);
            k_p_r_1 = h * f_p_r(r[i - 1], z[i - 1], p_psi[i - 1]);
            k_p_z_1 = h * f_p_z(r[i - 1], z[i - 1]);

            k_r_2 = h * f_r(p_r[i - 1] + 0.5 * k_p_r_1);
            k_psi_2 = h * f_psi(r[i - 1] + 0.5 * k_r_1, p_psi[i - 1]);
            k_z_2 = h * f_z(p_z[i - 1] + 0.5 * k_p_z_1);
            k_p_r_2 = h * f_p_r(r[i - 1] + 0.5 * k_r_1, z[i - 1] + 0.5 * k_z_1, p_psi[i - 1]);
            k_p_z_2 = h * f_p_z(r[i - 1] + 0.5 * k_r_1, z[i - 1] + 0.5 * k_z_1);

            k_r_3 = h * f_r(p_r[i - 1] + 0.5 * k_p_r_2);
            k_psi_3 = h * f_psi(r[i - 1] + 0.5 * k_r_2, p_psi[i - 1]);
            k_z_3 = h * f_z(p_z[i - 1] + 0.5 * k_p_z_2);
            k_p_r_3 = h * f_p_r(r[i - 1] + 0.5 * k_r_2, z[i - 1] + 0.5 * k_z_2, p_psi[i - 1]);
            k_p_z_3 = h * f_p_z(r[i - 1] + 0.5 * k_r_2, z[i - 1] + 0.5 * k_z_2);

            k_r_4 = h * f_r(p_r[i - 1] + k_p_r_3);
            k_psi_4 = h * f_psi(r[i - 1] + k_r_3, p_psi[i - 1]);
            k_z_4 = h * f_z(p_z[i - 1] + k_p_z_3);
            k_p_r_4 = h * f_p_r(r[i - 1] + k_r_3, z[i - 1] + k_z_3, p_psi[i - 1]);
            k_p_z_4 = h * f_p_z(r[i - 1] + k_r_3, z[i - 1] + k_z_3);

            println!(
                "{:.18e}",
                ((-(p_psi[i - 1] / r[i - 1].powi(2)).cos() * psi[i - 1].cos() * p_r[i - 1]
                    + (p_psi[i - 1] / r[i - 1].powi(2)).sin() * psi[i - 1].sin() * p_r[i - 1])
                    .powi(2)
                    + (-(p_psi[i - 1] / r[i - 1].powi(2)).cos() * psi[i - 1].sin() * p_r[i - 1]
                        + (p_psi[i - 1] / r[i - 1].powi(2)).sin() * psi[i - 1].cos() * p_r[i - 1])
                        .powi(2)
                    + (p_z[i - 1]).powi(2))
                    / 2.0
                    + phi(r[i - 1], z[i - 1]) * 100.0
            );

            r.push(r[i - 1] + 1.0 / 6.0 * (k_r_1 + 2.0 * k_r_2 + 2.0 * k_r_3 + k_r_4));
            psi.push(psi[i - 1] + 1.0 / 6.0 * (k_psi_1 + 2.0 * k_psi_2 + 2.0 * k_psi_3 + k_psi_4));
            z.push(z[i - 1] + 1.0 / 6.0 * (k_z_1 + 2.0 * k_z_2 + 2.0 * k_z_3 + k_z_4));
            p_r.push(p_r[i - 1] + 1.0 / 6.0 * (k_p_r_1 + 2.0 * k_p_r_2 + 2.0 * k_p_r_3 + k_p_r_4));
            p_psi.push(p_psi[i - 1]);
            p_z.push(p_z[i - 1] + 1.0 / 6.0 * (k_p_z_1 + 2.0 * k_p_z_2 + 2.0 * k_p_z_3 + k_p_z_4));
        }

        // println!(
        //     "{:.18e}",
        //     ((-(p_psi[0] / r[0].powi(2)).cos() * psi[0].cos() * p_r[0]
        //         + (p_psi[0] / r[0].powi(2)).sin() * psi[0].sin() * p_r[0])
        //         .powi(2)
        //         + (-(p_psi[0] / r[0].powi(2)).cos() * psi[0].sin() * p_r[0]
        //             + (p_psi[0] / r[0].powi(2)).sin() * psi[0].cos() * p_r[0])
        //             .powi(2)
        //         + (p_z[0]).powi(2))
        //         / 2.0
        //         + phi(r[0], z[0]) * 100.0
        // );

        // println!(
        //     "{:.18e}",
        //     ((-(p_psi[1] / r[1].powi(2)).cos() * psi[1].cos() * p_r[1]
        //         + (p_psi[1] / r[1].powi(2)).sin() * psi[1].sin() * p_r[1])
        //         .powi(2)
        //         + (-(p_psi[1] / r[1].powi(2)).cos() * psi[1].sin() * p_r[1]
        //             + (p_psi[1] / r[1].powi(2)).sin() * psi[1].cos() * p_r[1])
        //             .powi(2)
        //         + (p_z[1]).powi(2))
        //         / 2.0
        //         + phi(r[1], z[1]) * 100.0
        // );

        // let system = LagrangianEquations {};

        // let y0 = State::new(
        //     self.gcy_initials.r,
        //     self.gcy_initials.psi,
        //     self.gcy_initials.z,
        //     self.gcy_initials.r_vel * 3.240_779_289_666_4e-17,
        //     self.gcy_initials.r.powi(2) * self.gcy_initials.psi_vel,
        //     self.gcy_initials.z_vel * 3.240_779_289_666_4e-17,
        // );

        // let mut stepper = Rk4::new(system, 0.0, y0, n as F * h, h);
        // let res = stepper.integrate();

        // match res {
        //     Ok(_) => {
        //         let states = stepper.y_out();
        //         // Save the results
        //         self.gcy_integrated.r = states.iter().map(|x| x[0]).collect();
        //         self.gcy_integrated.psi = states.iter().map(|x| x[1]).collect();
        //         self.gcy_integrated.z = states.iter().map(|x| x[2]).collect();
        //     }
        //     Err(_) => eprintln!(":shrug:"),
        // }

        // Save the results
        self.gcy_integrated.r = r;
        self.gcy_integrated.psi = psi;
        self.gcy_integrated.z = z;
    }
}

// struct LagrangianEquations {}

// impl ode_solvers::System<State> for LagrangianEquations {
//     fn system(&self, _t: Time, y: &State, dy: &mut State) {
//         // r
//         dy[0] = f_r(y[3]);
//         // psi
//         dy[1] = f_psi(y[0], y[4]);
//         // z
//         dy[2] = f_z(y[5]);
//         // p_r
//         dy[3] = f_p_r(y[0], y[2], y[4]);
//         // p_psi
//         dy[4] = 0.0;
//         // p_z
//         dy[5] = f_p_z(y[0], y[2]);
//     }

//     // Stop the integration if x exceeds 25,500 km. Optional
//     // fn solout(&mut self, _t: Time, y: &State, _dy: &State) -> bool {
//     // y[0] > 25500.
//     // }
// }
