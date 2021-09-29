//! Integrate the orbit

// use ode_solvers::{Rk4, *};

use super::potentials::{miyamoto_nagai, navarro_frenk_white, plummer};
use crate::{orbit::Orbit, F, I, PADDING};

// type State = Vector6<F>;
// type Time = F;

/// Galactocentric distance of the Local Standard
/// of Rest around the Galactic center (kpc)
const G_R_SUN: F = 8.3;
/// Galactocentric linear velocity of the Local Standard
/// of Rest around the Galactic center (km/s)
const G_V_SUN: F = 244.0;
/// Height of the Sun above the Galactic plane (pc)
const G_H_SUN: F = 16.0;
/// The Sun’s peculiar velocity U with respect
/// to the Local Standard of Rest (km/s)
const P_U_SUN: F = 11.1;
/// The Sun’s peculiar velocity V with respect
/// to the Local Standard of Rest (km/s)
const P_V_SUN: F = 12.2;
/// The Sun’s peculiar velocity W with respect
/// to the Local Standard of Rest (km/s)
const P_W_SUN: F = 7.3;

/// `M` parameter of the Plummer potential
const M_B: F = 443.0;
/// `M` parameter of the Miyamoto & Nagai potential
const M_D: F = 2798.0;
/// `M` parameter of the Navarro-Frenk-White potential
const M_H: F = 12474.0;
/// `b` parameter of the Plummer potential
const B_B: F = 0.2672;
/// `a` parameter of the Miyamoto & Nagai potential
const A_D: F = 4.40;
/// `b` parameter of the Miyamoto & Nagai potential
const B_D: F = 0.3084;
/// `a` parameter of the Navarro-Frenk-White potential
const A_H: F = 7.7;

/// A constant to convert kpc/Myr to km/s
const KPC_PER_MYR_TO_KM_PER_S: F = 977.812_299_895_122;
/// A constant to convert km/s to kpc/Myr
const KM_PER_S_TO_KPC_PER_MYR: F = 0.001_022_73;
/// A constant to convert 100 km^2/s^2 to kpc^2/Myr^2
const HUNDREDS_KM_2_PER_S_2_TO_KPC_2_PER_MYR_2: F = 1.045_897_218_694_908e-4;
/// A constant to convert Myr to seconds
const MYR_TO_S: F = 3.155_695_2e13;
/// A constant to convert km to kpc
const KM_TO_KPC: F = 3.240_779_289_666_4e-17;

// const TEST: F = 1.0282;

/// Calculate the value of the Galactic potential
fn phi(r: F, z: F) -> F {
    plummer::phi(r, z, M_B, B_B)
        + miyamoto_nagai::phi(r, z, M_D, A_D, B_D)
        + navarro_frenk_white::phi(r, z, M_H, A_H)
    // * TEST
}

/// Calculate the value of the R derivative of the Galactic potential
fn phi_dr(r: F, z: F) -> F {
    (plummer::phi_dr(r, z, M_B, B_B)
        + miyamoto_nagai::phi_dr(r, z, M_D, A_D, B_D)
        + navarro_frenk_white::phi_dr(r, z, M_H, A_H))
        * HUNDREDS_KM_2_PER_S_2_TO_KPC_2_PER_MYR_2
    // * TEST
}

/// Calculate the value of the Z derivative of the Galactic potential
fn phi_dz(r: F, z: F) -> F {
    (plummer::phi_dz(r, z, M_B, B_B)
        + miyamoto_nagai::phi_dz(r, z, M_D, A_D, B_D)
        + navarro_frenk_white::phi_dz(r, z, M_H, A_H))
        * HUNDREDS_KM_2_PER_S_2_TO_KPC_2_PER_MYR_2
    // * TEST
}

/// Calculate the right-hand part of the
/// equation for the time derivative of `r`
fn f_r(p_r: F) -> F {
    p_r
}

/// Calculate the right-hand part of the
/// equation for the time derivative of `psi`
fn f_psi(r: F, p_psi: F) -> F {
    p_psi / r.powi(2)
}

/// Calculate the right-hand part of the
/// equation for the time derivative of `z`
fn f_z(p_z: F) -> F {
    p_z
}

/// Calculate the right-hand part of the
/// equation for the time derivative of `p_r`
fn f_p_r(r: F, z: F, p_psi: F) -> F {
    -phi_dr(r, z) + p_psi.powi(2) / r.powi(3)
}

/// Calculate the right-hand part of the
/// equation for the time derivative of `p_z`
fn f_p_z(r: F, z: F) -> F {
    -phi_dz(r, z)
}

/// Calculate total energy (km^2 / s^2)
fn total_energy(r: F, psi: F, z: F, p_r: F, p_psi: F, p_z: F) -> F {
    ((-(p_r * KPC_PER_MYR_TO_KM_PER_S * psi.cos()
        - p_psi / r * KPC_PER_MYR_TO_KM_PER_S * psi.sin())
        * psi.cos()
        + (p_r * KPC_PER_MYR_TO_KM_PER_S * psi.sin()
            + p_psi / r * KPC_PER_MYR_TO_KM_PER_S * psi.cos())
            * psi.sin())
    .powi(2)
        + ((p_r * KPC_PER_MYR_TO_KM_PER_S * psi.cos()
            - p_psi / r * KPC_PER_MYR_TO_KM_PER_S * psi.sin())
            * psi.sin()
            + (p_r * KPC_PER_MYR_TO_KM_PER_S * psi.sin()
                + p_psi / r * KPC_PER_MYR_TO_KM_PER_S * psi.cos())
                * psi.cos())
        .powi(2)
        + (p_z * KPC_PER_MYR_TO_KM_PER_S).powi(2))
        / 2.0
        + phi(r, z) * 100.0
}

impl Orbit {
    /// Integrate the orbit with the specified
    /// number of iterations `n` and the time step `h`
    /// using the fourth-order Runge-Kutta algorithm
    pub fn integrate(&mut self, n: I, h: F) {
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
        self.gcy_initials.dr = self.gca_initials.u * self.gcy_initials.psi.cos()
            + self.gca_initials.v * self.gcy_initials.psi.sin();
        self.gcy_initials.dpsi = (-self.gca_initials.u * self.gcy_initials.psi.sin()
            + self.gca_initials.v * self.gcy_initials.psi.cos())
            / self.gcy_initials.r
            * KM_TO_KPC;
        self.gcy_initials.dz = self.gca_initials.w;

        // {
        //     let r = 1.0;
        //     let z = 0.0;
        //     println!(
        //         "vel: {}",
        //         F::sqrt(
        //             (plummer::phi_dr(r, z, M_B, B_B)
        //                 + miyamoto_nagai::phi_dr(r, z, M_D, A_D, B_D)
        //                 + navarro_frenk_white::phi_dr(r, z, M_H, A_H))
        //                 * 100.0
        //         )
        //     );
        // }

        // println!(
        //     "{}",
        //     ((-(self.gcy_initials.r_vel * self.gcy_initials.psi.cos()
        //         - self.gcy_initials.psi_vel
        //             * (self.gcy_initials.r / KM_TO_PC)
        //             * self.gcy_initials.psi.sin())
        //         * self.gcy_initials.psi.cos()
        //         + (self.gcy_initials.r_vel * self.gcy_initials.psi.sin()
        //             + self.gcy_initials.psi_vel
        //                 * (self.gcy_initials.r / KM_TO_PC)
        //                 * self.gcy_initials.psi.cos())
        //             * self.gcy_initials.psi.sin())
        //     .powi(2)
        //         + ((self.gcy_initials.r_vel * self.gcy_initials.psi.cos()
        //             - self.gcy_initials.psi_vel
        //                 * (self.gcy_initials.r / KM_TO_PC)
        //                 * self.gcy_initials.psi.sin())
        //             * self.gcy_initials.psi.sin()
        //             + (self.gcy_initials.r_vel * self.gcy_initials.psi.sin()
        //                 + self.gcy_initials.psi_vel
        //                     * (self.gcy_initials.r / KM_TO_PC)
        //                     * self.gcy_initials.psi.cos())
        //                 * self.gcy_initials.psi.cos())
        //         .powi(2)
        //         + self.gcy_initials.z_vel.powi(2))
        //         / 2.0
        //         + phi(self.gcy_initials.r, self.gcy_initials.z) * 100.0
        // );

        // Prepare the vectors for the solutions

        // Galactic Cylindrical system
        let mut r = Vec::<F>::with_capacity(n + 1);
        let mut psi = Vec::<F>::with_capacity(n + 1);
        let mut z = Vec::<F>::with_capacity(n + 1);
        let mut p_r = Vec::<F>::with_capacity(n + 1);
        let mut p_psi = Vec::<F>::with_capacity(n + 1);
        let mut p_z = Vec::<F>::with_capacity(n + 1);

        // Galactic Cartesian system
        let mut x = Vec::<F>::with_capacity(n + 1);
        let mut y = Vec::<F>::with_capacity(n + 1);

        // Total energy
        let mut e = Vec::<F>::with_capacity(n + 1);
        let mut vel = Vec::<F>::with_capacity(n + 1);
        let mut phi_v = Vec::<F>::with_capacity(n + 1);

        // Push the initial values

        // Galactic Cylindrical system
        r.push(self.gcy_initials.r);
        psi.push(self.gcy_initials.psi);
        z.push(self.gcy_initials.z);
        p_r.push(self.gcy_initials.dr * KM_PER_S_TO_KPC_PER_MYR);
        p_psi.push(self.gcy_initials.r.powi(2) * self.gcy_initials.dpsi * MYR_TO_S);
        p_z.push(self.gcy_initials.dz * KM_PER_S_TO_KPC_PER_MYR);

        // Galactic Cartesian system
        x.push(r[0] * psi[0].cos());
        y.push(r[0] * psi[0].sin());

        // println!(
        //     "gca:\n{}\n{}\n{}\n{}\n{}\n{}\ngcy:\n{}\n{}\n{}\n{}\n{}\n{}",
        //     self.gca_initials.x,
        //     self.gca_initials.y,
        //     self.gca_initials.z,
        //     self.gca_initials.u,
        //     self.gca_initials.v,
        //     self.gca_initials.w,
        //     self.gcy_initials.r,
        //     self.gcy_initials.psi,
        //     self.gcy_initials.z,
        //     self.gcy_initials.r_vel,
        //     self.gcy_initials.psi_vel,
        //     self.gcy_initials.z_vel,
        // );

        // Total energy
        e.push(total_energy(r[0], psi[0], z[0], p_r[0], p_psi[0], p_z[0]));
        vel.push(
            ((-(p_r[0] * KPC_PER_MYR_TO_KM_PER_S * psi[0].cos()
                - p_psi[0] / r[0] * KPC_PER_MYR_TO_KM_PER_S * psi[0].sin())
                * psi[0].cos()
                + (p_r[0] * KPC_PER_MYR_TO_KM_PER_S * psi[0].sin()
                    + p_psi[0] / r[0] * KPC_PER_MYR_TO_KM_PER_S * psi[0].cos())
                    * psi[0].sin())
            .powi(2)
                + ((p_r[0] * KPC_PER_MYR_TO_KM_PER_S * psi[0].cos()
                    - p_psi[0] / r[0] * KPC_PER_MYR_TO_KM_PER_S * psi[0].sin())
                    * psi[0].sin()
                    + (p_r[0] * KPC_PER_MYR_TO_KM_PER_S * psi[0].sin()
                        + p_psi[0] / r[0] * KPC_PER_MYR_TO_KM_PER_S * psi[0].cos())
                        * psi[0].cos())
                .powi(2)
                + (p_z[0] * KPC_PER_MYR_TO_KM_PER_S).powi(2))
                / 2.0,
        );
        phi_v.push(phi_dz(r[0], z[0]));

        println!("\n{:1$}E at the start: {2:.16e}", "", PADDING + 2, e[0]);

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

        // Integrate the orbit using the fourth-order Runge-Kutta algorithm
        for i in 1..=n {
            // Calculate the increments for all integrable variables
            // (except for `p_psi`, since its derivative is 0)

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

            // println!(
            //     "{}",
            //     ((-(p_r[i - 1] * KPC_PER_MYR_TO_KM_PER_S * psi[i - 1].cos()
            //         - p_psi[i - 1] * S_TO_MYR / r[i - 1] / KM_TO_PC * psi[i - 1].sin())
            //         * psi[i - 1].cos()
            //         + (p_r[i - 1] * KPC_PER_MYR_TO_KM_PER_S * psi[i - 1].sin()
            //             + p_psi[i - 1] * S_TO_MYR / r[i - 1] / KM_TO_PC * psi[i - 1].cos())
            //             * psi[i - 1].sin())
            //     .powi(2)
            //         + ((p_r[i - 1] * KPC_PER_MYR_TO_KM_PER_S * psi[i - 1].cos()
            //             - p_psi[i - 1] * S_TO_MYR / r[i - 1] / KM_TO_PC * psi[i - 1].sin())
            //             * psi[i - 1].sin()
            //             + (p_r[i - 1] * KPC_PER_MYR_TO_KM_PER_S * psi[i - 1].sin()
            //                 + p_psi[i - 1] * S_TO_MYR / r[i - 1] / KM_TO_PC * psi[i - 1].cos())
            //                 * psi[i - 1].cos())
            //         .powi(2)
            //         + (p_z[i - 1] * KPC_PER_MYR_TO_KM_PER_S).powi(2))
            //         / 2.0
            //         + phi(r[i - 1], z[i - 1]) * 100.0
            // );

            // Calculate the next value

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
            e.push(total_energy(r[i], psi[i], z[i], p_r[i], p_psi[i], p_z[i]));
            vel.push(
                ((-(p_r[i] * KPC_PER_MYR_TO_KM_PER_S * psi[i].cos()
                    - p_psi[i] / r[i] * KPC_PER_MYR_TO_KM_PER_S * psi[i].sin())
                    * psi[i].cos()
                    + (p_r[i] * KPC_PER_MYR_TO_KM_PER_S * psi[i].sin()
                        + p_psi[i] / r[i] * KPC_PER_MYR_TO_KM_PER_S * psi[i].cos())
                        * psi[i].sin())
                .powi(2)
                    + ((p_r[i] * KPC_PER_MYR_TO_KM_PER_S * psi[i].cos()
                        - p_psi[i] / r[i] * KPC_PER_MYR_TO_KM_PER_S * psi[i].sin())
                        * psi[i].sin()
                        + (p_r[i] * KPC_PER_MYR_TO_KM_PER_S * psi[i].sin()
                            + p_psi[i] / r[i] * KPC_PER_MYR_TO_KM_PER_S * psi[i].cos())
                            * psi[i].cos())
                    .powi(2)
                    + (p_z[i] * KPC_PER_MYR_TO_KM_PER_S).powi(2))
                    / 2.0,
            );
            phi_v.push(phi_dz(r[i], z[i]));
        }

        println!("{:1$}E at the end:   {2:.16e}", "", PADDING + 2, e[n]);
        println!(
            "{:1$}Emax - Emin:    {2:.16e}\n",
            "",
            PADDING + 2,
            e.iter().copied().reduce(F::max).unwrap() - e.iter().copied().reduce(F::min).unwrap()
        );

        // let system = LagrangianEquations {};

        // let y0 = State::new(
        //     self.gcy_initials.r,
        //     self.gcy_initials.psi,
        //     self.gcy_initials.z,
        //     self.gcy_initials.r_vel * KM_PER_S_TO_KPC_PER_MYR,
        //     self.gcy_initials.r.powi(2) * self.gcy_initials.psi_vel * MYR_TO_S,
        //     self.gcy_initials.z_vel * KM_PER_S_TO_KPC_PER_MYR,
        // );

        // let mut stepper = Rk4::new(system, 0.0, y0, n as F * h, h);
        // let res = stepper.integrate();

        // match res {
        //     Ok(_) => {
        //         let states = stepper.y_out();
        //         // Save the results
        //         self.integrated.r = states.iter().map(|x| x[0]).collect();
        //         self.integrated.psi = states.iter().map(|x| x[1]).collect();
        //         self.integrated.z = states.iter().map(|x| x[2]).collect();
        //         let p_r: Vec<F> = states.iter().map(|x| x[3]).collect();
        //         let p_psi: Vec<F> = states.iter().map(|x| x[4]).collect();
        //         let p_z: Vec<F> = states.iter().map(|x| x[5]).collect();
        //         let mut e = Vec::<F>::with_capacity(n + 1);
        //         let mut x = Vec::<F>::with_capacity(n + 1);
        //         let mut y = Vec::<F>::with_capacity(n + 1);
        //         for i in 0..=n {
        //             e.push(total_energy(
        //                 self.integrated.r[i],
        //                 self.integrated.psi[i],
        //                 self.integrated.z[i],
        //                 p_r[i],
        //                 p_psi[i],
        //                 p_z[i],
        //             ));
        //             x.push(self.integrated.r[i] * self.integrated.psi[i].cos());
        //             y.push(self.integrated.r[i] * self.integrated.psi[i].sin());
        //         }
        //         self.integrated.x = x;
        //         self.integrated.y = y;
        //         self.integrated.e = e;
        //     }
        //     Err(_) => eprintln!(":shrug:"),
        // }

        // Save the results

        // Galactic Cylindrical system
        self.integrated.r = r;
        self.integrated.z = z;

        // Galactic Cartesian system
        self.integrated.x = x;
        self.integrated.y = y;

        // Total energy
        self.integrated.e = e;
        self.integrated.vel = vel;
        self.integrated.phi = phi_v;
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
