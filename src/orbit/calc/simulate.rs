//! This module provides routines for the Monte Carlo simulation

use rand::prelude::*;
use rand_distr::{Distribution, Normal};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::path::Path;

use crate::{Orbit, F, I};

/// Calculate the distance using R and Z coordinates
fn dist(r: F, z: F) -> F {
    F::sqrt(r.powi(2) + z.powi(2))
}

/// Find and return the apocentric and pericentric distances
fn dist_extrema(o: &Orbit) -> (F, F) {
    // Get the initial values
    let mut apo = dist(o.results.r[0], o.results.z[0]);
    let mut peri = apo;

    // For every other pair
    for j in 1..=o.n {
        // Calculate the distance
        let d = dist(o.results.r[j], o.results.z[j]);
        // Depending on the distance, update either
        // the apocentric or the pericentric distance
        if d > apo {
            apo = d;
        } else if d < peri {
            peri = d;
        }
    }

    (apo, peri)
}

impl Orbit {
    /// Apply the Monte Carlo method:
    /// simulate sampling from Normal distribution
    /// using the mean and standard deviation of the
    /// input data, use that as the input to the
    /// [`integrate`](super::integrate) routine `s`
    /// times with the number of iterations `n` and
    /// the time step `h`; write the coordinates and
    /// total energy in the `folder` if they're specified
    /// in the `fields`
    ///
    /// # Panics
    /// Can panic if standard deviations aren't finite.
    /// This should be handled in the [`load`](Orbit#method.load) routine.
    pub fn simulate(&mut self, folder: &Path, fields: &[&str], s: I, n: I, h: F) {
        // Initialize the random generator
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);

        // Initialize the Normal distributions
        let normal_x = Normal::new(self.hc_initials.x, self.hc_initials.x_err).unwrap();
        let normal_y = Normal::new(self.hc_initials.y, self.hc_initials.y_err).unwrap();
        let normal_z = Normal::new(self.hc_initials.z, self.hc_initials.z_err).unwrap();
        let normal_u = Normal::new(self.hc_initials.u, self.hc_initials.u_err).unwrap();
        let normal_v = Normal::new(self.hc_initials.v, self.hc_initials.v_err).unwrap();
        let normal_w = Normal::new(self.hc_initials.w, self.hc_initials.w_err).unwrap();

        // Prepare the results vectors
        self.results.apo = Vec::<F>::with_capacity(s + 1);
        self.results.peri = Vec::<F>::with_capacity(s + 1);

        // Define what fields will be sequentially saved
        let integrated_fields = fields
            .iter()
            .copied()
            .filter(|f| f == &"r" || f == &"z" || f == &"x" || f == &"y" || f == &"e")
            .collect::<Vec<&str>>();
        let distances_fields = fields
            .iter()
            .copied()
            .filter(|f| f == &"apo" || f == &"peri")
            .collect::<Vec<&str>>();

        // Integrate with initial values first
        self.integrate(n, h, false);

        // Find the apocentric and pericentric distances
        let (apo, peri) = dist_extrema(self);

        // Save the distances
        self.results.apo.push(apo);
        self.results.peri.push(peri);

        // Write the coordinates and total energy
        self.write(folder, &integrated_fields, false);

        // For each iteration
        for _ in 0..s {
            // Simulate the input data
            self.hc_initials.x = normal_x.sample(&mut rng);
            self.hc_initials.y = normal_y.sample(&mut rng);
            self.hc_initials.z = normal_z.sample(&mut rng);
            self.hc_initials.u = normal_u.sample(&mut rng);
            self.hc_initials.v = normal_v.sample(&mut rng);
            self.hc_initials.w = normal_w.sample(&mut rng);

            // Integrate the orbit
            self.integrate(n, h, false);

            // Find the apocentric and pericentric distances
            let (apo, peri) = dist_extrema(self);

            // Save the distances
            self.results.apo.push(apo);
            self.results.peri.push(peri);

            // Append the coordinates and total energy
            self.write(folder, &integrated_fields, true);
        }

        // Write the distances
        self.write(folder, &distances_fields, false);
    }
}
