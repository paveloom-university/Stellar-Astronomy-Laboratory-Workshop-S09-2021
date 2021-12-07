//! This binary crate is a command-line utility to integrate the orbits of
//! globular clusters by providing their initial positions and velocities
//! and choosing a model of the Galaxy potential. It is mostly based on the
//! paper of Bajkova & Bobylev ([2020, v1](https://arxiv.org/abs/2008.13624v1)).
//!
//! Module documentations:
//! - [Command-line interface](crate::cli)
//! - [Description of units conversions before and
//!   after the integration](crate::orbit::calc::integrate)
//! - [A list of available potentials and their
//!   derivatives](crate::orbit::calc::potentials)
//! - [A list of available models](crate::orbit::calc::models)

use std::path::PathBuf;

pub mod cli;
pub mod orbit;

use orbit::calc::models::MODELS;
use orbit::Orbit;

/// The floating point type used across the program
pub type F = f64;

/// The number $ \pi $ casted to the floating point type
pub const PI: F = std::f64::consts::PI;

/// The integer type used across the program
pub type I = usize;

/// Padding used when printing the output
const PADDING: usize = 5;

/// Run the program
fn main() {
    // Define the CLI of the program, get the matched arguments
    let matches = cli::get_matches();

    // Get the model
    let model = MODELS[matches.value_of("model").unwrap().parse::<I>().unwrap() - 1];
    // Get the number of iterations
    let n = matches.value_of("n").unwrap().parse::<I>().unwrap();
    // Get the time step
    let h = matches.value_of("h").unwrap().parse::<F>().unwrap();
    // Get path to the results directory
    let output = PathBuf::from(matches.value_of("output").unwrap());
    // Get the results fields to write
    let fields: Vec<&str> = matches.values_of("fields").unwrap().collect();
    // Get values of the input files
    let files: Vec<&str> = matches.values_of("file(s)").unwrap().collect();
    // Get switch value of the reverse mode
    let rev = matches.occurrences_of("rev") != 0;
    // Get switch value of the simulation mode
    let sim = matches.occurrences_of("sim") != 0;
    // Get the number of Monte Carlo simulations
    let s = if sim {
        matches.value_of("s").unwrap().parse::<I>().unwrap()
    } else {
        0
    };

    println!("\n{:1$}> Parsing the files...", "", PADDING);

    // Parse the files and get the orbits
    let (orbits, log) = Orbit::load(files);
    // Print the parsing log
    println!("{}", log.format(PADDING + 2));

    /// Print general information before the integration process begins
    macro_rules! info {
        () => {
            if rev {
                println!("{:1$}Reverse mode enabled.\n", "", PADDING + 2);
            } else if sim {
                println!("{:1$}Simulation mode enabled.\n", "", PADDING + 2);
            }
        };
    }

    // Check how many orbits (initial coordinates and velocities) were parsed
    if orbits.is_empty() {
        println!("{:1$}> No orbits were parsed. Exiting.", "", PADDING);
    } else if orbits.len() == 1 {
        println!("{:1$}> One orbit was parsed.\n", "", PADDING);
        info!();
    } else {
        println!(
            "{:1$}> {2} orbits were parsed.\n",
            "",
            PADDING,
            orbits.len()
        );
        info!();
    }

    // Integrate each orbit and write the results
    for ref mut orbit in orbits {
        if sim {
            println!(
                "{:1$}> Integrating the simulated orbits of {2}...",
                "",
                PADDING,
                orbit.id()
            );
            orbit.simulate(model, &output, &fields, s, n, h);
        } else {
            println!(
                "{:1$}> Integrating the orbit of {2}...",
                "",
                PADDING,
                orbit.id()
            );
            orbit.integrate(model, n, h, rev);
            orbit.write(&output, &fields, false);
        }
    }

    println!();
}
