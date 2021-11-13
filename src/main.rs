//! This binary crate is a command-line utility to integrate the orbits of
//! globular clusters by providing their initial positions and velocities
//! and choosing a model of the Galaxy potential. It is based on the paper
//! of Bajkova & Bobylev ([2020, v1](https://arxiv.org/abs/2008.13624v1)).
//!
//! Module documentations:
//! - [Description of the units conversions before and
//! after the integration](crate::orbit::calc::integrate)
//! - [A list of available potentials and their
//! derivatives](crate::orbit::calc::potentials)

use std::path::PathBuf;

mod cli;
pub mod orbit;

use orbit::calc::models::M1;
use orbit::Orbit;

/// The floating point type used across the program
pub type F = f64;

/// The integer type used across the program
pub type I = usize;

/// Padding used when printing the output
const PADDING: usize = 5;

/// Run the program
fn main() {
    // Define the CLI of the program, get the matched arguments
    let matches = cli::get_matches();

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
            println!("{:1$}h = {2}", "", PADDING + 2, h);
            println!("{:1$}n = {2}", "", PADDING + 2, n);
            if rev {
                println!("{:1$}\nReverse mode enabled.", "", PADDING + 2);
            } else if sim {
                println!("{:1$}s = {2}\n", "", PADDING + 2, s);
                println!("{:1$}Simulation mode enabled.", "", PADDING + 2);
            }
            println!();
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

    // Choose the model
    let m = M1 {};

    // Integrate each orbit and write the results
    for ref mut orbit in orbits {
        if sim {
            println!(
                "{:1$}> Integrating the simulated orbits of {2}...",
                "",
                PADDING,
                orbit.id()
            );
            orbit.simulate(&m, &output, &fields, s, n, h);
        } else {
            println!(
                "{:1$}> Integrating the orbit of {2}...",
                "",
                PADDING,
                orbit.id()
            );
            orbit.integrate(&m, n, h, rev);
            orbit.write(&output, &fields, false);
        }
    }

    println!();
}
