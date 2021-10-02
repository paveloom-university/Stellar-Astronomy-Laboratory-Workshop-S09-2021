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
    // Get values of the input files
    let files = matches.values_of("file(s)").unwrap();
    // Get switch value of the reverse mode
    let rev = matches.occurrences_of("rev") == 1;

    println!("\n{:1$}> Parsing the files...", "", PADDING);

    // Parse the files and get the orbits
    let (orbits, log) = Orbit::load(files);
    // Print the parsing log
    println!("{}", log.format(PADDING + 2));

    // Check how many orbits (initial coordinates and velocities) were parsed
    if orbits.is_empty() {
        println!("{:1$}> No orbits were parsed. Exiting.", "", PADDING);
    } else if orbits.len() == 1 {
        println!("{:1$}> One orbit was parsed.\n", "", PADDING);
        println!("{:1$}h = {2}", "", PADDING + 2, h);
        println!("{:1$}n = {2}\n", "", PADDING + 2, n);
        if rev {
            println!("{:1$}Reverse mode enabled.\n", "", PADDING + 2);
        }
    } else {
        println!(
            "{:1$}> {2} orbits were parsed.\n",
            "",
            PADDING,
            orbits.len()
        );
        println!("{:1$}h = {2}", "", PADDING + 2, h);
        println!("{:1$}n = {2}\n", "", PADDING + 2, n);
        if rev {
            println!("{:1$}Reverse mode enabled.\n", "", PADDING + 2);
        }
    }

    // Integrate each orbit and write the results
    for ref mut orbit in orbits {
        println!(
            "{:1$}> Integrating the orbit of {2}...",
            "",
            PADDING,
            orbit.id()
        );
        orbit.integrate(n, h, rev);
        orbit.write(&output);
    }

    println!();
}
