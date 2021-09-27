//! This binary crate is a command-line utility to integrate the orbits of
//! globular clusters by providing their initial positions and velocities
//! and choosing a model of the Galaxy potential. It is based on the paper
//! of Bajkova and Bobylev ([2020, v1](https://arxiv.org/abs/2008.13624v1)).

use clap::{crate_authors, crate_description, crate_name, crate_version, AppSettings, Arg};
use std::{
    ffi::OsString,
    path::{Path, PathBuf},
};

mod orbit;

use orbit::Orbit;

/// Padding used when printing the output
const PADDING: usize = 5;

/// Run the program
fn main() {
    // Define the CLI
    let matches = clap::app_from_crate!()
        .help_message("Print help information")
        .version_message("Print version information")
        .after_help("Documentation: ???\nReference: ???")
        .args(&[
            Arg::with_name("results")
                .help("Set the result directory to use")
                .required(true)
                .validator_os(|s| {
                    if Path::new(s).is_dir() {
                        Ok(())
                    } else {
                        Err(OsString::from(format!("\nNo such dir {:#?}.", s)))
                    }
                }),
            Arg::with_name("file(s)")
                .help("Set the input file(s) to use")
                .multiple(true)
                .required(true)
                .validator_os(|s| {
                    if Path::new(s).is_file() {
                        Ok(())
                    } else {
                        Err(OsString::from(format!("\nNo such file {:#?}.", s)))
                    }
                }),
        ])
        .settings(&[
            AppSettings::ArgRequiredElseHelp,
            AppSettings::TrailingVarArg,
        ])
        .get_matches();

    // Get path to the results directory
    let results = PathBuf::from(matches.value_of("results").unwrap());
    // Get values of the input files
    let files = matches.values_of("file(s)").unwrap();

    println!("\n{:1$}> Parsing the files...", "", PADDING);

    // Parse the files and get the orbits
    let (orbits, log) = Orbit::load(files);

    // Print the log
    println!("{}", log.format(PADDING + 2));

    // Check how many orbits (initial coordinates and velocities) were parsed
    if orbits.is_empty() {
        println!("{:1$}> No orbits were parsed. Exiting.", "", PADDING);
    } else if orbits.len() == 1 {
        println!("{:1$}> One orbit was parsed.\n", "", PADDING);
    } else {
        println!(
            "{:1$}> {2} orbits were parsed.\n",
            "",
            PADDING,
            orbits.len()
        );
    }

    // Integrate each orbit and write the results
    for ref mut orbit in orbits {
        println!(
            "{:1$}> Integrating the orbit of {2}...",
            "",
            PADDING,
            orbit.id()
        );
        orbit.integrate(5000000, -0.001);
        orbit.write(&results);
    }
}
