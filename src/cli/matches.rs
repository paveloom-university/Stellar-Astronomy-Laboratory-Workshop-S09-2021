//! Define the CLI of the program, return the matched arguments

use clap::{
    crate_authors, crate_description, crate_name, crate_version, AppSettings, Arg, ArgMatches,
};
use std::{ffi::OsString, path::Path};

use crate::{orbit::Results, F, I};

/// Define the `n` argument
fn n() -> Arg<'static, 'static> {
    Arg::with_name("n")
        .help("Set the number of iterations")
        .short("n")
        .long("iterations")
        .takes_value(true)
        .empty_values(false)
        .required(true)
        .validator_os(|s| {
            s.to_str().map_or(
                Err(OsString::from("Argument isn't a valid UTF-8 string.")),
                |sl| {
                    if sl.chars().all(|c| c.is_numeric() || c == '-') {
                        sl.parse::<I>().map_or(
                            Err(OsString::from("Argument isn't a positive integer.")),
                            |i| {
                                if i == 0 {
                                    Err(OsString::from("Argument isn't a positive integer."))
                                } else {
                                    Ok(())
                                }
                            },
                        )
                    } else {
                        Err(OsString::from("Argument isn't numeric."))
                    }
                },
            )
        })
}

/// Define the `h` argument
fn h() -> Arg<'static, 'static> {
    Arg::with_name("h")
        .help("Set the time step (in Myr)")
        .short("h")
        .long("step")
        .takes_value(true)
        .empty_values(false)
        .required(true)
        .validator_os(|s| {
            s.to_str().map_or(
                Err(OsString::from("Argument isn't a valid UTF-8 string.")),
                |sl| {
                    if sl
                        .chars()
                        .all(|c| c.is_numeric() || c == '-' || c == '.' || c == 'e' || c == 'E')
                    {
                        if sl.parse::<F>().is_ok() {
                            Ok(())
                        } else {
                            Err(OsString::from("Argument isn't a float number."))
                        }
                    } else {
                        Err(OsString::from("Argument isn't numeric."))
                    }
                },
            )
        })
}

/// Define the `s` argument
fn s() -> Arg<'static, 'static> {
    Arg::with_name("s")
        .help("[SM] Set the number of Monte Carlo simulations")
        .short("s")
        .long("simulations")
        .takes_value(true)
        .empty_values(false)
        .requires("sim")
        .validator_os(|s| {
            s.to_str().map_or(
                Err(OsString::from("Argument isn't a valid UTF-8 string.")),
                |sl| {
                    if sl.chars().all(|c| c.is_numeric() || c == '-') {
                        sl.parse::<I>().map_or(
                            Err(OsString::from("Argument isn't a positive integer.")),
                            |i| {
                                if i == 0 {
                                    Err(OsString::from("Argument isn't a positive integer."))
                                } else {
                                    Ok(())
                                }
                            },
                        )
                    } else {
                        Err(OsString::from("Argument isn't numeric."))
                    }
                },
            )
        })
}

/// Define the `output` argument
fn output() -> Arg<'static, 'static> {
    Arg::with_name("output")
        .help("Set the output directory to use")
        .short("o")
        .long("output")
        .takes_value(true)
        .empty_values(false)
        .required(true)
        .validator_os(|s| {
            if Path::new(s).is_dir() {
                Ok(())
            } else {
                Err(OsString::from(format!("\nNo such dir {:#?}.", s)))
            }
        })
}

// Define the `fields` argument
fn fields() -> Arg<'static, 'static> {
    Arg::with_name("fields")
        .help("Set the fields to write")
        .short("f")
        .long("fields")
        .takes_value(true)
        .empty_values(false)
        .required(true)
        .require_delimiter(true)
        .possible_values(Results::fields())
        .hide_possible_values(true)
}

/// Define the `file(s)` argument
fn files() -> Arg<'static, 'static> {
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
        })
}

/// Define the `rev` argument
fn rev() -> Arg<'static, 'static> {
    Arg::with_name("rev")
        .help("Enable reverse mode [RM]")
        .long("reverse")
        .conflicts_with("sim")
}

/// Define the `sim` argument
fn sim() -> Arg<'static, 'static> {
    Arg::with_name("sim")
        .help("Enable simulation mode [SM]")
        .long("simulate")
        .requires("s")
}

/// Define the CLI of the program, return the matched arguments
pub fn get() -> ArgMatches<'static> {
    clap::app_from_crate!()
        .help_message("Print help information")
        .version_message("Print version information")
        .after_help("Documentation: ???\nReference: ???")
        .args(&[n(), h(), s(), output(), fields(), files(), rev(), sim()])
        .settings(&[
            AppSettings::AllowNegativeNumbers,
            AppSettings::ArgRequiredElseHelp,
            AppSettings::TrailingVarArg,
        ])
        .get_matches()
}
