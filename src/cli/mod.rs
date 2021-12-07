//! This module provides the CLI routines
//!
//! The interface of the program can be obtained by
//! calling the program with the `--help` argument:
//!
//! ```bash
//! kidou 0.1.0
//! Pavel Sobolev <paveloom@riseup.net>
//! Calculate the orbits of globular clusters in a model of the Galaxy potential
//!
//! USAGE:
//! kidou [FLAGS] [OPTIONS] <file(s)>... --fields <fields> --step <h> --model <model> --iterations <n> --output <output>
//!
//! FLAGS:
//!         --help        Print help information
//!         --reverse     Enable reverse mode [RM]
//!         --simulate    Enable simulation mode [SM]
//!     -V, --version     Print version information
//!
//! OPTIONS:
//!     -f, --fields <fields>    Set the fields to write
//!     -h, --step <h>           Set the time step (in Myr)
//!         --model <model>      Set the model by index
//!     -n, --iterations <n>     Set the number of iterations
//!     -o, --output <output>    Set the output directory to use
//!     -s, --simulations <s>    [SM] Set the number of Monte Carlo simulations

//! ARGS:
//!     <file(s)>...    Set the input file(s) to use
//!
//! Documentation:
//! - https://codeberg.org/paveloom-university/Stellar-Astronomy-Laboratory-Workshop-S09-2021/src/branch/main/README.md
//! - https://github.com/paveloom-university/Stellar-Astronomy-Laboratory-Workshop-S09-2021/blob/main/README.md
//! - https://gitlab.com/paveloom-g/university/s09-2021/stellar-astronomy-laboratory-workshop/-/blob/main/README.md
//! - https://git.sr.ht/~paveloom/Stellar-Astronomy-Laboratory-Workshop-S09-2021/tree/main/item/README.md
//!
//! Reference:
//! - https://paveloom-university.github.io/Stellar-Astronomy-Laboratory-Workshop-S09-2021
//! - https://paveloom-g.gitlab.io/university/s09-2021/stellar-astronomy-laboratory-workshop
//! ```
//!
//! ## Examples
//!
//! 1. Integrate the orbits for 5 Gyr backward using the [first model](crate::orbit::calc::models#group-1):
//!
//! ```bash
//! cargo run --release -- --model 1 -n 500000 -h -0.01 -f=x,y,r,z -o data/output/M1 data/input/initial.dat
//! ```
//!
//! 2. Integrate the orbits for 5 Gyr backward and 5 Gyr forward (reverse mode):
//!
//! ```bash
//! cargo run --release -- --reverse --model 1 -n 500000 -h -0.01 -f=x,y,r,z -o data/output/M1R data/input/initial.dat
//! ```
//!
//! 3. Run 200 Monte Carlo simulations for 1 Gyr backward using the [second model](crate::orbit::calc::models#group-2):
//!
//! ```bash
//! cargo run --release -- --model 2 -n 100000 -h -0.01 -s 200 --simulate -f=r,z,x,y,apo,peri -o data/output/M2 data/input/initial.dat
//! ```

mod matches;

pub use matches::get as get_matches;
