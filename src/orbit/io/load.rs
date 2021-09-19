//! Try to load the file with initial positions and velocities

use clap::Values;
use csv::{ByteRecord, ReaderBuilder};

use super::super::{Initials, Orbit};

/// Status codes for the parsing process
enum Status {
    /// A CSV reader couldn't be built
    ReaderBuilderFailed,
    /// Couldn't parse the headers
    ParsingHeadersFailed,
    /// The headers are wrong
    WrongHeader,
    /// Parsed successfully this amount of orbits
    Parsed(i32),
}

/// A representation of the results after parsing the files
pub struct Log<'a> {
    /// Number of files parsed
    len: usize,
    /// Paths to the parsed files (as specified by a user)
    files: Vec<&'a str>,
    /// A vector of status codes for each file
    statuses: Vec<Status>,
}

impl Orbit {
    /// Parse the file(s) with initial positions and velocities
    ///
    /// Return all of the proper orbits found in the file(s) and
    /// the log of status codes.
    pub fn load(files: Values) -> (Vec<Orbit>, Log) {
        // Get the length of files
        let len = files.len();
        let files: Vec<&str> = files.collect();

        // Prepare an array of orbits
        let mut orbits = Vec::<Orbit>::new();
        // Prepare a mask of status codes
        let mut statuses = Vec::<Status>::with_capacity(len);
        // Prepare two byte records for headers and data
        let mut headers = ByteRecord::new();
        let mut record = ByteRecord::new();
        // Prepare a counter of records that were successfully deserialized and parsed
        let mut counter;

        // For each file in the files
        for file in &files {
            // If it's possible to build a CSV reader off the file
            if let Ok(mut rdr) = ReaderBuilder::new().delimiter(b' ').from_path(file) {
                // Parse the first row of the file (this operation borrows mutably)
                if let Ok(hs) = rdr.byte_headers() {
                    // Clone the headers, end the mutable borrow
                    headers = hs.clone();
                };
                // If the headers are empty
                if headers.is_empty() {
                    // Put the appropriate status and skip this file
                    statuses.push(Status::ParsingHeadersFailed);
                // Or, if the header isn't correct
                } else if headers.as_slice()
                    != &b"x_0x_0_erry_0y_0_errz_0z_0_erru_0u_0_errv_0v_0_errw_0w_0_err"[..]
                    || headers.len() != 12
                {
                    // Put the appropriate status and skip this file
                    statuses.push(Status::WrongHeader);
                // Otherwise,
                } else {
                    // Reset the count of records that were successfully deserialized and parsed
                    counter = 0;
                    // While there are records that could be read
                    while let Ok(read) = rdr.read_byte_record(&mut record) {
                        // If a record was read
                        if read {
                            // If the record could be deserialized
                            if let Ok(initials) = record.deserialize::<Initials>(Some(&headers)) {
                                // Create and save an orbit from initial coordinates and velocities
                                orbits.push(Orbit::from(initials));
                                counter += 1;
                            }
                        // Otherwise,
                        } else {
                            break;
                        }
                    }
                    // Put the appropriate status and finish with this file
                    statuses.push(Status::Parsed(counter));
                }
            // Otherwise,
            } else {
                // Put the appropriate status and skip this file
                statuses.push(Status::ReaderBuilderFailed);
            }
        }

        (
            orbits,
            Log {
                len,
                files,
                statuses,
            },
        )
    }
}

impl<'a> Log<'a> {
    /// Format the contents of the log with the specified padding, return a string
    pub fn format(&self, padding: usize) -> String {
        // Prepare strings for the output and the status codes
        let mut s = String::new();
        let mut code = String::new();
        // For each file
        for i in 0..self.len {
            // Get its status
            if let Some(status) = self.statuses.get(i) {
                // Define the string representation of the status code
                match status {
                    Status::ReaderBuilderFailed => code = "RF".to_owned(),
                    Status::ParsingHeadersFailed => code = "PF".to_owned(),
                    Status::WrongHeader => code = "WH".to_owned(),
                    Status::Parsed(i) => code = format!("{}", i),
                }
            }
            // Append a string with the stringified status code and the file path
            s.push_str(&format!(
                "{:1$}({2}) {3:?}\n",
                "", padding, &code, self.files[i]
            ));
        }
        s
    }
}
