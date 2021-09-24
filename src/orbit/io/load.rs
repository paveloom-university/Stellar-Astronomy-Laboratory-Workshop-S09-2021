//! Try to load the file with initial positions and velocities

use clap::Values;
use csv::{ByteRecord, ReaderBuilder};

use crate::orbit::{HeliocentricCartesianInitials, Orbit};

/// A representation of the results after parsing the files
pub struct Log<'a> {
    /// Number of files parsed
    len: usize,
    /// Paths to the parsed files (as specified by a user)
    files: Vec<&'a str>,
    /// A vector of status codes for each file
    statuses: Vec<String>,
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
        let mut statuses = Vec::<String>::with_capacity(len);
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
                    statuses.push("PF".to_string());
                // Or, if the header isn't correct
                } else if headers.as_slice()
                    != &b"x_0x_0_erry_0y_0_errz_0z_0_erru_0u_0_errv_0v_0_errw_0w_0_err"[..]
                    || headers.len() != 12
                {
                    // Put the appropriate status and skip this file
                    statuses.push("WH".to_string());
                // Otherwise,
                } else {
                    // Reset the count of records that were successfully deserialized and parsed
                    counter = 0;
                    // While there are records that could be read
                    while let Ok(read) = rdr.read_byte_record(&mut record) {
                        // If a record was read
                        if read {
                            // If the record could be deserialized
                            if let Ok(hc_initials) =
                                record.deserialize::<HeliocentricCartesianInitials>(Some(&headers))
                            {
                                // Create and save an orbit from initial coordinates and velocities
                                orbits.push(Orbit::from(hc_initials));
                                counter += 1;
                            }
                        // Otherwise,
                        } else {
                            break;
                        }
                    }
                    // Put the appropriate status and finish with this file
                    statuses.push(counter.to_string());
                }
            // Otherwise,
            } else {
                // Put the appropriate status and skip this file
                statuses.push("RF".to_string());
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
        // For each file
        for i in 0..self.len {
            // Append a string with the stringified status code and the file path
            s.push_str(&format!(
                "{:1$}({2}) {3:?}\n",
                "", padding, &self.statuses[i], self.files[i]
            ));
        }
        s
    }
}
