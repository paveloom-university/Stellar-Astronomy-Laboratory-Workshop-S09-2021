//! Try to load the file with initial positions and velocities

use csv::{ByteRecord, ReaderBuilder};
use serde::Deserialize;

use crate::orbit::{HCInitials, Orbit, F};

/// A representation of the results after parsing the files
pub struct Log<'a> {
    /// Number of files parsed
    len: usize,
    /// Paths to the parsed files (as specified by a user)
    files: Vec<&'a str>,
    /// A vector of status codes for each file
    statuses: Vec<String>,
}

/// A representation of a file record
#[derive(Deserialize)]
struct Record<'a> {
    /// ID of the object
    id: &'a str,
    /// X component of the radius vector $\[ \text{kpc} \]$
    x_0: F,
    /// Standard deviation of the X component of the radius vector $\[ \text{kpc} \]$
    x_0_err: F,
    /// Y component of the radius vector $\[ \text{kpc} \]$
    y_0: F,
    /// Standard deviation of the Y component of the radius vector $\[ \text{kpc} \]$
    y_0_err: F,
    /// Z component of the radius vector $\[ \text{kpc} \]$
    z_0: F,
    /// Standard deviation of the Z component of the radius vector $\[ \text{kpc} \]$
    z_0_err: F,
    /// U component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    u_0: F,
    /// Standard deviation of the U component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    u_0_err: F,
    /// V component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    v_0: F,
    /// Standard deviation of the V component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    v_0_err: F,
    /// W component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    w_0: F,
    /// Standard deviation of the W component of the velocity vector $\[ \text{km} \\, \text{s}^{-1} \]$
    w_0_err: F,
}

impl Orbit {
    /// Parse the file(s) with initial positions and velocities
    ///
    /// Return all of the proper orbits found in the file(s) and
    /// the log of status codes.
    #[must_use]
    pub fn load(files: Vec<&str>) -> (Vec<Orbit>, Log) {
        // Get the length of files
        let len = files.len();

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
                    != &b"idx_0x_0_erry_0y_0_errz_0z_0_erru_0u_0_errv_0v_0_errw_0w_0_err"[..]
                    || headers.len() != 13
                {
                    // Put the appropriate status and skip this file
                    statuses.push("WH".to_string());
                // Otherwise,
                } else {
                    // Reset the count of records that were successfully deserialized and parsed
                    counter = 0;
                    // While there are records that could be read
                    loop {
                        // Try to read a record
                        if let Ok(read) = rdr.read_byte_record(&mut record) {
                            // If a record was read successfully
                            if read {
                                // If the record could be deserialized, process
                                // the results; otherwise, skip this line
                                if let Ok(record) = record.deserialize::<Record>(Some(&headers)) {
                                    // Use the data of the record to create and save an orbit
                                    orbits.push(Orbit::from(
                                        record.id.to_string(),
                                        HCInitials {
                                            x: record.x_0,
                                            x_err: record.x_0_err,
                                            y: record.y_0,
                                            y_err: record.y_0_err,
                                            z: record.z_0,
                                            z_err: record.z_0_err,
                                            u: record.u_0,
                                            u_err: record.u_0_err,
                                            v: record.v_0,
                                            v_err: record.v_0_err,
                                            w: record.w_0,
                                            w_err: record.w_0_err,
                                        },
                                    ));
                                    counter += 1;
                                }
                            // Otherwise, break since this is the end of file
                            } else {
                                break;
                            }
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
