//! Write the results of integration to the specified file

use byteorder::{LittleEndian, WriteBytesExt};
use std::{
    fs::{self, OpenOptions},
    path::Path,
};

use crate::orbit::Orbit;

impl Orbit {
    /// Write the results of integration to the specified file
    pub fn write(&self, folder: &Path, fields: &[&str], append: bool) {
        // Define the path to the object directory
        let object_dir = folder.join(&self.id);

        // Prepare a result buffer
        let mut r: Result<_, _> = Ok(());

        // If the object directory doesn't exist, create it
        if !object_dir.exists() {
            r = fs::create_dir(&object_dir);
        }
        // If the object directory is ready
        if r.is_ok() {
            /// Do a write call for this field
            macro_rules! write_field {
                ($field:ident) => {
                    if fields.contains(&stringify!($field)) {
                        if let Ok(mut f) = OpenOptions::new()
                            .write(true)
                            .truncate(!append)
                            .append(append)
                            .create(true)
                            .open(object_dir.join(format!("{}.bin", stringify!($field))))
                        {
                            for v in &self.results.$field {
                                f.write_f64::<LittleEndian>(*v).unwrap()
                            }
                        }
                    }
                };
            }

            // Write the results

            // Solutions of the Lagrangian equations
            write_field!(r);
            write_field!(psi);
            write_field!(z);
            write_field!(p_r);
            write_field!(p_psi);
            write_field!(p_z);

            // Galactic Cartesian system
            write_field!(x);
            write_field!(y);

            // Total energy
            write_field!(e);

            // Apocentric and pericentric distances
            write_field!(apo);
            write_field!(peri);
        }
    }
}
