use byteorder::{LittleEndian, WriteBytesExt};
use std::{
    fs::{self, File},
    path::Path,
};

use crate::orbit::Orbit;

impl Orbit {
    /// Write the results of integration to the specified file
    pub fn write(&self, folder: &Path) {
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
            // Define the macro for easier write calls
            macro_rules! write_vector {
                ($field:ident) => {
                    if let Ok(mut f) =
                        File::create(object_dir.join(format!("{}.bin", stringify!($field))))
                    {
                        for v in &self.integrated.$field {
                            f.write_f64::<LittleEndian>(*v).ok();
                        }
                    }
                };
            }

            // Write the results

            // Galactic Cylindrical system
            write_vector!(r);
            write_vector!(e);
            write_vector!(z);

            // Galactic Cartesian system
            write_vector!(x);
            write_vector!(y);

            // Total energy
            write_vector!(e);
            write_vector!(vel);
            write_vector!(phi);
        }
    }
}
