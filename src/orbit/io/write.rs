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
            // Write the `r` vector
            if let Ok(mut f_r) = File::create(object_dir.join("r.bin")) {
                for v in &self.integrated.r {
                    f_r.write_f64::<LittleEndian>(*v).ok();
                }
            }
            // Write the `z` vector
            if let Ok(mut f_z) = File::create(object_dir.join("z.bin")) {
                for v in &self.integrated.z {
                    f_z.write_f64::<LittleEndian>(*v).ok();
                }
            }
            // Write the `e` vector
            if let Ok(mut f_e) = File::create(object_dir.join("e.bin")) {
                for v in &self.integrated.e {
                    f_e.write_f64::<LittleEndian>(*v).ok();
                }
            }
        }
    }
}
