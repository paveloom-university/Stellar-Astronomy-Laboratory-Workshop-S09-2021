use std::{fs::File, io::Write, path::Path};

use crate::orbit::Orbit;

impl Orbit {
    /// Write the results of integration to the specified file
    pub fn write(&self, folder: &Path) {
        if let Ok(mut f) = File::create(folder.join(format!("{}.dat", self.id))) {
            if writeln!(f, "r z").is_ok() {
                for i in 0..=self.n {
                    writeln!(
                        f,
                        "{:.18} {:.18}",
                        self.gcy_integrated.r[i], self.gcy_integrated.z[i]
                    )
                    .ok();
                }
            }
        }
    }
}
