use super::physics::dem::DemDiscrete;
use std::fs;
use std::fs::File;
use std::io::prelude::*;

#[macro_export]
macro_rules! create_directory_return_name {
    () => {{
        let mut dir = "".to_string();
        dir.push_str(&env!("CARGO_MANIFEST_DIR").to_string());
        dir.push_str("/");
        let file_name = file!().to_string();
        dir.push_str(&file_name[0..file_name.len() - 3]);
        dir.push_str("_output");
        create_output_directory(&dir);

        dir
    }};
}

pub fn get_output_directory_name(crate_root: String, file_name: String) -> String {
    let mut dir = "".to_string();
    dir.push_str(&crate_root);
    dir.push_str("/");
    dir.push_str(&file_name[0..file_name.len() - 3]);
    dir.push_str("_output");
    dir
}

pub fn create_output_directory(dir_name: &str) {
    let _ = fs::create_dir_all(dir_name).expect("Couldn't create directory");
}

pub trait DumpData {
    fn save_data(&self, &str, usize);
    fn write_vtk_xml(&self, &str, usize);
}

impl DumpData for DemDiscrete {
    fn save_data(&self, output_folder_name: &str, time_step_number: usize) {
        let file_name = format!(
            "{}/{}_{}.vtk",
            output_folder_name, self.name, time_step_number
        );

        // create the file
        let mut file = File::create(file_name).expect("Could not create file!");
        writeln!(&mut file, "# vtk DataFile Version 2.0").unwrap();
        writeln!(&mut file, "Data values of grains").unwrap();
        writeln!(&mut file, "ASCII").unwrap();
        writeln!(&mut file, "").unwrap();
        writeln!(&mut file, "DATASET POLYDATA").unwrap();

        // write header of positions of the entity
        let np = self.x.len();
        writeln!(&mut file, "POINTS {} float", np).unwrap();

        // write the positions
        for i in 0..np {
            writeln!(&mut file, "{} {} {}", self.x[i], self.y[i], 0.).unwrap();
        }
        // write diameter
        writeln!(&mut file, "Point_DATA {}", np).unwrap();
        writeln!(&mut file, "SCALARS Diameter float 1").unwrap();
        writeln!(&mut file, "LOOKUP_TABLE default").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{}", 2. * self.rad[i]).unwrap();
        }
        // write mass
        writeln!(&mut file, "SCALARS Mass float 1").unwrap();
        writeln!(&mut file, "LOOKUP_TABLE default").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{}", self.m[i]).unwrap();
        }
        // write velocity
        writeln!(&mut file, "VECTORS Velocity float").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{} {} {}", self.u[i], self.v[i], 0.).unwrap();
        }
        // write forces
        writeln!(&mut file, "VECTORS Force float").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{} {} {}", self.fx[i], self.fy[i], 0.).unwrap();
        }
    }

    fn write_vtk_xml(&self, output_folder_name: &str, time_step_number: usize) {
        let file_name = format!(
            "{}/{}_{}.vtu",
            output_folder_name, self.name, time_step_number
        );

        // create the file
        let mut file = File::create(file_name).expect("Could not create file!");
        writeln!(&mut file, "# vtk DataFile Version 2.0").unwrap();
        writeln!(&mut file, "Data values of grains").unwrap();
        writeln!(&mut file, "ASCII").unwrap();
        writeln!(&mut file, "").unwrap();
        writeln!(&mut file, "DATASET POLYDATA").unwrap();

        // write header of positions of the entity
        let np = self.x.len();
        writeln!(&mut file, "POINTS {} float", np).unwrap();

        // write the positions
        for i in 0..np {
            writeln!(&mut file, "{} {} 0", self.x[i], self.y[i]).unwrap();
        }
        // write diameter
        writeln!(&mut file, "Point_DATA {}", np).unwrap();
        writeln!(&mut file, "SCALARS Diameter float 1").unwrap();
        writeln!(&mut file, "LOOKUP_TABLE default").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{}", 2. * self.rad[i]).unwrap();
        }
        // write mass
        writeln!(&mut file, "SCALARS Mass float 1").unwrap();
        writeln!(&mut file, "LOOKUP_TABLE default").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{}", self.m[i]).unwrap();
        }
        // write velocity
        writeln!(&mut file, "VECTORS Velocity float").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{} {} 0.0", self.u[i], self.v[i]).unwrap();
        }
        // write forces
        writeln!(&mut file, "VECTORS Force float").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{} {} 0.0", self.fx[i], self.fy[i]).unwrap();
        }
    }
}

pub fn dump_output<T: DumpData>(
    entities: &mut Vec<&mut T>,
    time_step_number: usize,
    dir_name: &str,
) {
    for entity in entities {
        entity.save_data(&dir_name, time_step_number);
    }
}
