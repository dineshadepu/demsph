use super::physics::dem::DemDiscrete;
use super::physics::rigid_body::RigidBody;
use std;
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

pub fn get_output_directory_name() -> String {
    let name = std::env::current_exe().unwrap();
    let name = Path::new(&name)
        .file_name()
        .unwrap()
        .to_os_string()
        .into_string()
        .unwrap();
    let crate_dir = env!("CARGO_MANIFEST_DIR");
    let dir_name = format!("{}/examples/{}_output", crate_dir, name,);
    dir_name.to_string()
}

pub fn create_output_directory() {
    let dir_name = get_output_directory_name();
    let _ = fs::create_dir_all(dir_name).expect("Couldn't create directory");
}

pub trait DumpData {
    fn save_data(&self, &str, usize);
    fn write_vtk_xml(&self, &str, usize);
}

impl DumpData for RigidBody {
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
            writeln!(&mut file, "{} {} {}", self.x[i], self.y[i], self.z[i]).unwrap();
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
            writeln!(&mut file, "{} {} {}", self.u[i], self.v[i], self.w[i]).unwrap();
        }
        // write forces
        writeln!(&mut file, "VECTORS Force float").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{} {} {}", self.fx[i], self.fy[i], self.fz[i]).unwrap();
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
            writeln!(&mut file, "{} {} {}", self.x[i], self.y[i], self.z[i]).unwrap();
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
            writeln!(&mut file, "{} {} {}", self.u[i], self.v[i], self.w[i]).unwrap();
        }
        // write forces
        writeln!(&mut file, "VECTORS Force float").unwrap();
        for i in 0..np {
            writeln!(&mut file, "{} {} {}", self.fx[i], self.fy[i], self.fz[i]).unwrap();
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

pub fn dump_output<T: DumpData>(entities: &mut Vec<&mut T>, time_step_number: usize) {
    let output_folder = get_output_directory_name();
    for entity in entities {
        entity.save_data(&output_folder, time_step_number);
    }
}
