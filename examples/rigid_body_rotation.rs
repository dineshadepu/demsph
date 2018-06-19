extern crate granules;

use granules::contact_search::LinkedListGrid;
use granules::geometry::get_3d_block;
use granules::physics::rigid_body::RigidBody;
use granules::physics::rigid_body::equations::{body_force_rigid_body, make_forces_zero_rigid_body};
use granules::save_data::{create_output_directory, dump_output};
use granules::integrate::{RK2, integrate_stage1};

pub struct SimulationData {
    pub cube_spacing: f32,
    pub cube_length: f32,
    pub cube_height: f32,
    pub cube_width: f32,
}

impl SimulationData {
    fn new() -> Self {
        SimulationData {
            cube_spacing: 0.5,
            cube_length: 4.,
            cube_height: 5.,
            cube_width: 5.,
        }
    }
}

fn setup_particle_properties(
    part1: &mut RigidBody,
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    h: f32,
    mass: f32,
) {
    for i in 0..part1.len {
        part1.x[i] = x[i];
        part1.y[i] = y[i];
        part1.z[i] = z[i];
        part1.h[i] = h;
        part1.rad[i] = h;
        part1.m[i] = mass;
    }
}

fn main() {
    let sim_data = SimulationData::new();

    let (xg, yg, zg) = get_3d_block(
        sim_data.cube_length,
        sim_data.cube_height,
        sim_data.cube_width,
        sim_data.cube_spacing,
    );

    let mut cube = RigidBody::new(xg.len(), 0, "cube".to_string());

    setup_particle_properties(
        &mut cube,
        xg,
        yg,
        zg,
        sim_data.cube_spacing / 2.,
        1000. * sim_data.cube_spacing.powf(2.),
    );
    // calculate all the properties required by rigid body simulation

    // Assign angular velocity to the body, precomputed function will
    // evaluate its angular momentum
    cube.omega = vec![0., 0., 3.];
    cube.pre_compute_values();

    // shift grains up
    for i in 0..cube.x.len() {
        cube.x[i] += 2.5;
        cube.y[i] += 1.5;
        cube.z[i] += 1.5;
    }

    let dt = 5e-4;
    let dim = 3;
    let tf = 5.;
    let mut time_step_number = 0;
    let mut t = 0.;
    let scale = 2.;

    create_output_directory();

    while t < tf {
        // let grid = LinkedListGrid::new(&mut vec![&mut cube], scale);
        make_forces_zero_rigid_body(&mut vec![&mut cube]);
        // body_force_rigid_body(&mut vec![&mut cube], 0., -9.81, 0.);
        integrate_stage1(&mut vec![&mut cube], dt * 2.);
        t = t + dt;
        if time_step_number % 100 == 0 {
            println!("{:?}", time_step_number);
            dump_output(&mut vec![&mut cube], time_step_number);
        }
        time_step_number += 1;
    }
}
