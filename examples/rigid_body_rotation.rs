extern crate granules;
extern crate ndarray;

use granules::contact_search::LinkedListGrid;
use granules::geometry::get_3d_block;
use granules::physics::rigid_body::RigidBody;
use granules::save_data::{create_output_directory, dump_output};
use ndarray::prelude::*;

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
    x: Array1<f32>,
    y: Array1<f32>,
    z: Array1<f32>,
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
        Array::from_vec(xg),
        Array::from_vec(yg),
        Array::from_vec(zg),
        sim_data.cube_spacing / 2.,
        1000. * sim_data.cube_spacing.powf(2.),
    );
    // calculate all the properties required by rigid body simulation
    cube.pre_compute_values();
    // shift cube up
    cube.x += 2.5;
    cube.y += 1.5;
    cube.z += 1.5;

    let dt = 5e-4;
    let dim = 3;
    let tf = 2.;
    let mut time_step_number = 0;
    let mut t = 0.;
    let scale = 2.;
    println!("Total steps are {}", tf / dt);

    create_output_directory();

    while t < tf {
        let grid = LinkedListGrid::new(&mut vec![&mut cube], scale);
        make_forces_zero(&mut cube);
        body_force_dem(&mut cube, 0., -9.81, 0.);
        spring_force_self(&mut cube, 1e7, &grid, dim);
        spring_force_other(&mut cube, &mut tank, 1e7, &grid, dim);
        integrate(&mut cube, dt);
        t = t + dt;
        if time_step_number % 100 == 0 {
            println!("{:?}", time_step_number);
            dump_output(&mut vec![&mut cube, &mut tank], time_step_number);
        }
        time_step_number += 1;
    }
}
