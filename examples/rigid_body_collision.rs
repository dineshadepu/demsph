extern crate granules;

use granules::contact_search::LinkedListGrid;
use granules::geometry::dam_break_3d_geometry;
use granules::integrate::integrate_stage1;
use granules::physics::rigid_body::equations::{
    aggregate_forces_moments, body_force_rigid_body, make_forces_zero_rigid_body,
    spring_force_rigid_other,
};
use granules::physics::rigid_body::RigidBody;
use granules::save_data::{create_output_directory, dump_output};

pub struct SimulationData {
    pub body_spacing: f32,
    pub body_length: f32,
    pub body_height: f32,
    pub body_width: f32,
    pub tank_spacing: f32,
    pub tank_length: f32,
    pub tank_height: f32,
    pub tank_width: f32,
    pub tank_layers: usize,
    pub layers_outside: bool,
}

impl SimulationData {
    fn new() -> Self {
        SimulationData {
            body_spacing: 0.2,
            body_length: 2.,
            body_height: 2.,
            body_width: 2.,
            tank_spacing: 0.5,
            tank_length: 10.,
            tank_height: 7.,
            tank_width: 7.,
            tank_layers: 1,
            layers_outside: true,
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

    let (xg, yg, zg, xt, yt, zt) = dam_break_3d_geometry(
        sim_data.body_length,
        sim_data.body_height,
        sim_data.body_width,
        sim_data.body_spacing,
        sim_data.tank_length,
        sim_data.tank_height,
        sim_data.tank_width,
        sim_data.tank_spacing,
        sim_data.tank_layers,
        sim_data.layers_outside,
    );

    let mut body = RigidBody::new(xg.len(), 0, "body".to_string());
    let mut tank = RigidBody::new(xt.len(), 1, "tank".to_string());

    setup_particle_properties(
        &mut body,
        xg,
        yg,
        zg,
        sim_data.body_spacing / 2.,
        1000. * sim_data.body_spacing.powf(2.),
    );
    setup_particle_properties(
        &mut tank,
        xt,
        yt,
        zt,
        sim_data.tank_spacing / 2.,
        1000. * sim_data.tank_spacing.powf(2.),
    );

    for i in 0..body.x.len() {
        body.x[i] += 3.5;
        body.y[i] += 2.5;
        body.z[i] += 1.5;
    }
    body.omega = vec![0., 0., 0.];
    body.pre_compute_values();

    // shift body up

    let dt = 1e-3;
    let dim = 3;
    let tf = 5.;
    let mut time_step_number = 0;
    let mut t = 0.;
    let scale = 2.;
    println!("Total steps are {}", tf / dt);

    create_output_directory();

    while t < tf {
        let grid = LinkedListGrid::new(&mut vec![&mut body, &mut tank], scale);
        make_forces_zero_rigid_body(&mut vec![&mut body]);
        body_force_rigid_body(&mut vec![&mut body], 0., -9.81, 0.);
        spring_force_rigid_other(&mut body, &mut tank, 1e8, &grid, dim);
        aggregate_forces_moments(&mut vec![&mut body]);
        integrate_stage1(&mut vec![&mut body], dt * 2.);
        t = t + dt;
        if time_step_number % 100 == 0 {
            println!("{:?}", time_step_number);
            dump_output(&mut vec![&mut body, &mut tank], time_step_number);
        }
        time_step_number += 1;
    }
}
