extern crate granules;

use granules::contact_search::LinkedListGrid;
use granules::geometry::dam_break_3d_geometry;
use granules::integrate::{integrate_initialize, integrate_stage1, integrate_stage2};
use granules::physics::dem::equations::{
    body_force_dem, make_forces_zero, spring_force_other, spring_force_self,
};
use granules::physics::dem::DemDiscrete;
use granules::save_data::{create_output_directory, dump_output};

pub struct SimulationData {
    pub grains_spacing: f32,
    pub grains_length: f32,
    pub grains_height: f32,
    pub grains_width: f32,
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
            grains_spacing: 0.5,
            grains_length: 4.,
            grains_height: 5.,
            grains_width: 5.,
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
    part1: &mut DemDiscrete,
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    h: f32,
    mass: f32,
) {
    let m_inv = 1. / mass;
    for i in 0..part1.len {
        part1.x[i] = x[i];
        part1.y[i] = y[i];
        part1.z[i] = z[i];
        part1.h[i] = h;
        part1.rad[i] = h;
        part1.m[i] = mass;
        part1.m_inv[i] = m_inv;
    }
}

fn main() {
    let sim_data = SimulationData::new();

    let (xg, yg, zg, xt, yt, zt) = dam_break_3d_geometry(
        sim_data.grains_length,
        sim_data.grains_height,
        sim_data.grains_width,
        sim_data.grains_spacing,
        sim_data.tank_length,
        sim_data.tank_height,
        sim_data.tank_width,
        sim_data.tank_spacing,
        sim_data.tank_layers,
        sim_data.layers_outside,
    );

    let mut grains = DemDiscrete::new(xg.len(), 0, "grains".to_string());
    let mut tank = DemDiscrete::new(xt.len(), 1, "tank".to_string());

    setup_particle_properties(
        &mut grains,
        xg,
        yg,
        zg,
        sim_data.grains_spacing / 2.,
        1000. * sim_data.grains_spacing.powf(2.),
    );
    setup_particle_properties(
        &mut tank,
        xt,
        yt,
        zt,
        sim_data.tank_spacing / 2.,
        1000. * sim_data.tank_spacing.powf(2.),
    );
    // shift grains up
    for i in 0..grains.x.len() {
        grains.x[i] += 2.5;
        grains.y[i] += 1.5;
        grains.z[i] += 1.5;
    }

    let dt = 5e-4;
    let dim = 3;
    let tf = 2.;
    let mut time_step_number = 0;
    let mut t = 0.;
    let scale = 2.;
    println!("Total steps are {}", tf / dt);

    create_output_directory();

    while t < tf {
        let grid = LinkedListGrid::new(&mut vec![&mut grains, &mut tank], scale);
        make_forces_zero(&mut grains);
        body_force_dem(&mut grains, 0., -9.81, 0.);
        spring_force_self(&mut grains, 1e7, 0.5, dt, &grid, dim);
        spring_force_other(&mut grains, &mut tank, 1e7, &grid, dim);
        integrate_stage1(&mut vec![&mut grains], dt);
        t = t + dt;
        if time_step_number % 100 == 0 {
            println!("{:?}", time_step_number);
            dump_output(&mut vec![&mut grains, &mut tank], time_step_number);
        }
        time_step_number += 1;
    }
}
