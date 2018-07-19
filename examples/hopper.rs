#[macro_use]
extern crate demsph;

use demsph::contact_search::LinkedListGrid;
use demsph::geometry::{grid_2d, hopper_2d};
use demsph::integrate::{integrate_initialize, integrate_stage1, integrate_stage2};
use demsph::physics::dem::DemDiscrete;
use demsph::physics::dem::equations::{body_force_dem, linear_viscoelastic_model_dem_other,
                                      linear_viscoelastic_model_dem_self, make_forces_zero};
use demsph::save_data::{create_output_directory, dump_output};

pub struct SimulationData {
    pub grains_spacing: f32,
    pub grains_length: f32,
    pub grains_height: f32,
    pub hopper_spacing: f32,
    pub hopper_br: f32,
    pub hopper_tr: f32,
    pub hopper_height: f32,
}

impl SimulationData {
    fn new() -> Self {
        SimulationData {
            grains_spacing: 0.3,
            grains_length: 4.,
            grains_height: 5.,
            hopper_spacing: 0.3,
            hopper_tr: 5.,
            hopper_br: 1.,
            hopper_height: 7.,
        }
    }
}

fn setup_particle_properties(part1: &mut DemDiscrete, x: Vec<f32>, y: Vec<f32>, h: f32, mass: f32) {
    let m_inv = 1. / mass;
    for i in 0..part1.len {
        part1.x[i] = x[i];
        part1.y[i] = y[i];
        part1.h[i] = h;
        part1.rad[i] = h;
        part1.m[i] = mass;
        part1.m_inv[i] = m_inv;
    }
}

fn main() {
    let sim_data = SimulationData::new();

    let (xh, yh) = hopper_2d(
        sim_data.hopper_br,
        sim_data.hopper_tr,
        sim_data.hopper_height,
        sim_data.hopper_spacing,
    );

    let (xg, yg) = grid_2d(
        sim_data.grains_length,
        sim_data.grains_height,
        sim_data.grains_spacing,
    );
    let mut grains = DemDiscrete::new(xg.len(), 0, "grains".to_string());
    let mut hopper = DemDiscrete::new(xh.len(), 1, "hopper".to_string());

    setup_particle_properties(
        &mut grains,
        xg,
        yg,
        sim_data.grains_spacing / 2.,
        1000. * sim_data.grains_spacing.powf(2.),
    );
    setup_particle_properties(
        &mut hopper,
        xh,
        yh,
        sim_data.hopper_spacing / 2.,
        1000. * sim_data.hopper_spacing.powf(2.),
    );

    // move the grains left
    for i in 0..grains.len{
        grains.x[i] -= 2.;
        grains.y[i] += 3.;
    }

    let dt = 1e-4;
    let dim = 2;
    let tf = 2.;
    let mu = 0.4;
    let mut time_step_number = 0;
    let mut t = 0.;
    let scale = 2.;
    let stage1 = 1;
    let stage2 = 2;

    let dir_name = create_directory_return_name![];
    let pfreq = 100;

    while t < tf {
        let grid = LinkedListGrid::new(&mut vec![&mut grains, &mut hopper], scale);
        // initialize the components
        integrate_initialize(&mut vec![&mut grains], dt);

        // Compute the forces
        make_forces_zero(&mut grains);
        body_force_dem(&mut grains, 0., -9.81, 0.0);
        linear_viscoelastic_model_dem_other(
            &mut grains,
            &mut hopper,
            1e7,
            mu,
            dt,
            stage1,
            &grid,
            dim,
        );
        linear_viscoelastic_model_dem_self(&mut grains, 1e7, 0.0, dt, stage1, &grid, dim);

        // Execulte stage 1
        integrate_stage1(&mut vec![&mut grains], dt);

        // Compute the forces
        make_forces_zero(&mut grains);
        body_force_dem(&mut grains, 0., -9.81, 0.0);
        linear_viscoelastic_model_dem_other(
            &mut grains,
            &mut hopper,
            1e7,
            mu,
            dt,
            stage2,
            &grid,
            dim,
        );
        linear_viscoelastic_model_dem_self(&mut grains, 1e7, 0.0, dt, stage2, &grid, dim);

        // Execulte stage 2
        integrate_stage2(&mut vec![&mut grains], dt);

        // increase the time
        t = t + dt;
        if time_step_number % pfreq == 0 {
            println!("{:?}", time_step_number);
            dump_output(
                &mut vec![&mut grains, &mut hopper],
                time_step_number,
                &dir_name,
            );
        }
        time_step_number += 1;
    }
}
