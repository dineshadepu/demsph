use ndarray::prelude::*;
use contact_search::{NNPS, NNPSMutParts};
pub mod equations;

pub struct DemDiscrete {
    pub len: usize,
    pub m: Array1<f32>,
    pub x: Array1<f32>,
    pub y: Array1<f32>,
    pub z: Array1<f32>,
    pub u: Array1<f32>,
    pub v: Array1<f32>,
    pub w: Array1<f32>,
    pub omega_x: Array1<f32>,
    pub omega_y: Array1<f32>,
    pub omega_z: Array1<f32>,
    pub inertia: Array1<f32>,
    pub h: Array1<f32>,
    pub m_inv: Array1<f32>,
    pub rad: Array1<f32>,
    pub fx: Array1<f32>,
    pub fy: Array1<f32>,
    pub fz: Array1<f32>,
    pub taux: Array1<f32>,
    pub tauy: Array1<f32>,
    pub tauz: Array1<f32>,
    pub id: usize,
    pub name: String,
}

impl DemDiscrete {
    pub fn new(len: usize, id: usize, name: String) -> Self {
        DemDiscrete {
            len,
            name,
            id,
            m: Array1::zeros(len),
            x: Array1::zeros(len),
            y: Array1::zeros(len),
            z: Array1::zeros(len),
            u: Array1::zeros(len),
            v: Array1::zeros(len),
            w: Array1::zeros(len),
            omega_x: Array1::zeros(len),
            omega_y: Array1::zeros(len),
            omega_z: Array1::zeros(len),
            inertia: Array1::zeros(len),
            h: Array1::zeros(len),
            m_inv: Array1::zeros(len),
            rad: Array1::zeros(len),
            fx: Array1::zeros(len),
            fy: Array1::zeros(len),
            fz: Array1::zeros(len),
            tauz: Array1::zeros(len),
            taux: Array1::zeros(len),
            tauy: Array1::zeros(len),
        }
    }
}


impl NNPS for DemDiscrete{
    fn get_parts_mut(&mut self) -> NNPSMutParts {
        NNPSMutParts {
            len: &mut self.len,
            x: &mut self.x,
            y: &mut self.y,
            z: &mut self.z,
            h: &mut self.h,
            id: &mut self.id,
        }
    }

    fn get_x(&self) -> &Array1<f32> {
        &self.x
    }
    fn get_y(&self) -> &Array1<f32> {
        &self.y
    }
    fn get_z(&self) -> &Array1<f32> {
        &self.z
    }
}
