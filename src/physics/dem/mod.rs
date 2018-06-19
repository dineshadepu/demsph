use contact_search::{NNPS, NNPSMutParts};
pub mod equations;

pub struct DemDiscrete {
    pub len: usize,
    pub m: Vec<f32>,
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Vec<f32>,
    pub u: Vec<f32>,
    pub v: Vec<f32>,
    pub w: Vec<f32>,
    pub omega_x: Vec<f32>,
    pub omega_y: Vec<f32>,
    pub omega_z: Vec<f32>,
    pub inertia: Vec<f32>,
    pub h: Vec<f32>,
    pub m_inv: Vec<f32>,
    pub rad: Vec<f32>,
    pub fx: Vec<f32>,
    pub fy: Vec<f32>,
    pub fz: Vec<f32>,
    pub taux: Vec<f32>,
    pub tauy: Vec<f32>,
    pub tauz: Vec<f32>,
    pub id: usize,
    pub name: String,
}

impl DemDiscrete {
    pub fn new(len: usize, id: usize, name: String) -> Self {
        DemDiscrete {
            len,
            name,
            id,
            m: vec![0.; len],
            x: vec![0.; len],
            y: vec![0.; len],
            z: vec![0.; len],
            u: vec![0.; len],
            v: vec![0.; len],
            w: vec![0.; len],
            omega_x: vec![0.; len],
            omega_y: vec![0.; len],
            omega_z: vec![0.; len],
            inertia: vec![0.; len],
            h: vec![0.; len],
            m_inv: vec![0.; len],
            rad: vec![0.; len],
            fx: vec![0.; len],
            fy: vec![0.; len],
            fz: vec![0.; len],
            tauz: vec![0.; len],
            taux: vec![0.; len],
            tauy: vec![0.; len],
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

    fn get_x(&self) -> &Vec<f32> {
        &self.x
    }
    fn get_y(&self) -> &Vec<f32> {
        &self.y
    }
    fn get_z(&self) -> &Vec<f32> {
        &self.z
    }
}
