use ndarray::prelude::*;
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

// A mutable struct which is required by all dem equations as a destination
pub struct DEMDiscreteMutDest<'a> {
    pub len: &'a mut usize,
    pub m: &'a mut Array1<f32>,
    pub x: &'a mut Array1<f32>,
    pub y: &'a mut Array1<f32>,
    pub z: &'a mut Array1<f32>,
    pub u: &'a mut Array1<f32>,
    pub v: &'a mut Array1<f32>,
    pub w: &'a mut Array1<f32>,
    pub omega_x: &'a mut Array1<f32>,
    pub omega_y: &'a mut Array1<f32>,
    pub omega_z: &'a mut Array1<f32>,
    pub inertia: &'a mut Array1<f32>,
    pub h: &'a mut Array1<f32>,
    pub m_inv: &'a mut Array1<f32>,
    pub rad: &'a mut Array1<f32>,
    pub fx: &'a mut Array1<f32>,
    pub fy: &'a mut Array1<f32>,
    pub fz: &'a mut Array1<f32>,
    pub tauz: &'a mut Array1<f32>,
    pub taux: &'a mut Array1<f32>,
    pub tauy: &'a mut Array1<f32>,
    pub id: &'a mut usize,
}

pub struct DEMDiscreteMutSource<'a> {
    pub len: &'a mut usize,
    pub m: &'a mut Array1<f32>,
    pub x: &'a mut Array1<f32>,
    pub y: &'a mut Array1<f32>,
    pub z: &'a mut Array1<f32>,
    pub u: &'a mut Array1<f32>,
    pub v: &'a mut Array1<f32>,
    pub w: &'a mut Array1<f32>,
    pub omega_x: &'a mut Array1<f32>,
    pub omega_y: &'a mut Array1<f32>,
    pub omega_z: &'a mut Array1<f32>,
    pub inertia: &'a mut Array1<f32>,
    pub h: &'a mut Array1<f32>,
    pub m_inv: &'a mut Array1<f32>,
    pub rad: &'a mut Array1<f32>,
    pub id: &'a mut usize,
}

pub trait GetMutDEMDest{
    fn get_mut_parts(&mut self) -> DEMDiscreteMutDest;
}

impl GetMutDEMDest for DemDiscrete {
    fn get_mut_parts(&mut self) -> DEMDiscreteMutDest {
        DEMDiscreteMutDest {
            len: &mut self.len,
            m: &mut self.m,
            x: &mut self.x,
            y: &mut self.y,
            z: &mut self.z,
            u: &mut self.u,
            v: &mut self.v,
            w: &mut self.w,
            omega_x: &mut self.omega_x,
            omega_y: &mut self.omega_y,
            omega_z: &mut self.omega_z,
            inertia: &mut self.inertia,
            h: &mut self.h,
            m_inv: &mut self.m_inv,
            rad: &mut self.rad,
            fx: &mut self.fx,
            fy: &mut self.fy,
            fz: &mut self.fz,
            taux: &mut self.taux,
            tauy: &mut self.tauy,
            tauz: &mut self.tauz,
            id: &mut self.id,
        }
    }
}
