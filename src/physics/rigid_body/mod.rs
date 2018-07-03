use contact_search::{NNPSMutParts, NNPS};
use rulinalg::matrix::Matrix;
pub mod equations;

pub struct RigidBody {
    pub len: usize,
    pub m: Vec<f32>,
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Vec<f32>,
    pub x_body: Vec<f32>,
    pub y_body: Vec<f32>,
    pub z_body: Vec<f32>,
    pub x_body_global: Vec<f32>,
    pub y_body_global: Vec<f32>,
    pub z_body_global: Vec<f32>,
    pub u: Vec<f32>,
    pub v: Vec<f32>,
    pub w: Vec<f32>,
    pub fx: Vec<f32>,
    pub fy: Vec<f32>,
    pub fz: Vec<f32>,
    pub h: Vec<f32>,
    pub rad: Vec<f32>,
    // properties of centre of mass
    pub r_com: Vec<f32>,
    pub v_com: Vec<f32>,
    pub omega: Vec<f32>,
    pub omega_matrix: Matrix<f32>,
    pub ang_mom: Vec<f32>,
    // this wont change
    pub mass_matrix_body: Matrix<f32>,
    pub mass_matrix_body_inverse: Matrix<f32>,
    pub mass_matrix_global: Matrix<f32>,
    pub mass_matrix_global_inverse: Matrix<f32>,
    pub orientation: Matrix<f32>,
    pub orientation_transpose: Matrix<f32>,
    pub m_total: f32,
    pub frc: Vec<f32>,
    pub tau: Vec<f32>,
    pub id: usize,
    pub name: String,
}

impl RigidBody {
    pub fn new(len: usize, id: usize, name: String) -> Self {
        RigidBody {
            len,
            name,
            id,
            m: vec![0.; len],
            x: vec![0.; len],
            y: vec![0.; len],
            z: vec![0.; len],
            x_body: vec![0.; len],
            y_body: vec![0.; len],
            z_body: vec![0.; len],
            x_body_global: vec![0.; len],
            y_body_global: vec![0.; len],
            z_body_global: vec![0.; len],
            u: vec![0.; len],
            v: vec![0.; len],
            w: vec![0.; len],
            fx: vec![0.; len],
            fy: vec![0.; len],
            fz: vec![0.; len],
            h: vec![0.; len],
            rad: vec![0.; len],
            r_com: vec![0., 0., 0.],
            v_com: vec![0., 0., 0.],
            omega: vec![0., 0., 0.],
            ang_mom: vec![0., 0., 0.],
            omega_matrix: Matrix::<f32>::zeros(3, 3),
            orientation: Matrix::<f32>::zeros(3, 3),
            orientation_transpose: Matrix::<f32>::zeros(3, 3),
            mass_matrix_body: Matrix::<f32>::zeros(3, 3),
            mass_matrix_body_inverse: Matrix::<f32>::zeros(3, 3),
            mass_matrix_global: Matrix::<f32>::zeros(3, 3),
            mass_matrix_global_inverse: Matrix::<f32>::zeros(3, 3),
            m_total: 0.,
            frc: vec![0., 0., 0.],
            tau: vec![0., 0., 0.],
        }
    }
    pub fn pre_compute_values(&mut self) {
        // total mass
        self.m_total = self.m.iter().fold(0., |total_mass, x| total_mass + x);

        // total mass has to be greater than zero
        assert!(self.m_total > 0.);

        // center of mass
        for i in 0..self.x.len() {
            self.r_com[0] += self.m[i] * self.x[i];
            self.r_com[1] += self.m[i] * self.y[i];
            self.r_com[2] += self.m[i] * self.z[i];
        }
        // divide by total mass to get the com
        self.r_com[0] /= self.m_total;
        self.r_com[1] /= self.m_total;
        self.r_com[2] /= self.m_total;

        // orientation is initially similar to the global axis
        self.orientation = matrix![1., 0., 0.;
                                   0., 1., 0.;
                                   0., 0., 1.];

        self.orientation_transpose = matrix![1., 0., 0.;
                                             0., 1., 0.;
                                             0., 0., 1.];
        // initialize position vector of points in body frame
        for i in 0..self.x.len() {
            self.x_body[i] = self.x[i] - self.r_com[0];
            self.y_body[i] = self.y[i] - self.r_com[1];
            self.z_body[i] = self.z[i] - self.r_com[2];
            self.x_body_global[i] = self.x[i] - self.r_com[0];
            self.y_body_global[i] = self.y[i] - self.r_com[1];
            self.z_body_global[i] = self.z[i] - self.r_com[2];
        }

        // initialize body moment of inertia
        let mut i_11 = 0.;
        let mut i_22 = 0.;
        let mut i_33 = 0.;
        let mut i_12 = 0.;
        let mut i_13 = 0.;
        let mut i_23 = 0.;

        for i in 0..self.x.len() {
            i_11 += self.m[i] * (self.y_body[i].powf(2.) + self.z_body[i].powf(2.));
            i_22 += self.m[i] * (self.x_body[i].powf(2.) + self.z_body[i].powf(2.));
            i_33 += self.m[i] * (self.x_body[i].powf(2.) + self.y_body[i].powf(2.));
            i_12 += -self.m[i] * self.x_body[i] * self.y_body[i];
            i_13 += -self.m[i] * self.x_body[i] * self.z_body[i];
            i_23 += -self.m[i] * self.y_body[i] * self.z_body[i];
        }

        {
            let mom_of_ine = &mut self.mass_matrix_body;
            mom_of_ine[[0, 0]] = i_11;
            mom_of_ine[[1, 1]] = i_22;
            mom_of_ine[[2, 2]] = i_33;
            mom_of_ine[[0, 1]] = i_12;
            mom_of_ine[[0, 2]] = i_13;
            mom_of_ine[[1, 2]] = i_23;

            // symmetrix elements
            mom_of_ine[[1, 0]] = mom_of_ine[[0, 1]];
            mom_of_ine[[2, 0]] = mom_of_ine[[0, 2]];
            mom_of_ine[[2, 1]] = mom_of_ine[[1, 2]];
            // TODO: take inverse of the matrix
        }
        self.mass_matrix_body_inverse = self
            .mass_matrix_body
            .clone()
            .inverse()
            .expect("This matrix should have an inverse!");
        // depending on the orientation calculate the moment
        // of inertia in global axis
        self.mass_matrix_global_inverse =
            &self.orientation * &self.mass_matrix_body_inverse * &self.orientation_transpose;

        // evaluate angular momentum from angular velocity and moment of inertia
        // Here I should use moment of inertial in global axis but at time t = 0
        // they both are equal
        self.ang_mom = (&self.mass_matrix_body
            * matrix![self.omega[0]; self.omega[1]; self.omega[2]])
            .into_vec();
    }
}

impl NNPS for RigidBody {
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
