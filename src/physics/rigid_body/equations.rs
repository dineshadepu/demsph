use super::RigidBody;
use contact_search::{get_neighbours_ll_2d, get_neighbours_ll_3d, LinkedListGrid};
use integrate::RK2;

pub fn make_forces_zero_rigid_body(entities: &mut Vec<&mut RigidBody>) {
    for entity in entities {
        entity.frc[0] = 0.;
        entity.frc[1] = 0.;
        entity.frc[2] = 0.;
        entity.tau[0] = 0.;
        entity.tau[1] = 0.;
        entity.tau[2] = 0.;

        for i in 0..entity.x.len() {
            entity.fx[i] = 0.;
            entity.fy[i] = 0.;
            entity.fz[i] = 0.;
        }
    }
}

pub fn body_force_rigid_body(entities: &mut Vec<&mut RigidBody>, gx: f32, gy: f32, gz: f32) {
    for entity in entities {
        // since body is concentrated at centre of mass.
        // we can directly add the force to the centre of mass
        // position, as there would be no torque
        entity.frc[0] = entity.m_total * gx;
        entity.frc[1] = entity.m_total * gy;
        entity.frc[2] = entity.m_total * gz;
    }
}

pub fn spring_force_rigid_other(
    destination: &mut RigidBody,
    source: &mut RigidBody,
    kn: f32,
    grid: &LinkedListGrid,
    dim: usize,
) {
    let get_nbrs = if dim == 3 {
        get_neighbours_ll_3d
    } else {
        get_neighbours_ll_2d
    };

    for i in 0..destination.len {
        let nbrs = get_nbrs(
            [destination.x[i], destination.y[i], destination.z[i]],
            &grid,
            &source.id,
        );

        for sub_view in nbrs {
            for &j in sub_view {
                let dx = destination.x[i] - source.x[j];
                let dy = destination.y[i] - source.y[j];
                let dz = destination.z[i] - source.z[j];
                let dist = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).powf(0.5);
                let overlap = destination.rad[i] + source.rad[j] - dist;
                if overlap > 0. {
                    let nx = dx / dist;
                    let ny = dy / dist;
                    let nz = dz / dist;
                    // damping force in normal direction
                    let du = destination.u[i] - source.u[j];
                    let dv = destination.v[i] - source.v[j];
                    let dw = destination.w[i] - source.w[j];
                    let v_n = -(du * nx + dv * ny + dw * nz);

                    destination.fx[i] += (kn * overlap + 1000. * v_n) * nx;
                    destination.fy[i] += (kn * overlap + 1000. * v_n) * ny;
                    destination.fz[i] += (kn * overlap + 1000. * v_n) * nz;
                }
            }
        }
    }
}

pub fn aggregate_forces_moments(entities: &mut Vec<&mut RigidBody>) {
    for entity in entities {
        for i in 0..entity.fx.len() {
            // add the force on each particle to global force
            entity.frc[0] += entity.fx[i];
            entity.frc[1] += entity.fy[i];
            entity.frc[2] += entity.fz[i];

            // find the torque due to force on each particle
            entity.tau[0] +=
                entity.y_body_global[i] * entity.fz[i] - entity.z_body_global[i] * entity.fy[i];
            entity.tau[1] +=
                entity.z_body_global[i] * entity.fx[i] - entity.x_body_global[i] * entity.fz[i];
            entity.tau[2] +=
                entity.x_body_global[i] * entity.fy[i] - entity.y_body_global[i] * entity.fx[i];
        }
    }
}

impl RK2 for RigidBody {
    fn initialize(&mut self, dt: f32) {}
    fn stage1(&mut self, dt: f32) {
        let dtb2 = dt / 2.;
        self.r_com[0] = self.r_com[0] + dtb2 * self.v_com[0];
        self.r_com[1] = self.r_com[1] + dtb2 * self.v_com[1];
        self.r_com[2] = self.r_com[2] + dtb2 * self.v_com[2];

        self.v_com[0] = self.v_com[0] + dtb2 * self.frc[0] / self.m_total;
        self.v_com[1] = self.v_com[1] + dtb2 * self.frc[1] / self.m_total;
        self.v_com[2] = self.v_com[2] + dtb2 * self.frc[2] / self.m_total;

        // rotate the axis to next time
        self.omega_matrix = matrix![0., -self.omega[2], self.omega[1];
                                    self.omega[2], 0., -self.omega[0];
                                    -self.omega[1], self.omega[0], 0.];
        self.orientation = &self.orientation + &self.omega_matrix * &self.orientation * dtb2;
        // update the angular momentum, it is a vector Vec<f32>
        self.ang_mom[0] = self.ang_mom[0] + self.tau[0] * dtb2;
        self.ang_mom[1] = self.ang_mom[1] + self.tau[1] * dtb2;
        self.ang_mom[2] = self.ang_mom[2] + self.tau[2] * dtb2;

        // Since orientation is updated, to compute the MoI global at
        // next time step we need to do
        self.mass_matrix_global_inverse =
            &self.orientation * &self.mass_matrix_body_inverse * &self.orientation_transpose;

        // We updated angular momentum to next time step, and We know MoI at
        // that time, Lets update the angular velocity to next step using
        // these two
        self.omega = (&self.mass_matrix_body_inverse
            * matrix![self.ang_mom[0]; self.ang_mom[1]; self.ang_mom[2]])
            .into_vec();

        // now update quantites of each particle
        for i in 0..self.x.len() {
            // rotated local position vector from center of mass
            // to particle in global axis would be
            let pos_global_i =
                (&self.orientation * matrix![self.x_body[i]; self.y_body[i]; self.z_body[i]])
                    .into_vec();
            self.x_body_global[i] = pos_global_i[0];
            self.y_body_global[i] = pos_global_i[1];
            self.z_body_global[i] = pos_global_i[2];

            // now update the position in global axis
            self.x[i] = self.r_com[0] + pos_global_i[0];
            self.y[i] = self.r_com[1] + pos_global_i[1];
            self.z[i] = self.r_com[2] + pos_global_i[2];

            // now update the position in global axis
            // angular velocity of the particle is
            self.u[i] =
                self.v_com[0] + (self.omega[1] * pos_global_i[2] - self.omega[2] * pos_global_i[1]);
            self.v[i] =
                self.v_com[1] + (self.omega[2] * pos_global_i[1] - self.omega[0] * pos_global_i[2]);
            self.w[i] =
                self.v_com[2] + (self.omega[0] * pos_global_i[1] - self.omega[1] * pos_global_i[0]);
        }
    }
    fn stage2(&mut self, dt: f32) {}
}
