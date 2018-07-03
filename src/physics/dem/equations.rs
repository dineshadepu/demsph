use super::DemDiscrete;
use cm::{InnerSpace, Vector3};
use contact_search::{get_neighbours_ll_2d, get_neighbours_ll_3d, LinkedListGrid};
use integrate::RK2;
use math::unit_vector_from_dx;
use std::collections::HashMap;

pub fn make_forces_zero(entity: &mut DemDiscrete) {
    for i in 0..entity.len {
        entity.fx[i] = 0.;
        entity.fy[i] = 0.;
        entity.fz[i] = 0.;
        entity.taux[i] = 0.;
        entity.tauy[i] = 0.;
        entity.tauz[i] = 0.;
    }
}

pub fn body_force_dem(entity: &mut DemDiscrete, gx: f32, gy: f32, gz: f32) {
    for i in 0..entity.len {
        entity.fx[i] += entity.m[i] * gx;
        entity.fy[i] += entity.m[i] * gy;
        entity.fz[i] += entity.m[i] * gz;
    }
}

/// Relative velocity between two particles
///
/// Find relative velocity of particle i with respect to particle
/// j at contact point
///
/// Relative velocity arises due to two components, one is
/// due to the linear velocity of the particle
/// and the other is due to angular velocity
///
/// $v_{ij} = v_i - v_j + (R_i \omega_i + R_j \omega_j) \times n_{ij}$
///
/// To find relative velocity due to angular component, we need
/// to take cross product of angular velocity with normal vector.
///
/// Given velocity of particle i and j with velocity and normal
/// passing from i to j, we find the relative velocity of particle i
/// with respect to particle j at contact point.

/// # Example
/// ```
/// # extern crate granules;
/// # extern crate cgmath;
/// # use cgmath::Vector3;
/// # use granules::physics::dem::equations::relative_velocity;
/// # use granules::math::{vec_compare};
/// let vi = Vector3::new(1., 0., 0.); // linear velocity of particle i
/// let vj = Vector3::new(-1., 0., 0.); // linear velocity of particle j
/// let ang_i = Vector3::new(0., 0., 1.); // angular velocity of particle i
/// let ang_j = Vector3::new(0., 0., 1.); // angular velocity of particle j
/// let rad_i = 1.; // radius of i
/// let rad_j = 1.; // radius of j
/// let n_ij = Vector3::new(1., 0., 0.); // normal vector from i to j
/// let rel_v = relative_velocity(vi, vj, ang_i, ang_j, n_ij, rad_i, rad_j);
/// let expected = Vector3::new(2., 2., 0.);
/// assert_eq!(vec_compare(&rel_v, &expected), true);
/// ```
pub fn relative_velocity(
    vi: Vector3<f32>,
    vj: Vector3<f32>,
    ang_v_i: Vector3<f32>,
    ang_v_j: Vector3<f32>,
    nij: Vector3<f32>,
    rad_i: f32,
    rad_j: f32,
) -> Vector3<f32> {
    vi - vj + (rad_i * ang_v_i + rad_j * ang_v_j).cross(nij)
}
/// Linear dashpot model introduced by Cundall and Strack.
pub fn linear_viscoelastic_model_other_dem(
    dst: &mut DemDiscrete,
    src: &mut DemDiscrete,
    kn: f32,
    mu: f32,
    dt: f32,
    grid: &LinkedListGrid,
    dim: usize,
) {
    // Select the neighbours function according to the dimention
    let get_nbrs = if dim == 3 {
        get_neighbours_ll_3d
    } else {
        get_neighbours_ll_2d
    };

    // Select the neighbours function according to the dimention
    for i in 0..dst.len {
        // position of particle i
        let pos_i = Vector3::new(dst.x[i], dst.y[i], dst.z[i]);
        // linear velocity of particle i
        let vel_i = Vector3::new(dst.u[i], dst.v[i], dst.w[i]);
        // angular velocity of particle i
        let ang_vel_i = Vector3::new(dst.omega_x[i], dst.omega_y[i], dst.omega_z[i]);

        let nbrs = get_nbrs([dst.x[i], dst.y[i], dst.z[i]], &grid, &src.id);

        for sub_view in nbrs {
            // neighbour indices j
            for &j in sub_view {
                // position of particle j in source
                let pos_j = Vector3::new(src.x[j], src.y[j], src.z[j]);
                // velocity of particle j
                let vel_j = Vector3::new(src.u[j], src.v[j], src.w[j]);
                // angular velocity of particle j
                let ang_vel_j = Vector3::new(src.omega_x[j], src.omega_y[j], src.omega_z[j]);

                // find the unit vector from i to j
                let dx = pos_j.x - pos_i.x;
                let dy = pos_j.y - pos_i.y;
                let dz = pos_j.z - pos_i.z;

                let distance = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).sqrt();
                let nij = unit_vector_from_dx(dx, dy, dz, distance);

                // relative velocity
                let v_ij = relative_velocity(
                    vel_i, vel_j, ang_vel_i, ang_vel_j, nij, dst.rad[i], src.rad[j],
                ); // this is vector

                // relative  normal velocity
                let v_n = v_ij.dot(nij) * nij; //this is vector
                // relative  tangential velocity
                let v_t = v_ij - v_n; //this is vector

                // radius sum
                let radsum = dst.rad[i] + src.rad[j];

                // overlap amount
                let delta_n = radsum - distance;

                // check if particles are in overlap
                if delta_n > 0. {

                }
                // if they are not overlapping, remove the particle j of src id
                // from history of particle i
                else {
                    // get the history of all particles being tracked by
                    // particle i
                    let hist = &mut dst.tang_overlap;

                    // the variable above (hist) is an array. Each element
                    // contains information of particles it is in contact

                    // The particle i's neighbours looks like

                    // hist[9] = {'0': {'2': Vector3, '31': Vector3, '7': Vector3},
                    //           '1': {'3': Vector3, '5': Vector3, '9': Vector3}}

                    // The meaning of above format is, particle i is is already
                    // contact with an entity with id '0' and '1'.

                    // diving little deep gives us the indices of those entites
                    // as, particle of index '9' has neighbours [2, 31, 7] of
                    // entity '0'. And also has neighbours [3, 5, 9] of entity
                    // '1'.

                }
            }
        }
    }
}

pub fn spring_force_self(
    entity: &mut DemDiscrete,
    kn: f32,
    mu: f32,
    dt: f32,
    grid: &LinkedListGrid,
    dim: usize,
) {
    let get_nbrs = if dim == 3 {
        get_neighbours_ll_3d
    } else {
        get_neighbours_ll_2d
    };

    for i in 0..entity.len {
        let nbrs = get_nbrs([entity.x[i], entity.y[i], entity.z[i]], &grid, &entity.id);

        for sub_view in nbrs {
            for &j in sub_view {
                if i != j {
                    let dx = entity.x[i] - entity.x[j];
                    let dy = entity.y[i] - entity.y[j];
                    let dz = entity.z[i] - entity.z[j];
                    let dist = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).powf(0.5);
                    let overlap = entity.rad[i] + entity.rad[j] - dist;

                    if overlap > 0. {
                        // force on i due to j
                        let mut fij_x = 0.;
                        let mut fij_y = 0.;
                        let mut fij_z = 0.;
                        let nx = dx / dist;
                        let ny = dy / dist;
                        let nz = dz / dist;

                        // relative velocity of contact point w.r.t particle j
                        let alpha_i = entity.rad[i] - overlap / 2.;
                        let alpha_j = entity.rad[j] - overlap / 2.;
                        let du = entity.u[i] - entity.u[j]
                            + alpha_i * (ny * entity.omega_z[i] - nz * entity.omega_y[i])
                            + alpha_j * (ny * entity.omega_z[j] - nz * entity.omega_y[j]);
                        let dv = entity.v[i] - entity.v[j]
                            + alpha_i * (nz * entity.omega_x[i] - nx * entity.omega_z[i])
                            + alpha_j * (nz * entity.omega_z[j] - nx * entity.omega_y[j]);
                        let dw = entity.w[i] - entity.w[j]
                            + alpha_i * (nx * entity.omega_y[i] - ny * entity.omega_x[i])
                            + alpha_j * (nx * entity.omega_y[j] - ny * entity.omega_x[j]);

                        let vndot = -(du * nx + dv * ny + dw * nz);

                        // compute normal force due to j on i
                        // use linear force model
                        // contact force in normal direction
                        let fdotn = kn * overlap;

                        // add normal contact force to total contact force
                        fij_x += fdotn * nx;
                        fij_y += fdotn * ny;
                        fij_z += fdotn * nz;

                        // ----------------------------
                        // Tangential force computation
                        // Check if the particles have friction.
                        // No friction means no tangential force computation
                        if mu != 0. {
                            // find the relative tangential velocity
                            let vt_x = du + nx * vndot;
                            let vt_y = dv + ny * vndot;
                            let vt_z = dw + nz * vndot;
                            // let the tangential force be
                            let mut ft = vec![0., 0., 0.];
                            // get the tangential spring connected from i to j.
                            // if no spring is available, that implies this is the
                            // first time interaction. So we will add a spring.

                            // check if the particle at index has the corresponding
                            // source id
                            match entity.tang_overlap[i].contains_key(&entity.id) {
                                true => {
                                    // now check if particle j is there in the list
                                    // of tangential overlaps
                                    match entity.tang_overlap[i][&entity.id].contains_key(&j) {
                                        true => {
                                            // this implies that the particle j is already in contact
                                            // with particle i
                                            // So we can simply increment the tangential overlap of particle
                                            // j spring connected to i

                                            // first project the tangential spring onto the current plane
                                            let mut delta_t = &mut entity.tang_overlap[i]
                                                .get_mut(&entity.id)
                                                .unwrap()
                                                .get_mut(&j)
                                                .unwrap();
                                            let delta_t_dot_n =
                                                delta_t[0] * nx + delta_t[1] * ny + delta_t[2] * nz;
                                            delta_t[0] -= delta_t_dot_n * nx;
                                            delta_t[1] -= delta_t_dot_n * ny;
                                            delta_t[2] -= delta_t_dot_n * nz;

                                            // find the tangential force using the spring projected
                                            // onto current tangential plane
                                            // FIXME: Dissipation is needed to add
                                            let mut ft_0 = vec![0., 0., 0.];
                                            // temporary dissipation coefficient
                                            let disp = 10.;
                                            ft_0[0] = -kn * delta_t[0] - disp * vt_x;
                                            ft_0[1] = -kn * delta_t[1] - disp * vt_y;
                                            ft_0[2] = -kn * delta_t[2] - disp * vt_z;
                                            let ft_0_magn = (ft_0[0].powf(2.)
                                                + ft_0[1].powf(2.)
                                                + ft_0[2].powf(2.))
                                                .sqrt();
                                            let norm_fn = fdotn.abs();

                                            // tangential direction will be
                                            let tx = ft_0[0] / ft_0_magn;
                                            let ty = ft_0[1] / ft_0_magn;
                                            let tz = ft_0[2] / ft_0_magn;

                                            let f_couloumb = mu * norm_fn;

                                            if ft_0_magn < f_couloumb {
                                                // Add the tangential force to totoal tangential
                                                // force of i due to j
                                                ft[0] += ft_0[0];
                                                ft[1] += ft_0[1];
                                                ft[2] += ft_0[2];

                                                // using the current velocity increment the the
                                                // spring for next iteration.
                                                // This is similar to integrating the position
                                                // dt will depend on the stage of integrator.
                                                // If we are in the first stage of RK2, then
                                                // the time step will be dt / 2. So an appropriate
                                                // time step is needed to be provided by the user.
                                                delta_t[0] += vt_x * dt;
                                                delta_t[1] += vt_y * dt;
                                                delta_t[2] += vt_z * dt;
                                            } else {
                                                // Add the tangential force to totoal tangential
                                                // force of i due to j
                                                ft[0] += f_couloumb * tx;
                                                ft[1] += f_couloumb * ty;
                                                ft[2] += f_couloumb * tz;

                                                // FIXME: Differentiate between slding and
                                                // dynamic friction

                                                // adjust the spring length to corresponding
                                                // tangential force in dynamic friction
                                                delta_t[0] = -f_couloumb * tx / kn;
                                                delta_t[1] = -f_couloumb * ty / kn;
                                                delta_t[2] = -f_couloumb * tz / kn;
                                            }

                                            // add
                                        }
                                        false => {
                                            // since tangential overlap particle i of entity id,
                                            // doesn't already have a contact with particle j of
                                            // another entity
                                            // We have to add such id to particle i's tangential
                                            // overlap attribute
                                            assert_eq!(
                                                entity.tang_overlap[i]
                                                    .get_mut(&entity.id)
                                                    .unwrap()
                                                    .insert(j, vec![0., 0., 0.]),
                                                None
                                            );
                                            // Note: No force is computed as this is first time overlap

                                            // increment the tangential overlap to next time step
                                            // using the current velocity

                                            let mut delta_t = &mut entity.tang_overlap[i]
                                                .get_mut(&entity.id)
                                                .unwrap()
                                                .get_mut(&j)
                                                .unwrap();
                                            delta_t[0] += vt_x * dt;
                                            delta_t[1] += vt_y * dt;
                                            delta_t[2] += vt_z * dt;
                                        }
                                    }
                                }
                                false => {
                                    // this case implies that for particle i, there was never
                                    // a connection to the entity

                                    // so first we need to add the outer hashmap
                                    assert_eq!(
                                        entity.tang_overlap[i].insert(entity.id, HashMap::new()),
                                        None
                                    );
                                    // now add the particle j to the list of
                                    // entity of particle i
                                    assert_eq!(
                                        entity.tang_overlap[i]
                                            .get_mut(&entity.id)
                                            .unwrap()
                                            .insert(j, vec![0., 0., 0.]),
                                        None
                                    );
                                    // Note: No force is computed as this is first time overlap

                                    // increment the tangential overlap to next time step
                                    // using the current velocity

                                    let mut delta_t = &mut entity.tang_overlap[i]
                                        .get_mut(&entity.id)
                                        .unwrap()
                                        .get_mut(&j)
                                        .unwrap();
                                    delta_t[0] += vt_x * dt;
                                    delta_t[1] += vt_y * dt;
                                    delta_t[2] += vt_z * dt;
                                }
                            }
                        }

                        entity.fx[i] += fij_x;
                        entity.fy[i] += fij_y;
                        entity.fz[i] += fij_z;
                    } else {
                        // Since this element is not in overlap,
                        // remove it from the tangential tracking of the
                        // particles, if it is been tracked

                        // check if the particle j is in the contact list
                        // of particle i
                        match entity.tang_overlap[i].contains_key(&entity.id) {
                            true => {
                                match entity.tang_overlap[i][&entity.id].contains_key(&j) {
                                    true => {
                                        // remove the particle j from tracking
                                        let _deleted = entity.tang_overlap[i]
                                            .get_mut(&entity.id)
                                            .unwrap()
                                            .remove(&j);
                                    }
                                    _ => {}
                                }
                            }
                            _ => {}
                        }
                    }
                }
            }
        }
    }
}

pub fn spring_force_other(
    destination: &mut DemDiscrete,
    source: &mut DemDiscrete,
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
                    destination.fx[i] += kn * overlap * nx;
                    destination.fy[i] += kn * overlap * ny;
                    destination.fz[i] += kn * overlap * nz;
                }
            }
        }
    }
}

impl RK2 for DemDiscrete {
    fn initialize(&mut self, dt: f32) {
        for i in 0..self.x.len() {
            self.x0[i] = self.x[i];
            self.y0[i] = self.y[i];
            self.z0[i] = self.z[i];
            self.u0[i] = self.u[i];
            self.v0[i] = self.v[i];
            self.w0[i] = self.w[i];
            self.omega_x0[i] = self.omega_x[i];
            self.omega_y0[i] = self.omega_y[i];
            self.omega_z0[i] = self.omega_z[i];
        }
    }
    fn stage1(&mut self, dt: f32) {
        let dtb2 = dt / 2.;
        for i in 0..self.x.len() {
            // propagate particles to next half time step
            self.x[i] = self.x0[i] + self.u[i] * dtb2;
            self.y[i] = self.y0[i] + self.v[i] * dtb2;
            self.z[i] = self.z0[i] + self.w[i] * dtb2;
            self.u[i] = self.u0[i] + self.fx[i] * self.m_inv[i] * dtb2;
            self.v[i] = self.v0[i] + self.fy[i] * self.m_inv[i] * dtb2;
            self.w[i] = self.w0[i] + self.fz[i] * self.m_inv[i] * dtb2;
            self.omega_x[i] = self.omega_x0[i] + self.taux[i] * self.i_inv[i] * dtb2;
            self.omega_y[i] = self.omega_y0[i] + self.tauy[i] * self.i_inv[i] * dtb2;
            self.omega_z[i] = self.omega_z0[i] + self.tauz[i] * self.i_inv[i] * dtb2;
        }
    }
    fn stage2(&mut self, dt: f32) {
        for i in 0..self.x.len() {
            // propagate particles to next time step
            self.x[i] = self.x0[i] + self.u[i] * dt;
            self.y[i] = self.y0[i] + self.v[i] * dt;
            self.z[i] = self.z0[i] + self.w[i] * dt;
            self.u[i] = self.u0[i] + self.fx[i] * self.m_inv[i] * dt;
            self.v[i] = self.v0[i] + self.fy[i] * self.m_inv[i] * dt;
            self.w[i] = self.w0[i] + self.fz[i] * self.m_inv[i] * dt;
            self.omega_x[i] = self.omega_x0[i] + self.taux[i] * self.i_inv[i] * dt;
            self.omega_y[i] = self.omega_y0[i] + self.tauy[i] * self.i_inv[i] * dt;
            self.omega_z[i] = self.omega_z0[i] + self.tauz[i] * self.i_inv[i] * dt;
        }
    }
}

#[test]
fn test_unit_vector() {}
