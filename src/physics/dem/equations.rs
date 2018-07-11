use super::DemDiscrete;
use cm::{InnerSpace, Vector3 as V3, dot, Zero};
use contact_search::{LinkedListGrid, get_neighbours_ll_2d, get_neighbours_ll_3d};
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
    vi: V3<f32>,
    vj: V3<f32>,
    ang_v_i: V3<f32>,
    ang_v_j: V3<f32>,
    nij: V3<f32>,
    rad_i: f32,
    rad_j: f32,
) -> V3<f32> {
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
        let pos_i = V3::new(dst.x[i], dst.y[i], dst.z[i]);
        // linear velocity of particle i
        let vel_i = V3::new(dst.u[i], dst.v[i], dst.w[i]);
        // angular velocity of particle i
        let ang_vel_i = V3::new(dst.omega_x[i], dst.omega_y[i], dst.omega_z[i]);

        let nbrs = get_nbrs([dst.x[i], dst.y[i], dst.z[i]], &grid, &src.id);

        for sub_view in nbrs {
            // neighbour indices j
            for &j in sub_view {
                // position of particle j in source
                let pos_j = V3::new(src.x[j], src.y[j], src.z[j]);
                // velocity of particle j
                let vel_j = V3::new(src.u[j], src.v[j], src.w[j]);
                // angular velocity of particle j
                let ang_vel_j = V3::new(src.omega_x[j], src.omega_y[j], src.omega_z[j]);

                // find the unit vector from i to j
                let dx = pos_j.x - pos_i.x;
                let dy = pos_j.y - pos_i.y;
                let dz = pos_j.z - pos_i.z;

                let distance = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).sqrt();
                let nij = unit_vector_from_dx(dx, dy, dz, distance);

                // relative velocity
                let v_ij = relative_velocity(
                    vel_i,
                    vel_j,
                    ang_vel_i,
                    ang_vel_j,
                    nij,
                    dst.rad[i],
                    src.rad[j],
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
                    let hist = &mut dst.tang_overlap[i];

                    // -------------------------------------------------
                    // Tangential overlap variable explanation
                    // -------------------------------------------------

                    // the variable above (hist) contains particles already in
                    // overlap.

                    // The particle i's neighbours looks like

                    // hist = {'0': {'2': Vector3, '31': Vector3, '7': Vector3},
                    //           '1': {'3': Vector3, '5': Vector3, '9': Vector3}}

                    // The meaning of above format is, particle i is is already
                    // contact with an entities with id's '0' and '1'.

                    // diving little deep gives us the indices of those entites
                    // as, particle of index '9' has neighbours [2, 31, 7] of
                    // entity '0'. And also has neighbours [3, 5, 9] of entity
                    // '1'.

                    // So the type of hist would be
                    // Vec<HashMap<usize, HashMap<usize, Vector3>>>
                    // -------------------------------------------------
                    // Tangential overlap variable explanation
                    // -------------------------------------------------

                    // Find the index j in hist of particle i with id of
                    // src.id

                    // If j is already been tracked then remove it
                    // If it is not been tracked then leave hist alone

                    // Note: To do this operatio I am using match
                    // match is provided by rust and it's awesome

                    // this leaf is to check if the particle i has history
                    // with the src
                    match hist.contains_key(&src.id) {
                        // If it has neighbours with src, then
                        // go ahead and check if it is tracking particle j
                        true => {
                            match hist[&src.id].contains_key(&j) {
                                // If it has particle j
                                // remove it
                                true => {
                                    hist.get_mut(&src.id).unwrap().remove(&j);
                                }

                                // if it doesn't have particle index j, then
                                // leave it alone
                                false => {}
                            };
                        }

                        // if it doesn't have src id, then leave it alone
                        false => {}
                    };
                }
            }
        }
    }
}

pub fn spring_force_self(
    e: &mut DemDiscrete,
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

    for i in 0..e.len {
        let nbrs = get_nbrs([e.x[i], e.y[i], e.z[i]], &grid, &e.id);

        for sub_view in nbrs {
            for &j in sub_view {
                if i != j {
                    let dx = e.x[i] - e.x[j];
                    let dy = e.y[i] - e.y[j];
                    let dz = e.z[i] - e.z[j];
                    let dist = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).powf(0.5);
                    let overlap = e.rad[i] + e.rad[j] - dist;

                    if overlap > 0. {
                        // force on i due to j
                        let mut fij = V3::new(0., 0., 0.);
                        // normal vector from j to i -> n_ij
                        let nij = V3::new(dx / dist, dy / dist, dz / dist);

                        // angular velocity of particle i and j
                        let ang_i = V3::new(e.omega_x[i], e.omega_y[i], e.omega_z[i]);
                        let ang_j = V3::new(e.omega_x[j], e.omega_y[j], e.omega_z[j]);

                        // -----------------------------
                        // relative velocity of contact point of i w.r.t particle j
                        // Reduced radius from overlap is
                        let alpha_i = e.rad[i] - overlap / 2.;
                        let alpha_j = e.rad[j] - overlap / 2.;

                        // cross product of angular velocity with normal
                        let ang_cross_n = nij.cross(alpha_i * ang_i + alpha_j * ang_j);

                        // relative velocity vij
                        let vij = V3::new(e.u[i], e.v[i], e.w[i])
                            - V3::new(e.u[j], e.v[j], e.w[j])
                            + ang_cross_n;

                        let vndot = - dot(vij, nij);

                        // compute normal force due to j on i
                        // use linear force model
                        // contact force in normal direction
                        let fdotn = kn * overlap + 0.001 * vndot;

                        // add normal contact force to total contact force
                        fij += fdotn * nij;

                        // ----------------------------
                        // Tangential force computation
                        // Check if the particles have friction.
                        // No friction means no tangential force computation
                        if mu != 0. {
                            // find the relative tangential velocity
                            let vt = vij + nij * vndot;
                            // let the tangential force be
                            let mut ft:V3<f32> = V3::zero();

                            // get the tangential spring connected from i to j.
                            // if no spring is available, that implies this is the
                            // first time interaction. So we will add a spring.

                            // check if the particle at index has the corresponding
                            // source (entities) id
                            match e.tang_overlap[i].contains_key(&e.id) {
                                true => {
                                    // now check if particle j is there in the list
                                    // of tangential overlaps
                                    match e.tang_overlap[i][&e.id].contains_key(&j) {
                                        true => {
                                            // this implies that the particle j is already in contact
                                            // with particle i
                                            // So we can simply increment the tangential overlap of particle
                                            // j spring connected to i

                                            // first project the tangential spring onto the current plane
                                            let mut delta_t = &mut e.tang_overlap[i]
                                                .get_mut(&e.id)
                                                .unwrap()
                                                .get_mut(&j)
                                                .unwrap();
                                            let delta_t_dot_n = dot(**delta_t, nij);
                                            **delta_t -= delta_t_dot_n * nij;

                                            // find the tangential force using the spring projected
                                            // onto current tangential plane
                                            // FIXME: Dissipation is needed to add
                                            let mut ft_0 = V3::zero();
                                            // temporary dissipation coefficient
                                            let disp = 10.;
                                            ft_0 = -kn * **delta_t - disp * vt;
                                            let ft_0_magn = ft_0.magnitude();
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
                                                **delta_t += vt * dt;
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
                                            // since tangential overlap particle i of e id,
                                            // doesn't already have a contact with particle j of
                                            // another e
                                            // We have to add such id to particle i's tangential
                                            // overlap attribute
                                            assert_eq!(
                                                e.tang_overlap[i]
                                                    .get_mut(&e.id)
                                                    .unwrap()
                                                    .insert(j, V3::new(0., 0., 0.)),
                                                None
                                            );
                                            // Note: No force is computed as this is first time overlap

                                            // increment the tangential overlap to next time step
                                            // using the current velocity

                                            let mut delta_t = &mut e.tang_overlap[i]
                                                .get_mut(&e.id)
                                                .unwrap()
                                                .get_mut(&j)
                                                .unwrap();
                                            **delta_t += vt * dt;
                                        }
                                    }
                                }
                                false => {
                                    // this case implies that for particle i, there was never
                                    // a connection to the e

                                    // so first we need to add the outer hashmap
                                    assert_eq!(
                                        e.tang_overlap[i].insert(e.id, HashMap::new()),
                                        None
                                    );
                                    // now add the particle j to the list of
                                    // e of particle i
                                    assert_eq!(
                                        e.tang_overlap[i]
                                            .get_mut(&e.id)
                                            .unwrap()
                                            .insert(j, V3::new(0., 0., 0.)),
                                        None
                                    );
                                    // Note: No force is computed as this is first time overlap

                                    // increment the tangential overlap to next time step
                                    // using the current velocity

                                    let mut delta_t = &mut e.tang_overlap[i]
                                        .get_mut(&e.id)
                                        .unwrap()
                                        .get_mut(&j)
                                        .unwrap();
                                    **delta_t += vt * dt;
                                }
                            }
                        }

                        e.fx[i] += fij[0];
                        e.fy[i] += fij[1];
                        e.fz[i] += fij[2];
                    } else {
                        // Since this element is not in overlap,
                        // remove it from the tangential tracking of the
                        // particles, if it is been tracked

                        // check if the particle j is in the contact list
                        // of particle i
                        match e.tang_overlap[i].contains_key(&e.id) {
                            true => {
                                match e.tang_overlap[i][&e.id].contains_key(&j) {
                                    true => {
                                        // remove the particle j from tracking
                                        let _deleted =
                                            e.tang_overlap[i].get_mut(&e.id).unwrap().remove(&j);
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
