use super::{DemDiscrete, DemDiscreteDstTrait, DemDiscreteSrcTrait};
use cm::{dot, InnerSpace, Vector3 as V3, Zero};
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
pub fn linear_viscoelastic_model_dem_other<T, U>(
    dst: &mut T,
    src: &mut U,
    kn: f32,
    mu: f32,
    dt: f32,
    stage: usize,
    grid: &LinkedListGrid,
    dim: usize,
) where
    T: DemDiscreteDstTrait,
    U: DemDiscreteSrcTrait,
{
    // Select the neighbours function according to the dimention
    let get_nbrs = if dim == 3 {
        get_neighbours_ll_3d
    } else {
        get_neighbours_ll_2d
    };

    let dest = dst.get_parts_mut();
    let srce = src.get_parts_mut();

    // Select the neighbours function according to the dimention
    for i in 0..*dest.len {
        // position of particle i
        let pos_i = V3::new(dest.x[i], dest.y[i], dest.z[i]);
        // linear velocity of particle i
        let vel_i = V3::new(dest.u[i], dest.v[i], dest.w[i]);
        // angular velocity of particle i
        let ang_vel_i = V3::new(dest.omega_x[i], dest.omega_y[i], dest.omega_z[i]);

        let nbrs = get_nbrs([dest.x[i], dest.y[i], dest.z[i]], &grid, &srce.id);

        for sub_view in nbrs {
            // neighbour indices j
            for &j in sub_view {
                // position of particle j in source
                let pos_j = V3::new(srce.x[j], srce.y[j], srce.z[j]);
                // velocity of particle j
                let vel_j = V3::new(srce.u[j], srce.v[j], srce.w[j]);
                // angular velocity of particle j
                let ang_vel_j = V3::new(srce.omega_x[j], srce.omega_y[j], srce.omega_z[j]);

                // find the unit vector from j to i
                let dx = pos_i.x - pos_j.x;
                let dy = pos_i.y - pos_j.y;
                let dz = pos_i.z - pos_j.z;

                let distance = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).sqrt();
                // radius sum
                let radsum = dest.rad[i] + srce.rad[j];

                // overlap amount
                let delta_n = radsum - distance;

                // check if particles are in overlap
                if delta_n > 0. {
                    // Define the force variables, total, tangential, torsion
                    let mut f = V3::zero();
                    let mut f_t: V3<f32> = V3::zero();
                    let mut f_r: V3<f32> = V3::zero();

                    // normal vector
                    let nij = unit_vector_from_dx(dx, dy, dz, distance);

                    // relative velocity
                    let v_ij = relative_velocity(
                        vel_i,
                        vel_j,
                        ang_vel_i,
                        ang_vel_j,
                        nij,
                        dest.rad[i],
                        srce.rad[j],
                    ); // this is vector

                    // relative  normal velocity
                    let v_n = v_ij.dot(nij) * nij; //this is vector
                    // relative  tangential velocity
                    let v_t = v_ij - v_n; //this is vector

                    // ----------------------------------------------------
                    // Normal force with damping
                    // FIX ME: Need to use real coefficients
                    let f_n = kn * delta_n * nij - v_n * 0.001;

                    // Add normal force to total force with damping in normal direction
                    f += f_n;

                    // ----------------------------------------------------
                    // ----------------Tangential force -------------------
                    // Check for tangential contacts only if there is friction
                    if mu != 0. {
                        // If tangential forces are present
                        // get the history of all particles being tracked by
                        // particle i
                        let hist = &mut dest.tang_history[i];

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
                        // srce.id

                        // If j is already been tracked then remove it
                        // If it is not been tracked then leave hist alone

                        // Note: To do this operatio I am using match
                        // match is provided by rust and it's awesome

                        // this leaf is to check if the particle i has history
                        // with the srce entity
                        match hist.contains_key(&srce.id) {
                            // If it has neighbours with srce, then
                            // go ahead and check if it is tracking particle j
                            true => {
                                match hist[&srce.id].contains_key(&j) {
                                    // If it already has particle j
                                    // Then find the force due to its current
                                    // tangential spring elongation
                                    true => {
                                        // Get the spring
                                        let mut tang_overlap =
                                            hist.get_mut(&srce.id).unwrap().get_mut(&j).unwrap();

                                        // Now project the spring onto current tangential plane
                                        // http://www.piko.ovgu.de/piko_media/aktuelles/Siegen/LUDING2012PIKO_Contacts.pdf
                                        let tang_overlap_rotated =
                                            *tang_overlap - nij * (dot(*tang_overlap, nij));

                                        // Find tangential test force from the rotated spring
                                        let f_t0 = -1e4 * tang_overlap_rotated - 0.001 * v_t;
                                        let f_t0_magn = f_t0.magnitude();

                                        // find the tangential unit vector
                                        let t_ij = if f_t0_magn > 0. {
                                            f_t0 / f_t0.magnitude()
                                        } else {
                                            V3::zero()
                                        };

                                        // Check for sliding
                                        let fn_norm = f_n.magnitude();
                                        if f_t0_magn <= mu * fn_norm {
                                            // Set the tangential force to test tangential
                                            // force
                                            f_t = f_t0;

                                            // Increment the tangential spring
                                            // for next time step
                                            // Note: dt changes for different stages
                                            if stage == 1 {
                                                *tang_overlap = tang_overlap_rotated + v_t * dt;
                                            } else if stage == 2 {
                                                // use the tangential overlap at time t i.e., tang_overlap0
                                                // project it onto current orientaton, i.e., t + dt / 2
                                                let hist0 = &mut dest.tang_history0[i];
                                                let tang_overlap0 = hist0
                                                    .get_mut(&srce.id)
                                                    .unwrap()
                                                    .get_mut(&j)
                                                    .unwrap();
                                                let tang_overlap_rotated0 = *tang_overlap0
                                                    - nij * (dot(*tang_overlap, nij));
                                                *tang_overlap = tang_overlap_rotated0 + v_t * dt;
                                                *tang_overlap0 = *tang_overlap;
                                            }
                                        } else {
                                            // So the particles slide.
                                            // Set the tangential force to the maximum force allowed
                                            // by Couloumb force

                                            // FIX ME: Change friction coefficient to sliding
                                            f_t = mu * fn_norm * t_ij;

                                            // Restrict the spring length such that the resultant
                                            // tangential force equals Couloumb force
                                            // Note: dt changes for different stages
                                            if stage == 1 {
                                                *tang_overlap = (f_t + 0.001 * v_t) / 1e4;
                                            } else if stage == 2 {
                                                // use the tangential overlap at time t i.e., tang_overlap0
                                                // project it onto current orientaton, i.e., t + dt / 2
                                                let hist0 = &mut dest.tang_history0[i];
                                                let tang_overlap0 = hist0
                                                    .get_mut(&srce.id)
                                                    .unwrap()
                                                    .get_mut(&j)
                                                    .unwrap();
                                                *tang_overlap = (f_t + 0.001 * v_t) / 1e4;
                                                *tang_overlap0 = *tang_overlap;
                                            }
                                        }
                                    }

                                    // if it doesn't have particle index j, then
                                    // add the particle to the history
                                    false => {
                                        // create tangential_overlap0 and tangential_overlap
                                        let tang_overlap0 = V3::zero();

                                        // Since this is first time contact, we will not
                                        // have any tangential force
                                        // ----------Skip force calculation-------------

                                        // Increment spring to next time step

                                        if stage == 1 {
                                            let tang_overlap = v_t * dt;
                                            hist.get_mut(&srce.id).unwrap().insert(j, tang_overlap);
                                            let hist0 = &mut dest.tang_history0[i];
                                            hist0
                                                .get_mut(&srce.id)
                                                .unwrap()
                                                .insert(j, tang_overlap0);
                                        } else if stage == 2 {
                                            // use the tangential overlap at time t i.e., tang_overlap0
                                            // since this is first time contact there won't
                                            // be any force calculation
                                            let tang_overlap = v_t * dt;
                                            hist.get_mut(&srce.id).unwrap().insert(j, tang_overlap);
                                            let hist0 = &mut dest.tang_history0[i];
                                            // Note the difference between stage 1 and stage 2
                                            hist0
                                                .get_mut(&srce.id)
                                                .unwrap()
                                                .insert(j, tang_overlap);
                                        }
                                    }
                                };
                            }

                            // if it doesn't have srce id, add the srce id
                            false => {
                                // add srce id as source to tangential history
                                hist.insert(*srce.id, HashMap::new());

                                // And also in the temporary history
                                let hist0 = &mut dest.tang_history0[i];
                                hist0.insert(*srce.id, HashMap::new());

                                // ------------------------------------
                                // now add the particle j in srce_d hashmap

                                // create tangential_overlap0 and tangential_overlap
                                let tang_overlap0 = V3::zero();

                                // Since this is first time contact, we will not
                                // have any tangential force
                                // ----------Skip force calculation-------------

                                // Increment spring to next time step

                                if stage == 1 {
                                    let tang_overlap = v_t * dt;
                                    hist.get_mut(&srce.id).unwrap().insert(j, tang_overlap);
                                    hist0.get_mut(&srce.id).unwrap().insert(j, tang_overlap0);
                                } else if stage == 2 {
                                    // use the tangential overlap at time t i.e., tang_overlap0
                                    // since this is first time contact there won't
                                    // be any force calculation
                                    let tang_overlap = v_t * dt;
                                    hist.get_mut(&srce.id).unwrap().insert(j, tang_overlap);
                                    // Note the difference between stage 1 and stage 2
                                    hist0.get_mut(&srce.id).unwrap().insert(j, tang_overlap);
                                }
                            }
                        };
                    }
                    dest.fx[i] += f[0];
                    dest.fy[i] += f[1];
                    dest.fz[i] += f[2];
                }
                // if they are not overlapping, remove the particle j of srce id
                // from history of particle i
                else {
                    // Check for tangential contacts only if there is friction
                    if mu != 0. {
                        // get the history of all particles being tracked by
                        // particle i
                        let hist = &mut dest.tang_history[i];
                        let hist0 = &mut dest.tang_history0[i];

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
                        // srce.id

                        // If j is already been tracked then remove it
                        // If it is not been tracked then leave hist alone

                        // Note: To do this operatio I am using match
                        // match is provided by rust and it's awesome

                        // this leaf is to check if the particle i has history
                        // with the srce
                        match hist.contains_key(&srce.id) {
                            // If it has neighbours with srce, then
                            // go ahead and check if it is tracking particle j
                            true => {
                                match hist[&srce.id].contains_key(&j) {
                                    // If it has particle j
                                    // remove it
                                    true => {
                                        hist.get_mut(&srce.id).unwrap().remove(&j);
                                        hist0.get_mut(&srce.id).unwrap().remove(&j);
                                    }

                                    // if it doesn't have particle index j, then
                                    // leave it alone
                                    false => {}
                                };
                            }

                            // if it doesn't have srce id, then leave it alone
                            false => {}
                        };
                    }
                }
            }
        }
    }
}

/// Linear dashpot model introduced by Cundall and Strack.
pub fn linear_viscoelastic_model_dem_self<T>(
    dst: &mut T,
    kn: f32,
    mu: f32,
    dt: f32,
    stage: usize,
    grid: &LinkedListGrid,
    dim: usize,
) where
    T: DemDiscreteDstTrait,
{
    // Select the neighbours function according to the dimention
    let get_nbrs = if dim == 3 {
        get_neighbours_ll_3d
    } else {
        get_neighbours_ll_2d
    };

    let dest = dst.get_parts_mut();

    // Select the neighbours function according to the dimention
    for i in 0..*dest.len {
        // position of particle i
        let pos_i = V3::new(dest.x[i], dest.y[i], dest.z[i]);
        // linear velocity of particle i
        let vel_i = V3::new(dest.u[i], dest.v[i], dest.w[i]);
        // angular velocity of particle i
        let ang_vel_i = V3::new(dest.omega_x[i], dest.omega_y[i], dest.omega_z[i]);

        let nbrs = get_nbrs([dest.x[i], dest.y[i], dest.z[i]], &grid, &dest.id);

        for sub_view in nbrs {
            // neighbour indices j
            for &j in sub_view {
                if i != j {
                    // position of particle j in source
                    let pos_j = V3::new(dest.x[j], dest.y[j], dest.z[j]);
                    // velocity of particle j
                    let vel_j = V3::new(dest.u[j], dest.v[j], dest.w[j]);
                    // angular velocity of particle j
                    let ang_vel_j = V3::new(dest.omega_x[j], dest.omega_y[j], dest.omega_z[j]);

                    // find the unit vector from j to i
                    let dx = pos_i.x - pos_j.x;
                    let dy = pos_i.y - pos_j.y;
                    let dz = pos_i.z - pos_j.z;

                    let distance = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).sqrt();
                    // radius sum
                    let radsum = dest.rad[i] + dest.rad[j];

                    // overlap amount
                    let delta_n = radsum - distance;

                    // check if particles are in overlap
                    if delta_n > 0. {
                        // Define the force variables, total, tangential, torsion
                        let mut f = V3::zero();
                        let mut f_t: V3<f32> = V3::zero();
                        let mut f_r: V3<f32> = V3::zero();

                        // normal vector
                        let nij = unit_vector_from_dx(dx, dy, dz, distance);

                        // relative velocity
                        let v_ij = relative_velocity(
                            vel_i,
                            vel_j,
                            ang_vel_i,
                            ang_vel_j,
                            nij,
                            dest.rad[i],
                            dest.rad[j],
                        ); // this is vector

                        // relative  normal velocity
                        let v_n = v_ij.dot(nij) * nij; //this is vector
                        // relative  tangential velocity
                        let v_t = v_ij - v_n; //this is vector

                        // ----------------------------------------------------
                        // Normal force with damping
                        // FIX ME: Need to use real coefficients
                        let f_n = kn * delta_n * nij - v_n * 0.001;

                        // Add normal force to total force with damping in normal direction
                        f += f_n;

                        // ----------------------------------------------------
                        // ----------------Tangential force -------------------
                        // Check for tangential contacts only if there is friction
                        if mu != 0. {
                            // If tangential forces are present
                            // get the history of all particles being tracked by
                            // particle i
                            let hist = &mut dest.tang_history[i];

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
                            // srce.id

                            // If j is already been tracked then remove it
                            // If it is not been tracked then leave hist alone

                            // Note: To do this operatio I am using match
                            // match is provided by rust and it's awesome

                            // this leaf is to check if the particle i has history
                            // with the srce entity
                            match hist.contains_key(&dest.id) {
                                // If it has neighbours with srce, then
                                // go ahead and check if it is tracking particle j
                                true => {
                                    match hist[&dest.id].contains_key(&j) {
                                        // If it already has particle j
                                        // Then find the force due to its current
                                        // tangential spring elongation
                                        true => {
                                            // Get the spring
                                            let mut tang_overlap = hist.get_mut(&dest.id)
                                                .unwrap()
                                                .get_mut(&j)
                                                .unwrap();

                                            // Now project the spring onto current tangential plane
                                            // http://www.piko.ovgu.de/piko_media/aktuelles/Siegen/LUDING2012PIKO_Contacts.pdf
                                            let tang_overlap_rotated =
                                                *tang_overlap - nij * (dot(*tang_overlap, nij));

                                            // Find tangential test force from the rotated spring
                                            let f_t0 = -1e4 * tang_overlap_rotated - 0.001 * v_t;
                                            let f_t0_magn = f_t0.magnitude();

                                            // find the tangential unit vector
                                            let t_ij = if f_t0_magn > 0. {
                                                f_t0 / f_t0.magnitude()
                                            } else {
                                                V3::zero()
                                            };

                                            // Check for sliding
                                            let fn_norm = f_n.magnitude();
                                            if f_t0_magn <= mu * fn_norm {
                                                // Set the tangential force to test tangential
                                                // force
                                                f_t = f_t0;

                                                // Increment the tangential spring
                                                // for next time step
                                                // Note: dt changes for different stages
                                                if stage == 1 {
                                                    *tang_overlap = tang_overlap_rotated + v_t * dt;
                                                } else if stage == 2 {
                                                    // use the tangential overlap at time t i.e., tang_overlap0
                                                    // project it onto current orientaton, i.e., t + dt / 2
                                                    let hist0 = &mut dest.tang_history0[i];
                                                    let tang_overlap0 = hist0
                                                        .get_mut(&dest.id)
                                                        .unwrap()
                                                        .get_mut(&j)
                                                        .unwrap();
                                                    let tang_overlap_rotated0 = *tang_overlap0
                                                        - nij * (dot(*tang_overlap, nij));
                                                    *tang_overlap =
                                                        tang_overlap_rotated0 + v_t * dt;
                                                    *tang_overlap0 = *tang_overlap;
                                                }
                                            } else {
                                                // So the particles slide.
                                                // Set the tangential force to the maximum force allowed
                                                // by Couloumb force

                                                // FIX ME: Change friction coefficient to sliding
                                                f_t = mu * fn_norm * t_ij;

                                                // Restrict the spring length such that the resultant
                                                // tangential force equals Couloumb force
                                                // Note: dt changes for different stages
                                                if stage == 1 {
                                                    *tang_overlap = (f_t + 0.001 * v_t) / 1e4;
                                                } else if stage == 2 {
                                                    // use the tangential overlap at time t i.e., tang_overlap0
                                                    // project it onto current orientaton, i.e., t + dt / 2
                                                    let hist0 = &mut dest.tang_history0[i];
                                                    let tang_overlap0 = hist0
                                                        .get_mut(&dest.id)
                                                        .unwrap()
                                                        .get_mut(&j)
                                                        .unwrap();
                                                    *tang_overlap = (f_t + 0.001 * v_t) / 1e4;
                                                    *tang_overlap0 = *tang_overlap;
                                                }
                                            }
                                        }

                                        // if it doesn't have particle index j, then
                                        // add the particle to the history
                                        false => {
                                            // create tangential_overlap0 and tangential_overlap
                                            let tang_overlap0 = V3::zero();

                                            // Since this is first time contact, we will not
                                            // have any tangential force
                                            // ----------Skip force calculation-------------

                                            // Increment spring to next time step

                                            if stage == 1 {
                                                let tang_overlap = v_t * dt;
                                                hist.get_mut(&dest.id)
                                                    .unwrap()
                                                    .insert(j, tang_overlap);
                                                let hist0 = &mut dest.tang_history0[i];
                                                hist0
                                                    .get_mut(&dest.id)
                                                    .unwrap()
                                                    .insert(j, tang_overlap0);
                                            } else if stage == 2 {
                                                // use the tangential overlap at time t i.e., tang_overlap0
                                                // since this is first time contact there won't
                                                // be any force calculation
                                                let tang_overlap = v_t * dt;
                                                hist.get_mut(&dest.id)
                                                    .unwrap()
                                                    .insert(j, tang_overlap);
                                                let hist0 = &mut dest.tang_history0[i];
                                                // Note the difference between stage 1 and stage 2
                                                hist0
                                                    .get_mut(&dest.id)
                                                    .unwrap()
                                                    .insert(j, tang_overlap);
                                            }
                                        }
                                    };
                                }

                                // if it doesn't have srce id, add the srce id
                                false => {
                                    // add srce id as source to tangential history
                                    hist.insert(*dest.id, HashMap::new());

                                    // And also in the temporary history
                                    let hist0 = &mut dest.tang_history0[i];
                                    hist0.insert(*dest.id, HashMap::new());

                                    // ------------------------------------
                                    // now add the particle j in srce_d hashmap

                                    // create tangential_overlap0 and tangential_overlap
                                    let tang_overlap0 = V3::zero();

                                    // Since this is first time contact, we will not
                                    // have any tangential force
                                    // ----------Skip force calculation-------------

                                    // Increment spring to next time step

                                    if stage == 1 {
                                        let tang_overlap = v_t * dt;
                                        hist.get_mut(&dest.id).unwrap().insert(j, tang_overlap);
                                        hist0.get_mut(&dest.id).unwrap().insert(j, tang_overlap0);
                                    } else if stage == 2 {
                                        // use the tangential overlap at time t i.e., tang_overlap0
                                        // since this is first time contact there won't
                                        // be any force calculation
                                        let tang_overlap = v_t * dt;
                                        hist.get_mut(&dest.id).unwrap().insert(j, tang_overlap);
                                        // Note the difference between stage 1 and stage 2
                                        hist0.get_mut(&dest.id).unwrap().insert(j, tang_overlap);
                                    }
                                }
                            };
                            // add tangential force to total force of particle
                            f += f_t;
                        }
                        // add total force to global force
                        dest.fx[i] += f[0];
                        dest.fy[i] += f[1];
                        dest.fz[i] += f[2];
                    }
                    // if they are not overlapping, remove the particle j of srce id
                    // from history of particle i
                    else {
                        // Check for tangential contacts only if there is friction
                        if mu != 0. {
                            // get the history of all particles being tracked by
                            // particle i
                            let hist = &mut dest.tang_history[i];
                            let hist0 = &mut dest.tang_history0[i];

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
                            // srce.id

                            // If j is already been tracked then remove it
                            // If it is not been tracked then leave hist alone

                            // Note: To do this operatio I am using match
                            // match is provided by rust and it's awesome

                            // this leaf is to check if the particle i has history
                            // with the srce
                            match hist.contains_key(&dest.id) {
                                // If it has neighbours with srce, then
                                // go ahead and check if it is tracking particle j
                                true => {
                                    match hist[&dest.id].contains_key(&j) {
                                        // If it has particle j
                                        // remove it
                                        true => {
                                            hist.get_mut(&dest.id).unwrap().remove(&j);
                                            hist0.get_mut(&dest.id).unwrap().remove(&j);
                                        }

                                        // if it doesn't have particle index j, then
                                        // leave it alone
                                        false => {}
                                    };
                                }

                                // if it doesn't have srce id, then leave it alone
                                false => {}
                            };
                        }
                    }
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
