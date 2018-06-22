use super::DemDiscrete;
use contact_search::{LinkedListGrid, get_neighbours_ll_2d, get_neighbours_ll_3d};
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

pub fn integrate(entity: &mut DemDiscrete, dt: f32) {
    for i in 0..entity.len {
        entity.u[i] += entity.fx[i] * entity.m_inv[i] * dt;
        entity.v[i] += entity.fy[i] * entity.m_inv[i] * dt;
        entity.w[i] += entity.fz[i] * entity.m_inv[i] * dt;
        entity.x[i] += entity.u[i] * dt;
        entity.y[i] += entity.v[i] * dt;
        entity.z[i] += entity.w[i] * dt;
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
                                            let ft_0_magn = (ft_0[0].powf(2.) + ft_0[1].powf(2.)
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
