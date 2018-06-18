use super::DemDiscrete;
use contact_search::{LinkedListGrid, get_neighbours_ll_2d, get_neighbours_ll_3d};

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

pub fn integrate(entity : &mut DemDiscrete, dt: f32) {
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
                        let nx = dx / dist;
                        let ny = dy / dist;
                        let nz = dz / dist;
                        entity.fx[i] += kn * overlap * nx;
                        entity.fy[i] += kn * overlap * ny;
                        entity.fz[i] += kn * overlap * nz;
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

        let nbrs = get_nbrs([destination.x[i], destination.y[i], destination.z[i]], &grid, &source.id);

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
