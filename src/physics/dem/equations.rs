use super::{GetMutDEMDest, GetMutDEMSrc};
use contact_search::{get_neighbours_ll_2d, get_neighbours_ll_3d, LinkedListGrid};

pub fn make_forces_zero<T: GetMutDEMDest>(entity: &mut T) {
    let ent1 = entity.get_mut_parts();

    for i in 0..*ent1.len {
        ent1.fx[i] = 0.;
        ent1.fy[i] = 0.;
        ent1.fz[i] = 0.;
        ent1.taux[i] = 0.;
        ent1.tauy[i] = 0.;
        ent1.tauz[i] = 0.;
    }
}

pub fn body_force_dem<T: GetMutDEMDest>(entity: &mut T, gx: f32, gy: f32, gz: f32) {
    let ent1 = entity.get_mut_parts();

    for i in 0..*ent1.len {
        ent1.fx[i] += ent1.m[i] * gx;
        ent1.fy[i] += ent1.m[i] * gy;
        ent1.fz[i] += ent1.m[i] * gz;
    }
}

pub fn integrate<T: GetMutDEMDest>(src: &mut T, dt: f32) {
    let entity = src.get_mut_parts();
    for i in 0..*entity.len {
        entity.u[i] += entity.fx[i] * entity.m_inv[i] * dt;
        entity.v[i] += entity.fy[i] * entity.m_inv[i] * dt;
        entity.w[i] += entity.fz[i] * entity.m_inv[i] * dt;
        entity.x[i] += entity.u[i] * dt;
        entity.y[i] += entity.v[i] * dt;
        entity.z[i] += entity.w[i] * dt;
    }
}

pub fn spring_force_self<T: GetMutDEMDest>(entity: &mut T, kn: f32, grid: &LinkedListGrid, dim: usize) {
    let dst = entity.get_mut_parts();

    let get_nbrs = if dim == 2 {
        get_neighbours_ll_2d
    }
    else {
        get_neighbours_ll_3d
    };

    for i in 0..*dst.len {
        // let nbrs = if dim == 2 {
        //     get_neighbours_ll_3d([dst.x[i], dst.y[i], dst.z[i]], &grid, &dst.id)
        // }
        // else {
        //     get_neighbours_ll_3d([dst.x[i], dst.y[i], dst.z[i]], &grid, &dst.id)
        // };
        let nbrs = get_nbrs([dst.x[i], dst.y[i], dst.z[i]], &grid, &dst.id);


        for j in nbrs {
            if i != j {
                let dx = dst.x[i] - dst.x[j];
                let dy = dst.y[i] - dst.y[j];
                let dz = dst.z[i] - dst.z[j];
                let dist = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).powf(0.5);
                let overlap = dst.rad[i] + dst.rad[j] - dist;

                if overlap > 0. {
                    let nx = dx / dist;
                    let ny = dy / dist;
                    let nz = dz / dist;
                    dst.fx[i] += kn * overlap * nx;
                    dst.fy[i] += kn * overlap * ny;
                    dst.fz[i] += kn * overlap * nz;
                }
            }
        }
    }
}

pub fn spring_force_other<T: GetMutDEMDest, U: GetMutDEMSrc>(
    destination: &mut T,
    source: &mut U,
    kn: f32,
    grid: &LinkedListGrid,
    dim: usize,
) {
    let dst = destination.get_mut_parts();
    let src = source.get_mut_parts();
    for i in 0..*dst.len {
        let nbrs = get_neighbours_ll_3d([dst.x[i], dst.y[i], dst.z[i]], &grid, &src.id);
        for j in nbrs {
            let dx = dst.x[i] - src.x[j];
            let dy = dst.y[i] - src.y[j];
            let dz = dst.z[i] - src.z[j];
            let dist = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).powf(0.5);
            let overlap = dst.rad[i] + src.rad[j] - dist;
            if overlap > 0. {
                let nx = dx / dist;
                let ny = dy / dist;
                let nz = dz / dist;
                dst.fx[i] += kn * overlap * nx;
                dst.fy[i] += kn * overlap * ny;
                dst.fz[i] += kn * overlap * nz;
            }
        }
    }
}
