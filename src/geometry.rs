pub fn linspace(start: f32, stop: f32, num: isize) -> Vec<f32> {
    if num < 0 {
        panic!("Number of samples, {}, must be non-negative.", num);
    }
    let range = stop - start;
    let divide = num as f32;
    let dx = range / divide;
    let mut tmp = 1;
    let mut y = vec![start];
    let mut next = start + dx;
    while tmp <= num {
        y.push(next);
        next = next + dx;
        tmp = tmp + 1;
    }
    return y;
}

pub fn arange(start: f32, stop: f32, spacing: f32) -> Vec<f32> {
    if spacing < 0. {
        panic!("Spacing {} can't be negative", spacing);
    }
    let mut y = Vec::new();
    let mut tmp = start;
    while tmp <= stop {
        y.push(tmp);
        tmp = tmp + spacing;
    }
    return y;
}

pub fn grid_2d(length: f32, height: f32, spacing: f32) -> (Vec<f32>, Vec<f32>) {
    let x = arange(0., length, spacing);
    let y = arange(0., height, spacing);

    // join them to get a 2d grid
    let mut x_grid = Vec::new();
    let mut y_grid = Vec::new();
    for i in &x {
        for j in &y {
            x_grid.push(*i);
            y_grid.push(*j);
        }
    }
    (x_grid, y_grid)
}

pub fn dam_break_2d_geometry(
    f_l: f32,
    f_h: f32,
    f_s: f32,
    t_l: f32,
    t_h: f32,
    t_s: f32,
    layers: usize,
) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
    if layers < 1 {
        panic!("Layers cant be less than 1");
    }
    let (mut xf, mut yf) = grid_2d(f_l, f_h, f_s);
    xf.iter_mut()
        .for_each(|x| *x = *x + (2. * layers as f32) * f_s);
    yf.iter_mut()
        .for_each(|y| *y = *y + (2. * layers as f32) * f_s);

    let (xtmp, ytmp) = grid_2d(t_l, t_h, t_s);

    // remove layers of the tank
    let x = arange(0., t_l, t_s);
    // find maximum of a vector
    let mut max = 0.;
    for i in &x {
        if max < *i {
            max = *i;
        }
    }

    // save all the particles in  the bottom of tank
    let mut xt = Vec::new();
    let mut yt = Vec::new();
    for (a, b) in xtmp.iter().zip(ytmp.iter()) {
        if *b < (layers - 1) as f32 * t_s + t_s / 2. {
            xt.push(*a);
            yt.push(*b);
        } else {
            if *a < (layers - 1) as f32 * t_s + t_s / 2.
                || *a > max - (layers - 1) as f32 * t_s - t_s / 2.
            {
                xt.push(*a);
                yt.push(*b);
            }
        }
    }
    (xf, yf, xt, yt)
}

pub fn get_3d_block(
    block_length: f32,
    block_height: f32,
    block_width: f32,
    block_spacing: f32,
) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let mut w = 0.;

    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut z = Vec::new();
    while w < block_width {
        let (mut x_tmp, mut y_tmp) = grid_2d(block_length, block_height, block_spacing);
        let mut z_tmp = vec![w; x_tmp.len()];

        // append the local slice at width w to global x, y, z
        x.append(&mut x_tmp);
        y.append(&mut y_tmp);
        z.append(&mut z_tmp);

        // increase the width to new layer
        w += block_spacing;
    }
    (x, y, z)
}

pub fn get_3d_tank(
    tank_length: f32,
    tank_height: f32,
    tank_width: f32,
    tank_spacing: f32,
    tank_layers: usize,
    layers_outside: bool,
) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    // change the type of tank layer
    let tank_layers = tank_layers as f32;

    let (x, y, z) = if layers_outside {
        // increase the length width and height of tank block
        let tank_height_extended = tank_height + tank_layers * tank_spacing;
        let tank_length_extended = tank_length + tank_layers * tank_spacing;
        let tank_width_extended = tank_width + tank_layers * 2. * tank_spacing;

        let (xt, yt, zt) = get_3d_block(
            tank_length_extended,
            tank_height_extended,
            tank_width_extended,
            tank_spacing,
        );

        // now remove the particles of the whole tank
        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut z = Vec::new();
        let x_left_limit = 0. + tank_spacing / 10. + (tank_layers - 1.) * tank_spacing;
        let y_left_limit = 0. + tank_spacing / 10. + (tank_layers - 1.) * tank_spacing;
        let z_left_limit = 0. + tank_spacing / 10. + (tank_layers - 1.) * tank_spacing;
        for (&x_tmp, &y_tmp, &z_tmp) in izip!(&xt, &yt, &zt) {
            if x_tmp < x_left_limit || x_tmp > tank_length {
                x.push(x_tmp);
                y.push(y_tmp);
                z.push(z_tmp);
            } else {
                if y_tmp < y_left_limit {
                    x.push(x_tmp);
                    y.push(y_tmp);
                    z.push(z_tmp);
                } else {
                    if z_tmp < z_left_limit || z_tmp > tank_width {
                        x.push(x_tmp);
                        y.push(y_tmp);
                        z.push(z_tmp);
                    }
                }
            }
        }
        (x, y, z)
    } else {
        let (xt, yt, zt) = get_3d_block(tank_length, tank_height, tank_width, tank_spacing);
        // now remove the particles of the whole tank
        let mut x: Vec<f32> = Vec::new();
        let mut y: Vec<f32> = Vec::new();
        let mut z: Vec<f32> = Vec::new();
        let x_left_limit = 0. + tank_spacing / 10. + (tank_layers - 1.) * tank_spacing;
        let x_right_limit = tank_length - tank_spacing / 10. - (tank_layers - 1.) * tank_spacing;
        let y_left_limit = 0. + tank_spacing / 10. + (tank_layers - 1.) * tank_spacing;
        let _y_right_limit = tank_height - tank_spacing / 10. - (tank_layers - 1.) * tank_spacing;
        let z_left_limit = 0. + tank_spacing / 10. + (tank_layers - 1.) * tank_spacing;
        let z_right_limit = tank_width - tank_spacing / 10. - (tank_layers - 1.) * tank_spacing;
        for (&x_tmp, &y_tmp, &z_tmp) in izip!(&xt, &yt, &zt) {
            if x_tmp < x_left_limit || x_tmp > x_right_limit {
                x.push(x_tmp);
                y.push(y_tmp);
                z.push(z_tmp);
            } else {
                if y_tmp < y_left_limit {
                    x.push(x_tmp);
                    y.push(y_tmp);
                    z.push(z_tmp);
                } else {
                    if z_tmp < z_left_limit || z_tmp > z_right_limit {
                        x.push(x_tmp);
                        y.push(y_tmp);
                        z.push(z_tmp);
                    }
                }
            }
        }
        (x, y, z)
    };
    (x, y, z)
}
pub fn dam_break_3d_geometry(
    f_l: f32,
    f_h: f32,
    f_w: f32,
    f_s: f32,
    t_l: f32,
    t_h: f32,
    t_w: f32,
    t_s: f32,
    layers: usize,
    layers_outside: bool,
) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
    // --------------------------------------------
    // --------------------------------------------
    // first create a 3d block of fluid
    let (xf, yf, zf) = get_3d_block(f_l, f_h, f_w, f_s);

    // now create tank block
    let (xt, yt, zt) = get_3d_tank(t_l, t_h, t_w, t_s, layers, layers_outside);

    // Now remove the layers
    (xf, yf, zf, xt, yt, zt)
}

pub fn zeros_like(x: &Vec<f32>) -> Vec<f32> {
    let y = vec![0.; x.len()];
    return y;
}

pub fn ones_like(x: &Vec<f32>) -> Vec<f32> {
    let y = vec![1.; x.len()];
    return y;
}
