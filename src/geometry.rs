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

pub fn zeros_like(x: &Vec<f32>) -> Vec<f32> {
    let y = vec![0.; x.len()];
    return y;
}

pub fn ones_like(x: &Vec<f32>) -> Vec<f32> {
    let y = vec![1.; x.len()];
    return y;
}
