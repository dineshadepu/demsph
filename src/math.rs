use cm::Vector3;


pub fn vec_compare(va: &Vector3<f32>, vb: &Vector3<f32>) -> bool {
    if {va.x != vb.x} || {va.y != vb.y} || {va.z != vb.z}{
        return false
    }
    true
}

/// Given dx, dy, dz and distance between two vectors
/// compute the unit vector passing through them.
///
/// This is used in discrete element method to find the normal passing from two
/// particles
///
/// # Examples

/// ```
/// # extern crate granules;
/// # extern crate cgmath;
/// # use granules::math::{unit_vector_from_dx, vec_compare};
/// # use cgmath::Vector3;
/// let vec1: Vec<f32> = vec![1., 4., 6.];
/// let vec2: Vec<f32> = vec![2., 5., 9.];
/// let dx = vec1[0] - vec2[0];
/// let dy = vec1[1] - vec2[1];
/// let dz = vec1[2] - vec2[2];
/// let magn = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).sqrt();
/// // Unit vector from vec2 to vec1
/// let unit_vec = unit_vector_from_dx(dx, dy, dz, magn);
/// let expected = Vector3::new(-0.30151135, -0.30151135, -0.904534);
/// assert_eq!(vec_compare(&unit_vec, &expected), true);
/// ```

pub fn unit_vector_from_dx(dx: f32, dy: f32, dz: f32, magn: f32) -> Vector3<f32> {
    Vector3::new(dx / magn, dy / magn, dz / magn)
}

/// Compute unit vector from vec1 to vec2.
///
/// This is used in discrete element method to find the normal passing from two
/// particles
///
/// # Examples

/// ```
/// # extern crate granules;
/// # extern crate cgmath;
/// # use granules::math::{unit_vector_from_dx, vec_compare,
/// #                      unit_vector_from_vector};
/// # use cgmath::Vector3;
/// let vec1: Vec<f32> = vec![1., 4., 6.];
/// let vec2: Vec<f32> = vec![2., 5., 9.];
/// // Unit vector from vec2 to vec1
/// let unit_vec = unit_vector_from_vector(vec2, vec1);
/// let expected = Vector3::new(-0.30151135, -0.30151135, -0.904534);
/// assert_eq!(vec_compare(&unit_vec, &expected), true);
/// ```
pub fn unit_vector_from_vector(from: Vec<f32>, to: Vec<f32>) -> Vector3<f32> {
    let dx = to[0] - from[0];
    let dy = to[1] - from[1];
    let dz = to[2] - from[2];
    let magn = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).sqrt();

    unit_vector_from_dx(dx, dy, dz, magn)
}

pub fn unit_vector_from_point(from: Vector3<f32>, to: Vector3<f32>) -> Vector3<f32> {
    let dx = to.x - from.x;
    let dy = to.y - from.y;
    let dz = to.z - from.z;
    let distance = (dx.powf(2.) + dy.powf(2.) + dz.powf(2.)).sqrt();

    unit_vector_from_dx(dx, dy, dz, distance)
}
