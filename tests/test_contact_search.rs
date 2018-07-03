extern crate granules;

use granules::contact_search::{get_neighbours_ll_2d, LinkedListGrid};
use granules::physics::dem::DemDiscrete;

fn extract_indices_from_reference(nbrs_reference: Vec<&Vec<usize>>) -> Vec<usize> {
    let mut nbrs = vec![];
    for subview in nbrs_reference {
        for &i in subview {
            nbrs.push(i);
        }
    }
    nbrs
}
#[test]
fn test_grid_cells_indices() {
    //
    //
    //
    //
    //           |-------------------------------------------|
    //           |          |         |          |           |
    //           |          |         |          |           |
    //           |          |         |          |           |
    //           |          |         |          |           |
    //           |__________|_________|__________|___________|
    //           |          |         |          |           |
    //           |          | rad 1   |          |           |
    //           |          |_________|          |           |
    //           |          |         |          |           |
    //           |__________|_________|__________|___________|
    //           |          |         |          |           |
    //           |          |         |          |           |
    //           | 2        |         | 2        |           |
    //           |. id 2    |         |. id 2    |           |
    //           |__________|_________|__________|___________|
    //           |          |         |          |           |
    //           |          |         |          |           |
    //           |          |         |          |           |
    //           | 1        |         | 1        |           |
    //           |. id 1    |         |.  id 2   |           |
    //           |-------------------------------------------|
    // check with dem entities
    let mut ent1 = DemDiscrete::new(2, 1, String::from("ent1"));
    ent1.h = vec![1., 1.];
    ent1.x = vec![0., 2.];
    ent1.y = vec![0., 0.];
    let mut ent2 = DemDiscrete::new(2, 2, String::from("ent2"));
    ent2.h = vec![1., 1.];
    ent2.x = vec![0., 2.];
    ent2.y = vec![2., 2.];

    // neighbours for a grid size of h=1., same as radius of particle
    let scale = 1.;
    let grid = LinkedListGrid::new(&mut vec![&mut ent1, &mut ent2], scale);

    let nbrs_reference = get_neighbours_ll_2d([ent1.x[0], ent1.y[0], 0.], &grid, &2);
    let nbrs = extract_indices_from_reference(nbrs_reference);
    assert_eq!(0, nbrs.len());
    let nbrs_reference = get_neighbours_ll_2d([ent1.x[0], ent1.y[0], 0.], &grid, &1);
    let nbrs = extract_indices_from_reference(nbrs_reference);
    assert_eq!(1, nbrs.len());

    // neighbours for a grid size of 2 * h=2., double the radius
    let scale = 2.;
    let grid = LinkedListGrid::new(&mut vec![&mut ent1, &mut ent2], scale);

    let nbrs_reference = get_neighbours_ll_2d([ent1.x[0], ent1.y[0], 0.], &grid, &2);
    let nbrs = extract_indices_from_reference(nbrs_reference);
    assert_eq!(2, nbrs.len());
    let nbrs_reference = get_neighbours_ll_2d([ent1.x[0], ent1.y[0], 0.], &grid, &1);
    let nbrs = extract_indices_from_reference(nbrs_reference);
    assert_eq!(2, nbrs.len());
}

#[test]
fn test_grid_cells_indices_a_corner_case() {
    // check a corner case:
    // every particle is in single cell. Make sure that indices are not repeating
    // check with dem entities
    let mut ent1 = DemDiscrete::new(3, 1, String::from("ent1"));
    ent1.h = vec![1., 1., 1.];
    ent1.x = vec![0., 0.1, 0.];
    ent1.y = vec![0., 0., 0.1];
    let mut ent2 = DemDiscrete::new(3, 2, String::from("ent2"));
    ent2.h = vec![1., 1., 1.];
    ent2.x = vec![-0.1, 0., 0.1];
    ent2.y = vec![0., -0.1, 0.1];

    // neighbours for a grid size of h=1., same as radius of particle
    let scale = 1.;
    let grid = LinkedListGrid::new(&mut vec![&mut ent1, &mut ent2], scale);

    let nbrs_reference = get_neighbours_ll_2d([ent1.x[0], ent1.y[0], 0.], &grid, &2);
    let nbrs = extract_indices_from_reference(nbrs_reference);
    assert_eq!(3, nbrs.len());
    let nbrs_reference = get_neighbours_ll_2d([ent1.x[0], ent1.y[0], 0.], &grid, &1);
    let nbrs = extract_indices_from_reference(nbrs_reference);
    assert_eq!(3, nbrs.len());

    let expected_nbrs = vec![0, 1, 2];
    for i in expected_nbrs {
        assert_eq!(true, nbrs.contains(&i));
    }
}
