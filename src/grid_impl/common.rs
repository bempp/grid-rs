//! Functionality common to multiple grid implementations
use rlst_dense::traits::{UnsafeRandomAccessByRef, Shape};
use bempp_traits::element::FiniteElement;
use rlst_common::types::Scalar;
use crate::grid_impl::traits::Geometry;

pub fn compute_point<T: Scalar, Table: UnsafeRandomAccessByRef<4, Item = T> + Shape<4>>(
    geometry: &impl Geometry<T=T>, table: Table, cell_index: usize, point_index: usize, point: &mut [T]
) {
    let gdim = geometry.dim();
    let cell = geometry.index_map()[cell_index];
    let element_npts = geometry.cell_element(cell).unwrap().dim();
    assert_eq!(point.len(), gdim);

    for component in point.iter_mut() {
        *component = T::from(0.0).unwrap();
    }
    for (i, v) in geometry.cell_points(cell).unwrap()
        .iter()
        .enumerate()
    {
        let t = unsafe { *table.get_unchecked([0, point_index, i, 0]) };
        for (j, component) in point.iter_mut().enumerate() {
            *component += *geometry.coordinate(*v, j).unwrap() * t;
        }
    }

}
