//! Functionality common to multiple grid implementations

use crate::grid::traits::Geometry;
use num::Float;
use rlst_common::types::Scalar;
use rlst_dense::{
    array::Array,
    traits::{Shape, UnsafeRandomAccessByRef, UnsafeRandomAccessByValue},
};

/// Compute a physical point
pub fn compute_point<T: Scalar, Table: UnsafeRandomAccessByRef<4, Item = T> + Shape<4>>(
    geometry: &impl Geometry<T = T>,
    table: Table,
    cell_index: usize,
    point_index: usize,
    point: &mut [T],
) {
    assert_eq!(geometry.dim(), point.len());

    let cell = geometry.index_map()[cell_index];

    for component in point.iter_mut() {
        *component = T::from(0.0).unwrap();
    }
    for (i, v) in geometry.cell_points(cell).unwrap().iter().enumerate() {
        let t = unsafe { *table.get_unchecked([0, point_index, i, 0]) };
        for (j, component) in point.iter_mut().enumerate() {
            *component += *geometry.coordinate(*v, j).unwrap() * t;
        }
    }
}

/// Compute a Jacobian
pub fn compute_jacobian<T: Scalar, Table: UnsafeRandomAccessByRef<4, Item = T> + Shape<4>>(
    geometry: &impl Geometry<T = T>,
    table: Table,
    tdim: usize,
    cell_index: usize,
    point_index: usize,
    jacobian: &mut [T],
) {
    let gdim = geometry.dim();
    assert_eq!(jacobian.len(), gdim * tdim);

    let cell = geometry.index_map()[cell_index];

    for component in jacobian.iter_mut() {
        *component = T::from(0.0).unwrap();
    }
    for (i, v) in geometry.cell_points(cell).unwrap().iter().enumerate() {
        for gd in 0..gdim {
            for td in 0..tdim {
                jacobian[td * gdim + gd] += *geometry.coordinate(*v, gd).unwrap()
                    * unsafe { *table.get_unchecked([1 + td, point_index, i, 0]) };
            }
        }
    }
}

/// Compute a normal from a Jacobian of a cell with topological dimension 2 and geometric dimension 3
pub fn compute_normal_from_jacobian23<T: Scalar>(jacobian: &[T], normal: &mut [T]) {
    assert_eq!(jacobian.len(), 6);
    assert_eq!(normal.len(), 3);

    for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
        normal[i] = jacobian[j] * jacobian[3 + k] - jacobian[k] * jacobian[3 + j];
    }
    let size = Scalar::sqrt(normal.iter().map(|&x| Scalar::powi(x, 2)).sum::<T>());
    for i in normal.iter_mut() {
        *i /= size;
    }
}

/// Compute the diameter of a triangle
pub fn compute_diameter_triangle<
    T: Scalar<Real = T> + Float,
    ArrayImpl: UnsafeRandomAccessByValue<1, Item = T> + Shape<1>,
>(
    v0: Array<T, ArrayImpl, 1>,
    v1: Array<T, ArrayImpl, 1>,
    v2: Array<T, ArrayImpl, 1>,
) -> T {
    let a = (v0.view() - v1.view()).norm_2();
    let b = (v0 - v2.view()).norm_2();
    let c = (v1 - v2).norm_2();
    Scalar::sqrt((b + c - a) * (a + c - b) * (a + b - c) / (a + b + c))
}

/// Compute the diameter of a quadrilateral
pub fn compute_diameter_quadrilateral<
    T: Scalar<Real = T> + Float,
    ArrayImpl: UnsafeRandomAccessByValue<1, Item = T> + Shape<1>,
>(
    v0: Array<T, ArrayImpl, 1>,
    v1: Array<T, ArrayImpl, 1>,
    v2: Array<T, ArrayImpl, 1>,
    v3: Array<T, ArrayImpl, 1>,
) -> T {
    T::max((v0 - v3).norm_2(), (v1 - v2).norm_2())
}
