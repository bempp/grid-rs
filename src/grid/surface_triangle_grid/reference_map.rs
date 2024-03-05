//! Implementation of the reference Map

use num::Float;
use rlst_common::types::Scalar;
use rlst_dense::{
    array::SliceArray,
    rlst_array_from_slice2,
    traits::{DefaultIterator, Shape},
};

use crate::traits::{GridType, ReferenceMapType};

use super::grid::TriangleSurfaceGrid;

pub struct TriangleReferenceMap<'a, T: Float + Scalar> {
    reference_points: SliceArray<'a, T, 2>,
    grid: &'a TriangleSurfaceGrid<T>,
}

impl<'a, T: Float + Scalar> TriangleReferenceMap<'a, T> {
    pub fn new(
        reference_points: &'a [T],
        grid: &'a TriangleSurfaceGrid<T>,
    ) -> Self {
        assert_eq!(reference_points.len() % 2, 0);
        let npoints = reference_points.len() / 2;
        Self {
            reference_points: rlst_array_from_slice2!(T, reference_points, [2, npoints]),
            grid,
        }
    }
}

impl<'a, T: Float + Scalar> ReferenceMapType for TriangleReferenceMap<'a, T> {
    type Grid = TriangleSurfaceGrid<T>;

    fn domain_dimension(&self) -> usize {
        2
    }

    fn physical_dimension(&self) -> usize {
        3
    }

    fn number_of_reference_points(&self) -> usize {
        self.reference_points.shape()[1]
    }

    fn reference_to_physical(&self, cell_index: usize, point_index: usize, value: &mut [<Self::Grid as GridType>::T]) {
        assert_eq!(value.len(), 3);
        let v0 = self.grid.vertices[self.grid.cells[cell_index][0]];
        let jacobian = &self.grid.jacobians[cell_index];
        for (index, val_out) in value.iter_mut().enumerate() {
            *val_out = v0[index]
                + jacobian[[index, 0]] * self.reference_points[[0, point_index]]
                + jacobian[[index, 1]] * self.reference_points[[1, point_index]];
        }
    }

    fn jacobian(&self, cell_index: usize, _point_index: usize, value: &mut [<Self::Grid as GridType>::T]) {
        for (elem_in, val_out) in self.grid.jacobians[cell_index]
            .iter()
            .zip(value.iter_mut())
        {
            *val_out = elem_in
        }
    }

    fn normal(&self, cell_index: usize, _point_index: usize, value: &mut [<Self::Grid as GridType>::T]) {
        for (elem_in, val_out) in self.grid.normals[cell_index]
            .iter()
            .zip(value.iter_mut())
        {
            *val_out = elem_in
        }
    }
}
