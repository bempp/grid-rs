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
    cell_index: usize,
    grid: &'a TriangleSurfaceGrid<T>,
}

impl<'a, T: Float + Scalar> TriangleReferenceMap<'a, T> {
    pub fn new(
        reference_points: &'a [T],
        cell_index: usize,
        grid: &'a TriangleSurfaceGrid<T>,
    ) -> Self {
        assert_eq!(reference_points.len() % 2, 0);
        let npoints = reference_points.len() / 2;
        Self {
            reference_points: rlst_array_from_slice2!(T, reference_points, [2, npoints]),
            cell_index,
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

    fn reference_to_physical(&self, point_index: usize, value: &mut [<Self::Grid as GridType>::T]) {
        assert_eq!(value.len(), 3);
        let v0 = self.grid.vertices[self.grid.cells[self.cell_index][0]];
        let jacobian = &self.grid.jacobians[self.cell_index];
        for (index, val_out) in value.iter_mut().enumerate() {
            *val_out = v0[index]
                + jacobian[[index, 0]] * self.reference_points[[0, point_index]]
                + jacobian[[index, 1]] * self.reference_points[[1, point_index]];
        }
    }

    fn jacobian(&self, _point_index: usize, value: &mut [<Self::Grid as GridType>::T]) {
        for (elem_in, val_out) in self.grid.jacobians[self.cell_index]
            .iter()
            .zip(value.iter_mut())
        {
            *val_out = elem_in
        }
    }

    fn normal(&self, _point_index: usize, value: &mut [<Self::Grid as GridType>::T]) {
        for (elem_in, val_out) in self.grid.normals[self.cell_index]
            .iter()
            .zip(value.iter_mut())
        {
            *val_out = elem_in
        }
    }
}

pub struct TriangleReferenceMapIterator<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>>
{
    iter: Iter,
    reference_points: &'a [Grid::T],
    grid: &'a Grid,
}

impl<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>>
    TriangleReferenceMapIterator<'a, Grid, Iter>
{
    pub fn new(iter: Iter, reference_points: &'a [Grid::T], grid: &'a Grid) -> Self {
        Self {
            iter,
            reference_points,
            grid,
        }
    }
}

impl<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> Iterator
    for TriangleReferenceMapIterator<'a, Grid, Iter>
{
    type Item = Grid::ReferenceMap<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(cell_index) = self.iter.next() {
            Some(
                self.grid
                   .reference_to_physical_map(self.reference_points, cell_index),
            )
        } else {
            None
        }
    }
}
