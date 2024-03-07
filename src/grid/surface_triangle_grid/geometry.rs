//! Definition of the Geometry

use std::iter::Copied;

use num::Float;
use rlst_common::types::Scalar;
use rlst_dense::traits::DefaultIterator;

use crate::{
    traits::{GeometryType, GridType},
    types::point_iterator::PointIterator,
};

use super::grid::TriangleSurfaceGrid;

pub struct TriangleGeometry<'a, T: Float + Scalar> {
    cell_index: usize,
    grid: &'a TriangleSurfaceGrid<T>,
}

impl<'a, T: Float + Scalar> TriangleGeometry<'a, T> {
    pub fn new(grid: &'a TriangleSurfaceGrid<T>, cell_index: usize) -> Self {
        Self { cell_index, grid }
    }
}

impl<'a, T: Float + Scalar> GeometryType for TriangleGeometry<'a, T> {
    type Grid = TriangleSurfaceGrid<T>;

    type VertexIterator<'iter> =
        PointIterator<'iter, Self::Grid, Copied<std::slice::Iter<'iter, usize>>> where Self: 'iter;

    type PointIterator<'iter> = Self::VertexIterator<'iter> where Self: 'iter;

    fn physical_dimension(&self) -> usize {
        3
    }

    fn midpoint(&self, point: &mut [<Self::Grid as GridType>::T]) {
        for (ind_out, ind_in) in point
            .iter_mut()
            .zip(self.grid.midpoints[self.cell_index].iter())
        {
            *ind_out = ind_in
        }
    }

    fn diameter(&self) -> <Self::Grid as GridType>::T {
        self.grid.diameters[self.cell_index]
    }

    fn volume(&self) -> <Self::Grid as GridType>::T {
        self.grid.volumes[self.cell_index]
    }

    fn vertices(&self) -> Self::VertexIterator<'_> {
        self.grid
            .iter_points(self.grid.cells[self.cell_index].iter().copied())
    }

    fn points(&self) -> Self::PointIterator<'_> {
        self.vertices()
    }
}
