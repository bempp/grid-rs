//! Triangle topology

use std::iter::Copied;
use std::marker::PhantomData;

use crate::types::Float;
use rlst_common::types::Scalar;

use crate::grid::surface_triangle_grid::cell::TriangleCell;
use crate::traits::*;

use super::grid::TriangleSurfaceGrid;

pub struct TriangleTopology<'a, T: Float + Scalar> {
    cell: &'a TriangleCell<'a, T>,
    _marker: PhantomData<T>,
}

impl<'a, T: Float + Scalar> TriangleTopology<'a, T> {
    pub fn new(cell: &'a TriangleCell<'a, T>) -> Self {
        Self {
            cell,
            _marker: PhantomData,
        }
    }
}

impl<'a, T: Float + Scalar> TopologyType for TriangleTopology<'a, T> {
    type Grid = TriangleSurfaceGrid<T>;

    type VertexIndexIter<'iter> = Copied<std::slice::Iter<'iter, usize>> where Self: 'iter;

    fn vertex_indices(&self) -> Self::VertexIndexIter<'_> {
        self.cell.grid().cells[self.cell.index()].iter().copied()
    }
}
