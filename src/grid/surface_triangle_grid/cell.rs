//! Definition of a triangle cell.

use crate::types::Float;
use rlst_common::types::Scalar;

use crate::traits::*;

use super::{geometry::TriangleGeometry, grid::TriangleSurfaceGrid, topology::TriangleTopology};

pub struct TriangleCell<'a, T: Float + Scalar> {
    index: usize,
    grid: &'a TriangleSurfaceGrid<T>,
}

impl<'a, T: Float + Scalar> TriangleCell<'a, T> {
    pub fn new(index: usize, grid: &'a TriangleSurfaceGrid<T>) -> Self {
        Self { index, grid }
    }
}

impl<'a, T: Float + Scalar> CellType for TriangleCell<'a, T> {
    type Grid = TriangleSurfaceGrid<T>;
    type Topology<'top> = TriangleTopology<'top, T> where Self: 'top;

    type Geometry<'geom>  = TriangleGeometry<'geom, T>
    where
        Self: 'geom;

    fn id(&self) -> usize {
        self.grid.cell_id_from_index(self.index)
    }

    fn index(&self) -> usize {
        self.index
    }

    fn topology(&self) -> Self::Topology<'_> {
        TriangleTopology::new(self)
    }

    fn grid(&self) -> &Self::Grid {
        self.grid
    }

    fn geometry(&self) -> Self::Geometry<'_> {
        TriangleGeometry::new(self.grid, self.index)
    }
}
