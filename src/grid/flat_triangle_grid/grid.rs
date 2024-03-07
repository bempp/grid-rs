//! Serial implementation of a grid

use crate::grid::flat_triangle_grid::{
    geometry::SerialFlatTriangleGeometry, topology::SerialFlatTriangleTopology,
};
use crate::grid::traits::Grid;
use num::Float;
use rlst_common::types::Scalar;
use rlst_dense::{
    array::{views::ArrayViewMut, Array},
    base_array::BaseArray,
    data_container::VectorContainer,
    traits::MatrixInverse,
};

/// A serial grid
pub struct SerialFlatTriangleGrid<T: Float + Scalar<Real = T>> {
    topology: SerialFlatTriangleTopology,
    geometry: SerialFlatTriangleGeometry<T>,
}

impl<T: Float + Scalar<Real = T>> SerialFlatTriangleGrid<T>
where
    for<'a> Array<T, ArrayViewMut<'a, T, BaseArray<T, VectorContainer<T>, 2>, 2>, 2>: MatrixInverse,
{
    pub fn new(points: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>, cells: &[usize]) -> Self {
        // TODO: cell and vertex ids
        let topology = SerialFlatTriangleTopology::new(cells);
        let geometry = SerialFlatTriangleGeometry::<T>::new(points, cells);
        Self { topology, geometry }
    }
}

impl<T: Float + Scalar<Real = T>> Grid for SerialFlatTriangleGrid<T> {
    type T = T;
    type Topology = SerialFlatTriangleTopology;
    type Geometry = SerialFlatTriangleGeometry<T>;

    fn topology(&self) -> &Self::Topology {
        &self.topology
    }

    fn geometry(&self) -> &Self::Geometry {
        &self.geometry
    }

    fn is_serial(&self) -> bool {
        true
    }
}
