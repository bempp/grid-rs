//! Serial implementation of a grid

use crate::grid::flat_triangle_grid::{
    geometry::SerialFlatTriangleGeometry, topology::SerialFlatTriangleTopology,
};
use crate::grid::traits::Grid;
use bempp_element::element::Inverse;
use num::Float;
use rlst_common::types::Scalar;
use rlst_dense::{array::Array, base_array::BaseArray, data_container::VectorContainer};

/// A serial grid
pub struct SerialFlatTriangleGrid<T: Float + Scalar<Real = T> + Inverse> {
    topology: SerialFlatTriangleTopology,
    geometry: SerialFlatTriangleGeometry<T>,
}

impl<T: Float + Scalar<Real = T> + Inverse> SerialFlatTriangleGrid<T> {
    pub fn new(points: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>, cells: &[usize]) -> Self {
        // TODO: cell and vertex ids
        let topology = SerialFlatTriangleTopology::new(cells);
        let geometry = SerialFlatTriangleGeometry::<T>::new(points, cells);
        Self { topology, geometry }
    }
}

impl<T: Float + Scalar<Real = T> + Inverse> Grid for SerialFlatTriangleGrid<T> {
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
