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
        // Create the topology
        let topology = SerialFlatTriangleTopology::new(cells);

        // Create the geometry
        let geometry = SerialFlatTriangleGeometry::<T>::new(points, cells);

        Self { topology, geometry }
    }
}

/// A grid
impl<T: Float + Scalar<Real = T> + Inverse> Grid for SerialFlatTriangleGrid<T> {
    type T = T;

    /// The type that implements [Topology]
    type Topology = SerialFlatTriangleTopology;

    /// The type that implements [Geometry]
    type Geometry = SerialFlatTriangleGeometry<T>;

    /// Get the grid topology (See [Topology])
    fn topology(&self) -> &Self::Topology {
        &self.topology
    }

    /// Get the grid geometry (See [Geometry])
    fn geometry(&self) -> &Self::Geometry {
        &self.geometry
    }

    // Check if the grid is stored in serial
    fn is_serial(&self) -> bool {
        true
    }
}
