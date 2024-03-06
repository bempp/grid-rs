//! Serial implementation of a grid

use crate::grid_impl::flat_triangle_grid::{
    geometry::SerialFlatTriangleGeometry, topology::SerialFlatTriangleTopology,
};
use crate::grid_impl::traits::Grid;
use crate::reference_cell;
use crate::reference_cell::ReferenceCellType;
use bempp_element::element::{create_element, ElementFamily, Inverse};
use bempp_traits::element::{Continuity, FiniteElement};
use log::warn;
use num::Float;
use rlst_common::types::Scalar;
use rlst_dense::{array::Array, base_array::BaseArray, data_container::VectorContainer};

/// A serial grid
pub struct SerialFlatTriangleGrid<T: Float + Scalar + Inverse> {
    topology: SerialFlatTriangleTopology,
    geometry: SerialFlatTriangleGeometry<T>,
}

impl<T: Float + Scalar + Inverse> SerialFlatTriangleGrid<T> {
    pub fn new(
        points: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
        cells: &[usize],
        cell_type: ReferenceCellType,
        cell_degree: usize,
    ) -> Self {
        panic!();
        /*        if cell_type == ReferenceCellType::Triangle && cell_degree == 1 {
                    warn!("Creating a single element grid with a P1 triangle. Using a FlatTriangleGrid would be faster.");
                }
                let element = create_element::<T>(
                    ElementFamily::Lagrange,
                    cell_type,
                    cell_degree,
                    Continuity::Continuous,
                );

                let mut cell_vertices = vec![];

                let mut start = 0;
                let nvertices = reference_cell::entity_counts(cell_type)[0];
                let npoints = element.dim();
                while start < cells.len() {
                    cell_vertices.extend_from_slice(&cells[start..start + nvertices]);
                    start += npoints;
                }

                // Create the topology
                let topology = SerialFlatTriangleTopology::new(&cell_vertices, cell_type);

                // Create the geometry
                let geometry = SerialFlatTriangleGeometry::<T>::new(points, cells, element);

                Self { topology, geometry }
        */
    }
}

/// A grid
impl<T: Float + Scalar + Inverse> Grid for SerialFlatTriangleGrid<T> {
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
