//! Serial implementation of a grid

use crate::grid_impl::single_element_grid::{
    geometry::SerialSingleElementGeometry, topology::SerialSingleElementTopology,
};
use crate::grid_impl::traits::Grid;
use crate::reference_cell;
use crate::reference_cell::ReferenceCellType;
use bempp_element::element::{create_element, ElementFamily, Inverse};
use bempp_traits::element::{Continuity, FiniteElement};
use num::Float;
use rlst_common::types::Scalar;
use rlst_dense::{array::Array, base_array::BaseArray, data_container::VectorContainer};

/// A serial grid
pub struct SerialSingleElementGrid<T: Float + Scalar> {
    topology: SerialSingleElementTopology,
    geometry: SerialSingleElementGeometry<T>,
}

impl<T: Float + Scalar + Inverse> SerialSingleElementGrid<T> {
    pub fn new(
        points: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
        cells: &[usize],
        cell_type: ReferenceCellType,
        cell_degree: usize,
    ) -> Self {
        let element = create_element::<T>(
            ElementFamily::Lagrange,
            cell_type,
            cell_degree,
            Continuity::Continuous,
        );

        let mut cell_vertices = vec![];

        let mut start = 0;
        let nvertices = reference_cell::entity_counts(cell_type)[0];
        println!("{nvertices}");
        let npoints = element.dim();
        while start < cells.len() {
            cell_vertices.extend_from_slice(&cells[start..start + nvertices]);
            start += npoints;
        }

        // Create the topology
        let topology = SerialSingleElementTopology::new(&cell_vertices, cell_type);

        // Create the geometry
        let geometry = SerialSingleElementGeometry::<T>::new(points, cells, element);

        Self { topology, geometry }
    }
}

/// A grid
impl<T: Float + Scalar + Inverse> Grid for SerialSingleElementGrid<T> {
    type T = T;

    /// The type that implements [Topology]
    type Topology = SerialSingleElementTopology;

    /// The type that implements [Geometry]
    type Geometry = SerialSingleElementGeometry<T>;

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
