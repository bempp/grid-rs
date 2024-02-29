//! Serial implementation of a grid

use crate::grid_impl::single_element_grid::{
    geometry::SerialSingleElementGeometry, topology::SerialSingleElementTopology,
};
use crate::grid_impl::traits::Grid;
use crate::reference_cell;
use crate::reference_cell::ReferenceCellType;
use bempp_element::element::{create_element, ElementFamily};
use bempp_traits::element::{Continuity, FiniteElement};
use num::Float;

/// A serial grid
pub struct SerialSingleElementGrid<T: Float> {
    topology: SerialSingleElementTopology,
    geometry: SerialSingleElementGeometry<T>,
}

impl<T: Float> SerialSingleElementGrid<T> {
    pub fn new(
        points: Vec<T>,
        gdim: usize,
        cells: &[usize],
        cell_type: ReferenceCellType,
        cell_degree: usize,
    ) -> Self {
        let element = create_element(
            ElementFamily::Lagrange,
            // TODO: remove this match once bempp-rs and grid-rs use the same ReferenceCellType
            match cell_type {
                ReferenceCellType::Interval => bempp_element::cell::ReferenceCellType::Interval,
                ReferenceCellType::Triangle => bempp_element::cell::ReferenceCellType::Triangle,
                ReferenceCellType::Quadrilateral => {
                    bempp_element::cell::ReferenceCellType::Quadrilateral
                }
                _ => {
                    panic!("Unsupported cell type: {:?}", cell_type);
                }
            },
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
        let topology = SerialSingleElementTopology::new(&cell_vertices, cell_type);

        // Create the geometry
        let geometry = SerialSingleElementGeometry::<T>::new(points, gdim, cells, element);

        Self { topology, geometry }
    }
}

/// A grid
impl<T: Float> Grid for SerialSingleElementGrid<T> {
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