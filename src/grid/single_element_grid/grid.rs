//! Serial implementation of a grid

use crate::grid::single_element_grid::{
    geometry::SerialSingleElementGeometry, topology::SerialSingleElementTopology,
};
use crate::grid::traits::Grid;
use crate::reference_cell;
use crate::reference_cell::ReferenceCellType;
use bempp_element::element::{create_element, ElementFamily, Inverse};
use bempp_traits::element::{Continuity, FiniteElement};
use log::warn;
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
        if cell_type == ReferenceCellType::Triangle && cell_degree == 1 {
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

        let topology = SerialSingleElementTopology::new(&cell_vertices, cell_type);
        let geometry = SerialSingleElementGeometry::<T>::new(points, cells, element);
        Self { topology, geometry }
    }
}

impl<T: Float + Scalar + Inverse> Grid for SerialSingleElementGrid<T> {
    type T = T;
    type Topology = SerialSingleElementTopology;
    type Geometry = SerialSingleElementGeometry<T>;

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