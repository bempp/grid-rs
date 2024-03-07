//! Serial implementation of a grid

use crate::grid::mixed_grid::{geometry::SerialMixedGeometry, topology::SerialMixedTopology};
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
pub struct SerialMixedGrid<T: Float + Scalar> {
    topology: SerialMixedTopology,
    geometry: SerialMixedGeometry<T>,
}

impl<T: Float + Scalar + Inverse> SerialMixedGrid<T> {
    pub fn new(
        points: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
        cells: &[usize],
        cell_types: &[ReferenceCellType],
        cell_degrees: &[usize],
    ) -> Self {
        let mut element_info = vec![];
        let mut element_numbers = vec![];

        for (cell, degree) in cell_types.iter().zip(cell_degrees) {
            let info = (*cell, *degree);
            if !element_info.contains(&info) {
                element_info.push(info);
            }
            element_numbers.push(element_info.iter().position(|&i| i == info).unwrap());
        }

        let elements = element_info
            .iter()
            .map(|(i, j)| {
                create_element::<T>(ElementFamily::Lagrange, *i, *j, Continuity::Continuous)
            })
            .collect::<Vec<_>>();

        if elements.len() == 1 {
            warn!("Creating a mixed grid with only one element. Using a SerialSingleElementGrid would be faster.");
        }

        let mut cell_vertices = vec![];

        let mut start = 0;
        for (cell_type, e_n) in cell_types.iter().zip(&element_numbers) {
            let nvertices = reference_cell::entity_counts(*cell_type)[0];
            let npoints = elements[*e_n].dim();
            cell_vertices.extend_from_slice(&cells[start..start + nvertices]);
            start += npoints;
        }

        let topology = SerialMixedTopology::new(&cell_vertices, cell_types);
        let geometry = SerialMixedGeometry::<T>::new(points, cells, elements, &element_numbers);

        Self { topology, geometry }
    }
}

impl<T: Float + Scalar + Inverse> Grid for SerialMixedGrid<T> {
    type T = T;
    type Topology = SerialMixedTopology;
    type Geometry = SerialMixedGeometry<T>;

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
