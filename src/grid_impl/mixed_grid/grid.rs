//! Serial implementation of a grid

use crate::grid_impl::mixed_grid::{geometry::SerialMixedGeometry, topology::SerialMixedTopology};
use crate::grid_impl::traits::Grid;
use crate::reference_cell;
use crate::reference_cell::ReferenceCellType;
use bempp_element::element::{create_element, ElementFamily, Inverse};
use bempp_traits::element::{Continuity, FiniteElement};
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
                create_element::<T>(
                    ElementFamily::Lagrange,
                    // TODO: remove this match once bempp-rs and grid-rs use the same ReferenceCellType
                    match i {
                        ReferenceCellType::Interval => {
                            bempp_element::cell::ReferenceCellType::Interval
                        }
                        ReferenceCellType::Triangle => {
                            bempp_element::cell::ReferenceCellType::Triangle
                        }
                        ReferenceCellType::Quadrilateral => {
                            bempp_element::cell::ReferenceCellType::Quadrilateral
                        }
                        _ => {
                            panic!("Unsupported cell type: {:?}", i);
                        }
                    },
                    *j,
                    Continuity::Continuous,
                )
            })
            .collect::<Vec<_>>();

        if elements.len() == 1 {
            println!("Creating a mixed grid with only one element. Using a SerialSingleElementGrid would be faster.");
            // TODO: make into a warning
        }

        let mut cell_vertices = vec![];

        let mut start = 0;
        for (cell_type, e_n) in cell_types.iter().zip(&element_numbers) {
            let nvertices = reference_cell::entity_counts(*cell_type)[0];
            let npoints = elements[*e_n].dim();
            cell_vertices.extend_from_slice(&cells[start..start + nvertices]);
            start += npoints;
        }

        // Create the topology
        let topology = SerialMixedTopology::new(&cell_vertices, cell_types);

        // Create the geometry
        let geometry = SerialMixedGeometry::<T>::new(points, cells, elements, &element_numbers);

        Self { topology, geometry }
    }
}

/// A grid
impl<T: Float + Scalar + Inverse> Grid for SerialMixedGrid<T> {
    type T = T;

    /// The type that implements [Topology]
    type Topology = SerialMixedTopology;

    /// The type that implements [Geometry]
    type Geometry = SerialMixedGeometry<T>;

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
