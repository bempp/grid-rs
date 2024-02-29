//! Serial implementation of a grid

use crate::grid_impl::traits::Grid;
use crate::grid_impl::{geometry::SerialGeometry, topology::SerialTopology};
use crate::reference_cell::ReferenceCellType;
use bempp_element::element::{create_element, ElementFamily};
use bempp_traits::element::Continuity;
use num::Float;

/// A serial grid
pub struct SerialGrid<T: Float> {
    topology: SerialTopology,
    geometry: SerialGeometry<T>,
}

impl<T: Float> SerialGrid<T> {
    pub fn new(
        _vertices: Vec<T>,
        _cells: &[usize],
        _cell_types: &[ReferenceCellType],
        _cell_degrees: &[usize],
    ) -> Self {
        // TODO
        let topology = SerialTopology::new(&[0, 1, 2, 2, 1, 3], &[ReferenceCellType::Triangle; 2]);
        let p1triangle = create_element(
            ElementFamily::Lagrange,
            bempp_element::cell::ReferenceCellType::Triangle,
            1,
            Continuity::Continuous,
        );
        let geometry = SerialGeometry::<T>::new(
            [0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0]
                .iter()
                .map(|x| T::from(*x).unwrap())
                .collect::<Vec<T>>(),
            2,
            &[0, 1, 2, 0, 2, 3],
            vec![p1triangle],
            &[0, 0],
        );

        Self { topology, geometry }
    }
}

/// A grid
impl<T: Float> Grid for SerialGrid<T> {
    type T = T;

    /// The type that implements [Topology]
    type Topology = SerialTopology;

    /// The type that implements [Geometry]
    type Geometry = SerialGeometry<T>;

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
