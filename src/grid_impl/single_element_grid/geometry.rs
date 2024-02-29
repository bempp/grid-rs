//! Implementation of grid geometry

use crate::grid_impl::traits::Geometry;
use bempp_element::element::CiarletElement;
use bempp_traits::element::FiniteElement;
use num::Float;

/// Geometry of a serial grid
pub struct SerialSingleElementGeometry<T: Float> {
    dim: usize,
    index_map: Vec<usize>,
    // TODO: change storage to rlst
    coordinates: Vec<T>,
    cells: Vec<usize>,
    element: CiarletElement,
    midpoints: Vec<Vec<T>>,
    diameters: Vec<T>,
    volumes: Vec<T>,
}

unsafe impl<T: Float> Sync for SerialSingleElementGeometry<T> {}

impl<T: Float> SerialSingleElementGeometry<T> {
    pub fn new(
        coordinates: Vec<T>,
        dim: usize,
        cells_input: &[usize],
        element: CiarletElement,
    ) -> Self {
        let size = element.dim();
        let ncells = cells_input.len() / size;

        let mut index_map = vec![0; ncells];
        let mut cells = vec![];
        let mut midpoints = vec![];
        let mut diameters = vec![];
        let mut volumes = vec![];

        for (cell_i, i) in index_map.iter_mut().enumerate() {
            *i = cell_i;
            midpoints.push(vec![T::from(0.0).unwrap(); dim]); // TODO
            diameters.push(T::from(0.0).unwrap()); // TODO
            volumes.push(T::from(0.0).unwrap()); // TODO
        }
        cells.extend_from_slice(cells_input);

        Self {
            dim,
            index_map,
            coordinates,
            cells,
            element,
            midpoints,
            diameters,
            volumes,
        }
    }
}

impl<T: Float> Geometry for SerialSingleElementGeometry<T> {
    type IndexType = usize;
    type T = T;
    type Element = CiarletElement;

    fn dim(&self) -> usize {
        self.dim
    }

    fn index_map(&self) -> &[usize] {
        &self.index_map
    }

    fn coordinate(&self, point_index: usize, coord_index: usize) -> Option<&Self::T> {
        if coord_index < self.dim && point_index * self.dim < self.coordinates.len() {
            Some(&self.coordinates[point_index * self.dim + coord_index])
        } else {
            None
        }
    }

    fn point_count(&self) -> usize {
        self.coordinates.len() / self.dim
    }

    fn cell_points(&self, index: usize) -> Option<&[usize]> {
        let npts = self.element.dim();
        if index * npts < self.cells.len() {
            Some(&self.cells[npts * index..npts * (index + 1)])
        } else {
            None
        }
    }

    fn cell_count(&self) -> usize {
        self.cells.len() / self.element.dim()
    }

    fn cell_element(&self, index: usize) -> Option<&Self::Element> {
        if index < self.cells.len() {
            Some(&self.element)
        } else {
            None
        }
    }

    fn element_count(&self) -> usize {
        1
    }
    fn element(&self, i: usize) -> Option<&Self::Element> {
        if i < self.cell_count() {
            Some(&self.element)
        } else {
            None
        }
    }
    fn cells(&self, i: usize) -> Option<&[usize]> {
        if i == 0 {
            Some(&self.cells)
        } else {
            None
        }
    }

    fn midpoint(&self, index: usize, point: &mut [Self::T]) {
        point.copy_from_slice(&self.midpoints[index]);
    }

    fn diameter(&self, index: usize) -> Self::T {
        self.diameters[index]
    }
    fn volume(&self, index: usize) -> Self::T {
        self.volumes[index]
    }
}

#[cfg(test)]
mod test {
    use crate::grid_impl::single_element_grid::geometry::*;
    use approx::*;
    use bempp_element::element::{create_element, ElementFamily};
    use bempp_traits::element::Continuity;

    fn example_geometry() -> SerialSingleElementGeometry<f64> {
        let p1triangle = create_element(
            ElementFamily::Lagrange,
            bempp_element::cell::ReferenceCellType::Triangle,
            1,
            Continuity::Continuous,
        );
        SerialSingleElementGeometry::new(
            vec![0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0],
            2,
            &[0, 1, 2, 0, 2, 3],
            p1triangle,
        )
    }

    #[test]
    fn test_counts() {
        let g = example_geometry();
        assert_eq!(g.point_count(), 4);
        assert_eq!(g.cell_count(), 2);
    }

    #[test]
    fn test_cell_points() {
        let g = example_geometry();
        for (cell_i, points) in [
            vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![1.0, 1.0]],
            vec![vec![0.0, 0.0], vec![1.0, 1.0], vec![0.0, 1.0]],
        ]
        .iter()
        .enumerate()
        {
            let vs = g.cell_points(cell_i).unwrap();
            for (p_i, point) in points.iter().enumerate() {
                for (c_i, coord) in point.iter().enumerate() {
                    assert_relative_eq!(*coord, *g.coordinate(vs[p_i], c_i).unwrap());
                }
            }
        }
    }
}
