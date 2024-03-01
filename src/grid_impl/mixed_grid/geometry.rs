//! Implementation of grid geometry

use crate::grid_impl::traits::Geometry;
use bempp_element::element::CiarletElement;
use bempp_traits::element::FiniteElement;
use num::Float;
use rlst_common::types::Scalar;

/// Geometry of a serial grid
pub struct SerialMixedGeometry<T: Float + Scalar> {
    dim: usize,
    index_map: Vec<(usize, usize)>,
    // TODO: change storage to rlst
    coordinates: Vec<T>,
    cells: Vec<Vec<usize>>,
    elements: Vec<CiarletElement<T>>,
    midpoints: Vec<Vec<Vec<T>>>,
    diameters: Vec<Vec<T>>,
    volumes: Vec<Vec<T>>,
}

unsafe impl<T: Float + Scalar> Sync for SerialMixedGeometry<T> {}

impl<T: Float + Scalar> SerialMixedGeometry<T> {
    pub fn new(
        coordinates: Vec<T>,
        dim: usize,
        cells_input: &[usize],
        elements: Vec<CiarletElement<T>>,
        cell_elements: &[usize],
    ) -> Self {
        let mut index_map = vec![(0, 0); cell_elements.len()];
        let mut cells = vec![];
        let mut midpoints = vec![];
        let mut diameters = vec![];
        let mut volumes = vec![];

        for (element_index, _e) in elements.iter().enumerate() {
            let mut e_cells = vec![];
            let mut start = 0;

            for (cell_i, element_i) in cell_elements.iter().enumerate() {
                let size = elements[*element_i].dim();
                if *element_i == element_index {
                    index_map[cell_i] = (element_index, e_cells.len() / size);
                    e_cells.extend_from_slice(&cells_input[start..start + size]);
                }
                start += size;
            }
            cells.push(e_cells);
        }

        for (element_index, e) in elements.iter().enumerate() {
            let mut e_midpoints = vec![];
            let mut e_diameters = vec![];
            let mut e_volumes = vec![];
            for _cell in 0..cells[element_index].len() / e.dim() {
                e_midpoints.push(vec![T::from(0.0).unwrap(); dim]); // TODO
                e_diameters.push(T::from(0.0).unwrap()); // TODO
                e_volumes.push(T::from(0.0).unwrap()); // TODO
            }
            midpoints.push(e_midpoints);
            diameters.push(e_diameters);
            volumes.push(e_volumes);
        }

        Self {
            dim,
            index_map,
            coordinates,
            cells,
            elements,
            midpoints,
            diameters,
            volumes,
        }
    }
}

impl<T: Float + Scalar> Geometry for SerialMixedGeometry<T> {
    type IndexType = (usize, usize);
    type T = T;
    type Element = CiarletElement<T>;

    fn dim(&self) -> usize {
        self.dim
    }

    fn index_map(&self) -> &[(usize, usize)] {
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

    fn cell_points(&self, index: (usize, usize)) -> Option<&[usize]> {
        if index.0 < self.cells.len() {
            let npts = self.elements[index.0].dim();
            if index.1 * npts < self.cells[index.0].len() {
                Some(&self.cells[index.0][npts * index.1..npts * (index.1 + 1)])
            } else {
                None
            }
        } else {
            None
        }
    }

    fn cell_count(&self) -> usize {
        self.elements
            .iter()
            .enumerate()
            .map(|(i, e)| self.cells[i].len() / e.dim())
            .sum()
    }

    fn cell_element(&self, index: (usize, usize)) -> Option<&Self::Element> {
        if index.0 < self.cells.len() {
            Some(&self.elements[index.0])
        } else {
            None
        }
    }

    fn element_count(&self) -> usize {
        self.elements.len()
    }
    fn element(&self, i: usize) -> Option<&Self::Element> {
        if i < self.elements.len() {
            Some(&self.elements[i])
        } else {
            None
        }
    }
    fn cells(&self, i: usize) -> Option<&[usize]> {
        if i < self.cells.len() {
            Some(&self.cells[i])
        } else {
            None
        }
    }

    fn midpoint(&self, index: (usize, usize), point: &mut [Self::T]) {
        point.copy_from_slice(&self.midpoints[index.0][index.1]);
    }

    fn diameter(&self, index: (usize, usize)) -> Self::T {
        self.diameters[index.0][index.1]
    }
    fn volume(&self, index: (usize, usize)) -> Self::T {
        self.volumes[index.0][index.1]
    }
}

#[cfg(test)]
mod test {
    use crate::grid_impl::mixed_grid::geometry::*;
    use approx::*;
    use bempp_element::element::{create_element, ElementFamily};
    use bempp_traits::element::Continuity;

    fn example_geometry() -> SerialMixedGeometry<f64> {
        let p1triangle = create_element(
            ElementFamily::Lagrange,
            bempp_element::cell::ReferenceCellType::Triangle,
            1,
            Continuity::Continuous,
        );
        SerialMixedGeometry::new(
            vec![0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0],
            2,
            &[0, 1, 2, 0, 2, 3],
            vec![p1triangle],
            &[0, 0],
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
            let vs = g.cell_points((0, cell_i)).unwrap();
            for (p_i, point) in points.iter().enumerate() {
                for (c_i, coord) in point.iter().enumerate() {
                    assert_relative_eq!(*coord, *g.coordinate(vs[p_i], c_i).unwrap());
                }
            }
        }
    }

    fn example_geometry_mixed() -> SerialMixedGeometry<f64> {
        let p1triangle = create_element(
            ElementFamily::Lagrange,
            bempp_element::cell::ReferenceCellType::Triangle,
            1,
            Continuity::Continuous,
        );
        let p1quad = create_element(
            ElementFamily::Lagrange,
            bempp_element::cell::ReferenceCellType::Quadrilateral,
            1,
            Continuity::Continuous,
        );
        SerialMixedGeometry::new(
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 0.0],
            2,
            &[0, 1, 2, 3, 1, 4, 3],
            vec![p1quad, p1triangle],
            &[0, 1],
        )
    }

    #[test]
    fn test_counts_mixed() {
        let g = example_geometry_mixed();
        assert_eq!(g.point_count(), 5);
        assert_eq!(g.cell_count(), 2);
    }

    #[test]
    fn test_cell_points_mixed() {
        let g = example_geometry_mixed();
        for (cell_i, points) in [
            (
                (0, 0),
                vec![
                    vec![0.0, 0.0],
                    vec![1.0, 0.0],
                    vec![0.0, 1.0],
                    vec![1.0, 1.0],
                ],
            ),
            ((1, 0), vec![vec![1.0, 0.0], vec![2.0, 0.0], vec![1.0, 1.0]]),
        ] {
            let vs = g.cell_points(cell_i).unwrap();
            for (p_i, point) in points.iter().enumerate() {
                for (c_i, coord) in point.iter().enumerate() {
                    assert_relative_eq!(*coord, *g.coordinate(vs[p_i], c_i).unwrap());
                }
            }
        }
    }
}
