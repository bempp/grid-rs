//! Implementation of grid geometry

use crate::grid::traits::Geometry;
use num::Float;
//use crate::grid::traits::{Ownership, Geometry};
//use crate::reference_cell;
//use crate::traits::cell::ReferenceCellType;
use bempp_element::element::CiarletElement;
use bempp_traits::element::FiniteElement;

//use std::collections::HashMap;

/// Geometry of a serial grid
pub struct SerialGeometry<T: Float> {
    dim: usize,
    index_map: Vec<usize>,
    coordinates: Vec<T>,
    cells: Vec<Vec<usize>>,
    elements: Vec<CiarletElement>,
}

unsafe impl<T: Float> Sync for SerialGeometry<T> {}

impl<T: Float> SerialGeometry<T> {
    pub fn new(
        coordinates: Vec<T>,
        dim: usize,
        cells_input: &[usize],
        elements: Vec<CiarletElement>,
        cell_elements: &[usize],
    ) -> Self {
        let mut index_map = vec![];
        //let mut element_changes = vec![];
        //let mut coordinate_elements = vec![];
        let mut cells = vec![];

        for (element_index, _e) in elements.iter().enumerate() {
            let mut e_cells = vec![];
            let mut start = 0;

            for (cell_i, element_i) in cell_elements.iter().enumerate() {
                let size = elements[*element_i].dim();
                println!("{element_i} {size}");
                if *element_i == element_index {
                    index_map.push(cell_i);
                    e_cells.extend_from_slice(&cells_input[start..start + size])
                }
                start += size;
            }
            cells.push(e_cells);
        }

        println!("{cells:?}");

        Self {
            dim,
            index_map,
            coordinates,
            cells,
            elements,
        }
    }
}

impl<T: Float> Geometry for SerialGeometry<T> {
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

    fn cell_vertices(&self, index: usize) -> Option<&[usize]> {
        let mut start = 0;
        for (element_index, e) in self.elements.iter().enumerate() {
            let npts = e.dim();
            let ncells = self.cells[element_index].len() / npts;
            if index < start + ncells {
                return Some(
                    &self.cells[element_index][npts * (index - start)..npts * (index - start + 1)],
                );
            }
            start += ncells;
        }
        None
    }

    fn cell_count(&self) -> usize {
        self.elements
            .iter()
            .enumerate()
            .map(|(i, e)| self.cells[i].len() / e.dim())
            .sum()
    }

    fn cell_element(&self, index: usize) -> Option<&Self::Element> {
        let mut start = 0;
        for (element_index, e) in self.elements.iter().enumerate() {
            let npts = e.dim();
            let ncells = self.cells[element_index].len() / npts;
            if index < start + ncells {
                return Some(e);
            }
            start += ncells;
        }
        None
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
}

#[cfg(test)]
mod test {
    use crate::grid::geometry::*;
    use approx::*;
    use bempp_element::element::{create_element, ElementFamily};
    use bempp_traits::element::Continuity;

    fn example_geometry() -> SerialGeometry<f64> {
        let p1triangle = create_element(
            ElementFamily::Lagrange,
            bempp_element::cell::ReferenceCellType::Triangle,
            1,
            Continuity::Continuous,
        );
        SerialGeometry::new(
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
    fn test_cell_vertices() {
        let g = example_geometry();
        for (cell_i, vertices) in [
            vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![1.0, 1.0]],
            vec![vec![0.0, 0.0], vec![1.0, 1.0], vec![0.0, 1.0]],
        ]
        .iter()
        .enumerate()
        {
            let vs = g.cell_vertices(cell_i).unwrap();
            for (p_i, point) in vertices.iter().enumerate() {
                for (c_i, coord) in point.iter().enumerate() {
                    assert_relative_eq!(*coord, *g.coordinate(vs[p_i], c_i).unwrap());
                }
            }
        }
    }

    fn example_geometry_mixed() -> SerialGeometry<f64> {
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
        SerialGeometry::new(
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
    fn test_cell_vertices_mixed() {
        let g = example_geometry_mixed();
        println!("{:?}", g.cell_vertices(0).unwrap());
        println!("{:?}", g.cell_vertices(1).unwrap());
        for (cell_i, vertices) in [
            vec![
                vec![0.0, 0.0],
                vec![1.0, 0.0],
                vec![0.0, 1.0],
                vec![1.0, 1.0],
            ],
            vec![vec![1.0, 0.0], vec![2.0, 0.0], vec![1.0, 1.0]],
        ]
        .iter()
        .enumerate()
        {
            let vs = g.cell_vertices(cell_i).unwrap();
            println!("{vs:?}");
            for (p_i, point) in vertices.iter().enumerate() {
                for (c_i, coord) in point.iter().enumerate() {
                    assert_relative_eq!(*coord, *g.coordinate(vs[p_i], c_i).unwrap());
                }
            }
        }
    }
}
