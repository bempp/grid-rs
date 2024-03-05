//! Implementation of grid geometry

use crate::grid_impl::common::{compute_jacobian, compute_normal_from_jacobian23, compute_point};
use crate::grid_impl::traits::{Geometry, GeometryEvaluator};
use crate::reference_cell;
use crate::types::ReferenceCellType;
use bempp_element::element::CiarletElement;
use bempp_traits::element::FiniteElement;
use num::Float;
use rlst_common::types::Scalar;
use rlst_dense::{
    array::Array,
    base_array::BaseArray,
    data_container::VectorContainer,
    rlst_dynamic_array4,
    traits::{RandomAccessByRef, Shape},
};

/// Geometry of a serial grid
pub struct SerialMixedGeometry<T: Float + Scalar> {
    dim: usize,
    index_map: Vec<(usize, usize)>,
    coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
    cells: Vec<Vec<usize>>,
    pub(crate) elements: Vec<CiarletElement<T>>,
    midpoints: Vec<Vec<Vec<T>>>,
    diameters: Vec<Vec<T>>,
    volumes: Vec<Vec<T>>,
}

unsafe impl<T: Float + Scalar> Sync for SerialMixedGeometry<T> {}

impl<T: Float + Scalar> SerialMixedGeometry<T> {
    pub fn new(
        coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
        cells_input: &[usize],
        elements: Vec<CiarletElement<T>>,
        cell_elements: &[usize],
    ) -> Self {
        let dim = coordinates.shape()[1];
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
    type Evaluator<'a> = GeometryEvaluatorMixed<'a, T>;

    fn dim(&self) -> usize {
        self.dim
    }

    fn index_map(&self) -> &[(usize, usize)] {
        &self.index_map
    }

    fn coordinate(&self, point_index: usize, coord_index: usize) -> Option<&Self::T> {
        self.coordinates.get([point_index, coord_index])
    }

    fn point_count(&self) -> usize {
        self.coordinates.shape()[0]
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

    fn get_evaluator<'a, Points: RandomAccessByRef<2, Item = T> + Shape<2>>(
        &'a self,
        points: &Points,
    ) -> Self::Evaluator<'a> {
        GeometryEvaluatorMixed::new(self, points)
    }
}

pub struct GeometryEvaluatorMixed<'a, T: Float + Scalar> {
    geometry: &'a SerialMixedGeometry<T>,
    tdim: usize,
    tables: Vec<Array<T, BaseArray<T, VectorContainer<T>, 4>, 4>>,
}

impl<'a, T: Float + Scalar> GeometryEvaluatorMixed<'a, T> {
    fn new<Points: RandomAccessByRef<2, Item = T> + Shape<2>>(
        geometry: &'a SerialMixedGeometry<T>,
        points: &Points,
    ) -> Self {
        let tdim = reference_cell::dim(
            // TODO: remove this match once bempp-rs and grid-rs use the same ReferenceCellType
            match geometry.elements[0].cell_type() {
                bempp_element::cell::ReferenceCellType::Interval => ReferenceCellType::Interval,
                bempp_element::cell::ReferenceCellType::Triangle => ReferenceCellType::Triangle,
                bempp_element::cell::ReferenceCellType::Quadrilateral => {
                    ReferenceCellType::Quadrilateral
                }
                _ => {
                    panic!(
                        "Unsupported cell type: {:?}",
                        geometry.elements[0].cell_type()
                    );
                }
            },
        );
        assert_eq!(tdim, points.shape()[1]);
        let mut tables = vec![];
        for e in &geometry.elements {
            assert_eq!(
                reference_cell::dim(
                    // TODO: remove this match once bempp-rs and grid-rs use the same ReferenceCellType
                    match e.cell_type() {
                        bempp_element::cell::ReferenceCellType::Interval =>
                            ReferenceCellType::Interval,
                        bempp_element::cell::ReferenceCellType::Triangle =>
                            ReferenceCellType::Triangle,
                        bempp_element::cell::ReferenceCellType::Quadrilateral => {
                            ReferenceCellType::Quadrilateral
                        }
                        _ => {
                            panic!("Unsupported cell type: {:?}", e.cell_type());
                        }
                    },
                ),
                tdim
            );
            let mut table = rlst_dynamic_array4!(T, e.tabulate_array_shape(1, points.shape()[0]));
            e.tabulate(points, 1, &mut table);
            tables.push(table);
        }
        Self {
            geometry,
            tdim,
            tables,
        }
    }
}

impl<'a, T: Float + Scalar> GeometryEvaluator for GeometryEvaluatorMixed<'a, T> {
    type T = T;

    fn point_count(&self) -> usize {
        self.tables[0].shape()[1]
    }

    fn compute_point(&self, cell_index: usize, point_index: usize, point: &mut [T]) {
        let cell = self.geometry.index_map()[cell_index];
        compute_point(
            self.geometry,
            self.tables[cell.0].view(),
            cell_index,
            point_index,
            point,
        );
    }

    fn compute_jacobian(&self, cell_index: usize, point_index: usize, jacobian: &mut [T]) {
        let cell = self.geometry.index_map()[cell_index];
        compute_jacobian(
            self.geometry,
            self.tables[cell.0].view(),
            self.tdim,
            cell_index,
            point_index,
            jacobian,
        );
    }

    fn compute_normal(&self, cell_index: usize, point_index: usize, normal: &mut [T]) {
        let gdim = self.geometry.dim();
        let tdim = self.tdim;
        assert_eq!(tdim, 2);
        assert_eq!(tdim, gdim - 1);

        // TODO: remove memory assignment?
        let mut jacobian = vec![T::from(0.0).unwrap(); gdim * tdim];
        self.compute_jacobian(cell_index, point_index, &mut jacobian[..]);
        compute_normal_from_jacobian23(&jacobian, normal);
    }
}

#[cfg(test)]
mod test {
    use crate::grid_impl::mixed_grid::geometry::*;
    use approx::*;
    use bempp_element::element::{create_element, ElementFamily};
    use bempp_traits::element::Continuity;
    use rlst_dense::{rlst_dynamic_array2, traits::RandomAccessMut};

    fn example_geometry() -> SerialMixedGeometry<f64> {
        let p1triangle = create_element(
            ElementFamily::Lagrange,
            bempp_element::cell::ReferenceCellType::Triangle,
            1,
            Continuity::Continuous,
        );
        let mut points = rlst_dynamic_array2!(f64, [4, 2]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([2, 0]).unwrap() = 1.0;
        *points.get_mut([2, 1]).unwrap() = 1.0;
        *points.get_mut([3, 0]).unwrap() = 0.0;
        *points.get_mut([3, 1]).unwrap() = 1.0;
        SerialMixedGeometry::new(points, &[0, 1, 2, 0, 2, 3], vec![p1triangle], &[0, 0])
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
        let mut points = rlst_dynamic_array2!(f64, [5, 2]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([2, 0]).unwrap() = 0.0;
        *points.get_mut([2, 1]).unwrap() = 1.0;
        *points.get_mut([3, 0]).unwrap() = 1.0;
        *points.get_mut([3, 1]).unwrap() = 1.0;
        *points.get_mut([4, 0]).unwrap() = 2.0;
        *points.get_mut([4, 1]).unwrap() = 0.0;
        SerialMixedGeometry::new(
            points,
            &[0, 1, 2, 3, 1, 4, 3],
            vec![p1quad, p1triangle],
            &[0, 1],
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
                    assert_relative_eq!(
                        *coord,
                        *g.coordinate(vs[p_i], c_i).unwrap(),
                        epsilon = 1e-12
                    );
                }
            }
        }
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
                    assert_relative_eq!(
                        *coord,
                        *g.coordinate(vs[p_i], c_i).unwrap(),
                        epsilon = 1e-12
                    );
                }
            }
        }
    }
}
