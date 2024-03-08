//! Implementation of grid geometry

use crate::grid::common::{
    compute_diameter_quadrilateral, compute_diameter_triangle, compute_jacobian,
    compute_normal_from_jacobian23, compute_point,
};
use crate::grid::traits::{Geometry, GeometryEvaluator};
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
    rlst_array_from_slice2, rlst_dynamic_array1, rlst_dynamic_array4,
    traits::{DefaultIteratorMut, RandomAccessByRef, Shape, UnsafeRandomAccessByRef},
};
use std::collections::HashMap;

/// Geometry of a serial grid
pub struct SerialMixedGeometry<T: Float + Scalar<Real = T>> {
    dim: usize,
    index_map: Vec<(usize, usize)>,
    coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
    cells: Vec<Vec<usize>>,
    pub(crate) elements: Vec<CiarletElement<T>>,
    midpoints: Vec<Vec<Vec<T>>>,
    diameters: Vec<Vec<T>>,
    volumes: Vec<Vec<T>>,
    cell_indices: Vec<Vec<(usize, usize)>>,
    point_indices_to_ids: Vec<usize>,
    point_ids_to_indices: HashMap<usize, usize>,
    cell_indices_to_ids: HashMap<(usize, usize), usize>,
    cell_ids_to_indices: HashMap<usize, (usize, usize)>,
}

unsafe impl<T: Float + Scalar<Real = T>> Sync for SerialMixedGeometry<T> {}

impl<T: Float + Scalar<Real = T>> SerialMixedGeometry<T> {
    /// Create a geometry
    pub fn new(
        coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
        cells_input: &[usize],
        elements: Vec<CiarletElement<T>>,
        cell_elements: &[usize],
        point_indices_to_ids: Vec<usize>,
        point_ids_to_indices: HashMap<usize, usize>,
        grid_cell_indices_to_ids: &[usize],
    ) -> Self {
        let dim = coordinates.shape()[1];
        let mut index_map = vec![(0, 0); cell_elements.len()];
        let mut cells = vec![];
        let mut cell_indices = vec![];

        let mut cell_indices_to_ids = HashMap::new();
        let mut cell_ids_to_indices = HashMap::new();

        for (element_index, _e) in elements.iter().enumerate() {
            let mut start = 0;

            cells.push(vec![]);
            cell_indices.push(vec![]);

            for (cell_i, element_i) in cell_elements.iter().enumerate() {
                let size = elements[*element_i].dim();
                if *element_i == element_index {
                    let cell_index = (element_index, cell_indices[element_index].len());
                    index_map[cell_i] = cell_index;
                    cell_indices[element_index].push(cell_index);
                    cells[element_index].extend_from_slice(&cells_input[start..start + size]);
                    cell_indices_to_ids.insert(cell_index, grid_cell_indices_to_ids[cell_i]);
                    cell_ids_to_indices.insert(grid_cell_indices_to_ids[cell_i], cell_index);
                }
                start += size;
            }
        }

        let mut midpoints = vec![vec![]; cell_elements.len()];
        let mut diameters = vec![vec![]; cell_elements.len()];
        let mut volumes = vec![vec![]; cell_elements.len()];

        for (element_index, e) in elements.iter().enumerate() {
            let ncells = cells[element_index].len() / e.dim();
            let size = e.dim();

            midpoints[element_index] = vec![vec![T::from(0.0).unwrap(); dim]; ncells];
            diameters[element_index] = vec![T::from(0.0).unwrap(); ncells];
            volumes[element_index] = vec![T::from(0.0).unwrap(); ncells]; // TODO

            let mut table = rlst_dynamic_array4!(T, e.tabulate_array_shape(0, 1));
            e.tabulate(
                &rlst_array_from_slice2!(
                    T,
                    &reference_cell::midpoint(e.cell_type()),
                    [1, reference_cell::dim(e.cell_type())]
                ),
                0,
                &mut table,
            );

            let mut start = 0;
            for cell_i in 0..ncells {
                for (i, v) in cells[element_index][start..start + size].iter().enumerate() {
                    let t = unsafe { *table.get_unchecked([0, 0, i, 0]) };
                    for (j, component) in midpoints[element_index][cell_i].iter_mut().enumerate() {
                        *component += unsafe { *coordinates.get_unchecked([*v, j]) } * t;
                    }
                }

                start += size;
            }

            match e.cell_type() {
                ReferenceCellType::Triangle => {
                    let mut v0 = rlst_dynamic_array1!(T, [dim]);
                    let mut v1 = rlst_dynamic_array1!(T, [dim]);
                    let mut v2 = rlst_dynamic_array1!(T, [dim]);
                    for cell_i in 0..ncells {
                        for (j, v) in [&mut v0, &mut v1, &mut v2].iter_mut().enumerate() {
                            for (i, c) in v.iter_mut().enumerate() {
                                *c = unsafe {
                                    *coordinates
                                        .get_unchecked([cells[element_index][size * cell_i + j], i])
                                };
                            }
                        }
                        diameters[element_index][cell_i] =
                            compute_diameter_triangle(v0.view(), v1.view(), v2.view());
                    }
                }
                ReferenceCellType::Quadrilateral => {
                    let mut v0 = rlst_dynamic_array1!(T, [dim]);
                    let mut v1 = rlst_dynamic_array1!(T, [dim]);
                    let mut v2 = rlst_dynamic_array1!(T, [dim]);
                    let mut v3 = rlst_dynamic_array1!(T, [dim]);
                    for cell_i in 0..ncells {
                        for (j, v) in [&mut v0, &mut v1, &mut v2, &mut v3].iter_mut().enumerate() {
                            for (i, c) in v.iter_mut().enumerate() {
                                *c = unsafe {
                                    *coordinates
                                        .get_unchecked([cells[element_index][size * cell_i + j], i])
                                };
                            }
                        }
                        diameters[element_index][cell_i] = compute_diameter_quadrilateral(
                            v0.view(),
                            v1.view(),
                            v2.view(),
                            v3.view(),
                        );
                    }
                }
                _ => {
                    panic!("Unsupported cell type: {:?}", e.cell_type());
                }
            }
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
            cell_indices,
            point_indices_to_ids,
            point_ids_to_indices,
            cell_indices_to_ids,
            cell_ids_to_indices,
        }
    }
}

impl<T: Float + Scalar<Real = T>> Geometry for SerialMixedGeometry<T> {
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
    fn cell_indices(&self, i: usize) -> Option<&[Self::IndexType]> {
        if i < self.cells.len() {
            Some(&self.cell_indices[i])
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

    fn get_evaluator<'a>(&'a self, points: &'a [Self::T]) -> Self::Evaluator<'a> {
        GeometryEvaluatorMixed::new(self, points)
    }

    fn point_index_to_id(&self, index: usize) -> usize {
        self.point_indices_to_ids[index]
    }
    fn cell_index_to_id(&self, index: (usize, usize)) -> usize {
        self.cell_indices_to_ids[&index]
    }
    fn point_id_to_index(&self, id: usize) -> usize {
        self.point_ids_to_indices[&id]
    }
    fn cell_id_to_index(&self, id: usize) -> (usize, usize) {
        self.cell_ids_to_indices[&id]
    }
}

/// Geometry evaluator for a mixed grid
pub struct GeometryEvaluatorMixed<'a, T: Float + Scalar<Real = T>> {
    geometry: &'a SerialMixedGeometry<T>,
    tdim: usize,
    tables: Vec<Array<T, BaseArray<T, VectorContainer<T>, 4>, 4>>,
}

impl<'a, T: Float + Scalar<Real = T>> GeometryEvaluatorMixed<'a, T> {
    /// Create a geometry evaluator
    fn new(geometry: &'a SerialMixedGeometry<T>, points: &'a [T]) -> Self {
        let tdim = reference_cell::dim(geometry.elements[0].cell_type());
        assert_eq!(points.len() % tdim, 0);
        let npoints = points.len() / tdim;
        let rlst_points = rlst_array_from_slice2!(T, points, [tdim, npoints]);

        let mut tables = vec![];
        for e in &geometry.elements {
            assert_eq!(reference_cell::dim(e.cell_type()), tdim);
            let mut table = rlst_dynamic_array4!(T, e.tabulate_array_shape(1, npoints));
            e.tabulate(&rlst_points, 1, &mut table);
            tables.push(table);
        }
        Self {
            geometry,
            tdim,
            tables,
        }
    }
}

impl<'a, T: Float + Scalar<Real = T>> GeometryEvaluator for GeometryEvaluatorMixed<'a, T> {
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

        let mut jacobian = vec![T::from(0.0).unwrap(); gdim * tdim];
        self.compute_jacobian(cell_index, point_index, &mut jacobian[..]);
        compute_normal_from_jacobian23(&jacobian, normal);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::*;
    use bempp_element::element::{create_element, ElementFamily};
    use bempp_traits::element::Continuity;
    use rlst_dense::{rlst_dynamic_array2, traits::RandomAccessMut};

    fn example_geometry() -> SerialMixedGeometry<f64> {
        //! A geometry with a single cell type
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
        SerialMixedGeometry::new(
            points,
            &[0, 1, 2, 0, 2, 3],
            vec![p1triangle],
            &[0, 0],
            vec![0, 1, 2, 3],
            HashMap::from([(0, 0), (1, 1), (2, 2), (3, 3)]),
            &[0, 1],
        )
    }

    fn example_geometry_mixed() -> SerialMixedGeometry<f64> {
        //! A geometry with a mixture of cell types
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
            vec![0, 1, 2, 3, 4],
            HashMap::from([(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]),
            &[0, 1],
        )
    }

    #[test]
    fn test_counts() {
        //! Test the point and cell counts
        let g = example_geometry();
        assert_eq!(g.point_count(), 4);
        assert_eq!(g.cell_count(), 2);
    }

    #[test]
    fn test_cell_points() {
        //! Test the cell points
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
        //! Test the point and cell counts
        let g = example_geometry_mixed();
        assert_eq!(g.point_count(), 5);
        assert_eq!(g.cell_count(), 2);
    }

    #[test]
    fn test_cell_points_mixed() {
        //! Test the cell points
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

    #[test]
    fn test_midpoint() {
        //! Test midpoints
        let g = example_geometry();

        let mut midpoint = vec![0.0; 2];
        for (cell_i, point) in [vec![2.0 / 3.0, 1.0 / 3.0], vec![1.0 / 3.0, 2.0 / 3.0]]
            .iter()
            .enumerate()
        {
            g.midpoint(g.index_map()[cell_i], &mut midpoint);
            for (i, j) in midpoint.iter().zip(point) {
                assert_relative_eq!(*i, *j, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn test_midpoint_mixed() {
        //! Test midpoints
        let g = example_geometry_mixed();

        let mut midpoint = vec![0.0; 2];
        for (cell_i, point) in [vec![0.5, 0.5], vec![4.0 / 3.0, 1.0 / 3.0]]
            .iter()
            .enumerate()
        {
            g.midpoint(g.index_map()[cell_i], &mut midpoint);
            for (i, j) in midpoint.iter().zip(point) {
                assert_relative_eq!(*i, *j, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn test_diameter() {
        //! Test diameters
        let g = example_geometry();

        for cell_i in 0..2 {
            assert_relative_eq!(
                g.diameter((0, cell_i)),
                2.0 * f64::sqrt(1.5 - f64::sqrt(2.0)),
                epsilon = 1e-12
            );
        }

        let g = example_geometry_mixed();

        for (cell_i, d) in [f64::sqrt(2.0), 2.0 * f64::sqrt(1.5 - f64::sqrt(2.0))]
            .iter()
            .enumerate()
        {
            assert_relative_eq!(g.diameter((cell_i, 0)), d, epsilon = 1e-12);
        }
    }
}
