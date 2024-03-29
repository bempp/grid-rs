//! Implementation of grid geometry

use crate::grid::common::{
    compute_det, compute_diameter_quadrilateral, compute_diameter_triangle, compute_jacobian,
    compute_normal_from_jacobian23, compute_point,
};
use crate::grid::traits::{Geometry, GeometryEvaluator};
use crate::reference_cell;
use crate::types::ReferenceCellType;
use bempp_element::element::CiarletElement;
use bempp_quadrature::simplex_rules::simplex_rule;
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

/// Geometry of a single element grid
pub struct SerialSingleElementGeometry<T: Float + Scalar> {
    dim: usize,
    index_map: Vec<usize>,
    pub(crate) coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
    pub(crate) cells: Vec<usize>,
    pub(crate) element: CiarletElement<T>,
    midpoints: Vec<Vec<T>>,
    diameters: Vec<T>,
    volumes: Vec<T>,
    cell_indices: Vec<usize>,
    point_indices_to_ids: Vec<usize>,
    point_ids_to_indices: HashMap<usize, usize>,
    cell_indices_to_ids: Vec<usize>,
    cell_ids_to_indices: HashMap<usize, usize>,
}

unsafe impl<T: Float + Scalar<Real = T>> Sync for SerialSingleElementGeometry<T> {}

impl<T: Float + Scalar<Real = T>> SerialSingleElementGeometry<T> {
    /// Create a geometry
    pub fn new(
        coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
        cells_input: &[usize],
        element: CiarletElement<T>,
        point_indices_to_ids: Vec<usize>,
        point_ids_to_indices: HashMap<usize, usize>,
        cell_indices_to_ids: Vec<usize>,
        cell_ids_to_indices: HashMap<usize, usize>,
    ) -> Self {
        let dim = coordinates.shape()[1];
        let tdim = reference_cell::dim(element.cell_type());
        let size = element.dim();
        let ncells = cells_input.len() / size;

        let mut index_map = vec![0; ncells];
        let mut cells = vec![];
        let mut midpoints = vec![vec![T::from(0.0).unwrap(); dim]; ncells];
        let mut diameters = vec![T::from(0.0).unwrap(); ncells];
        let mut volumes = vec![T::from(0.0).unwrap(); ncells];

        let mut mpt_table = rlst_dynamic_array4!(T, element.tabulate_array_shape(0, 1));
        element.tabulate(
            &rlst_array_from_slice2!(T, &reference_cell::midpoint(element.cell_type()), [1, tdim]),
            0,
            &mut mpt_table,
        );

        // TODO: pick rule number of points sensibly
        // NOTE: 37 used for now as rules with 37 points exist for both a triangle and a quadrilateral
        let nqpts = 37;
        let qrule = simplex_rule(element.cell_type(), nqpts).unwrap();
        let qpoints = qrule
            .points
            .iter()
            .map(|x| T::from(*x).unwrap())
            .collect::<Vec<_>>();
        let qweights = qrule
            .weights
            .iter()
            .map(|x| T::from(*x).unwrap())
            .collect::<Vec<_>>();

        let mut jdet_table = rlst_dynamic_array4!(T, element.tabulate_array_shape(1, nqpts));
        element.tabulate(
            &rlst_array_from_slice2!(T, &qpoints, [nqpts, tdim], [tdim, 1]),
            1,
            &mut jdet_table,
        );

        let mut jacobian = vec![T::from(0.0).unwrap(); dim * tdim];

        for (cell_i, index) in index_map.iter_mut().enumerate() {
            *index = cell_i;

            for (i, v) in cells_input[size * cell_i..size * (cell_i + 1)]
                .iter()
                .enumerate()
            {
                let t = unsafe { *mpt_table.get_unchecked([0, 0, i, 0]) };
                for (j, component) in midpoints[cell_i].iter_mut().enumerate() {
                    *component += unsafe { *coordinates.get_unchecked([*v, j]) } * t;
                }
            }
            for (point_index, w) in qweights.iter().enumerate() {
                for component in jacobian.iter_mut() {
                    *component = T::from(0.0).unwrap();
                }
                for (i, v) in cells_input[size * cell_i..size * (cell_i + 1)]
                    .iter()
                    .enumerate()
                {
                    for gd in 0..dim {
                        for td in 0..tdim {
                            jacobian[td * dim + gd] +=
                                unsafe { *coordinates.get_unchecked([*v, gd]) }
                                    * unsafe {
                                        *jdet_table.get_unchecked([1 + td, point_index, i, 0])
                                    };
                        }
                    }
                }
                volumes[cell_i] += *w * compute_det(&jacobian, tdim, dim);
            }
        }
        cells.extend_from_slice(cells_input);

        match element.cell_type() {
            ReferenceCellType::Triangle => {
                let mut v0 = rlst_dynamic_array1!(T, [dim]);
                let mut v1 = rlst_dynamic_array1!(T, [dim]);
                let mut v2 = rlst_dynamic_array1!(T, [dim]);
                for cell_i in 0..ncells {
                    for (j, v) in [&mut v0, &mut v1, &mut v2].iter_mut().enumerate() {
                        for (i, c) in v.iter_mut().enumerate() {
                            *c = unsafe {
                                *coordinates.get_unchecked([cells[size * cell_i + j], i])
                            };
                        }
                    }
                    diameters[cell_i] = compute_diameter_triangle(v0.view(), v1.view(), v2.view());
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
                                *coordinates.get_unchecked([cells[size * cell_i + j], i])
                            };
                        }
                    }
                    diameters[cell_i] =
                        compute_diameter_quadrilateral(v0.view(), v1.view(), v2.view(), v3.view());
                }
            }
            _ => {
                panic!("Unsupported cell type: {:?}", element.cell_type());
            }
        }

        let cell_indices = (0..ncells).collect::<Vec<_>>();

        Self {
            dim,
            index_map,
            coordinates,
            cells,
            element,
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

impl<T: Float + Scalar> Geometry for SerialSingleElementGeometry<T> {
    type IndexType = usize;
    type T = T;
    type Element = CiarletElement<T>;
    type Evaluator<'a> = GeometryEvaluatorSingleElement<'a, T> where Self: 'a;

    fn dim(&self) -> usize {
        self.dim
    }

    fn index_map(&self) -> &[usize] {
        &self.index_map
    }

    fn coordinate(&self, point_index: usize, coord_index: usize) -> Option<&Self::T> {
        self.coordinates.get([point_index, coord_index])
    }

    fn point_count(&self) -> usize {
        self.coordinates.shape()[0]
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
        if i == 0 {
            Some(&self.element)
        } else {
            None
        }
    }
    fn cell_indices(&self, i: usize) -> Option<&[usize]> {
        if i == 0 {
            Some(&self.cell_indices)
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

    fn get_evaluator<'a>(&'a self, points: &'a [Self::T]) -> GeometryEvaluatorSingleElement<'a, T> {
        GeometryEvaluatorSingleElement::<T>::new(self, points)
    }

    fn point_index_to_id(&self, index: usize) -> usize {
        self.point_indices_to_ids[index]
    }
    fn cell_index_to_id(&self, index: usize) -> usize {
        self.cell_indices_to_ids[index]
    }
    fn point_id_to_index(&self, id: usize) -> usize {
        self.point_ids_to_indices[&id]
    }
    fn cell_id_to_index(&self, id: usize) -> usize {
        self.cell_ids_to_indices[&id]
    }
}

/// Geometry evaluator for a single element grid
pub struct GeometryEvaluatorSingleElement<'a, T: Float + Scalar> {
    geometry: &'a SerialSingleElementGeometry<T>,
    tdim: usize,
    table: Array<T, BaseArray<T, VectorContainer<T>, 4>, 4>,
}

impl<'a, T: Float + Scalar> GeometryEvaluatorSingleElement<'a, T> {
    /// Create a geometry evaluator
    fn new(geometry: &'a SerialSingleElementGeometry<T>, points: &'a [T]) -> Self {
        let tdim = reference_cell::dim(geometry.element.cell_type());
        assert_eq!(points.len() % tdim, 0);
        let npoints = points.len() / tdim;
        let rlst_points = rlst_array_from_slice2!(T, points, [tdim, npoints]);

        let mut table = rlst_dynamic_array4!(T, geometry.element.tabulate_array_shape(1, npoints));
        geometry.element.tabulate(&rlst_points, 1, &mut table);
        Self {
            geometry,
            tdim,
            table,
        }
    }
}

impl<'a, T: Float + Scalar> GeometryEvaluator for GeometryEvaluatorSingleElement<'a, T> {
    type T = T;

    fn point_count(&self) -> usize {
        self.table.shape()[1]
    }

    fn compute_point(&self, cell_index: usize, point_index: usize, point: &mut [T]) {
        compute_point(
            self.geometry,
            self.table.view(),
            cell_index,
            point_index,
            point,
        );
    }

    fn compute_jacobian(&self, cell_index: usize, point_index: usize, jacobian: &mut [T]) {
        compute_jacobian(
            self.geometry,
            self.table.view(),
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
    use crate::types::ReferenceCellType;
    use approx::*;
    use bempp_element::element::{create_element, ElementFamily};
    use bempp_traits::element::Continuity;
    use rlst_dense::{
        rlst_dynamic_array2,
        traits::{RandomAccessMut, RawAccess, RawAccessMut},
    };

    fn example_geometry_2d() -> SerialSingleElementGeometry<f64> {
        //! A 2D geometry
        let p1triangle = create_element(
            ElementFamily::Lagrange,
            ReferenceCellType::Triangle,
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
        SerialSingleElementGeometry::new(
            points,
            &[0, 1, 2, 0, 2, 3],
            p1triangle,
            vec![0, 1, 2, 3],
            HashMap::from([(0, 0), (1, 1), (2, 2), (3, 3)]),
            vec![0, 1],
            HashMap::from([(0, 0), (1, 1)]),
        )
    }

    fn example_geometry_3d() -> SerialSingleElementGeometry<f64> {
        //! A 3D geometry
        let p2triangle = create_element(
            ElementFamily::Lagrange,
            ReferenceCellType::Triangle,
            2,
            Continuity::Continuous,
        );
        let mut points = rlst_dynamic_array2!(f64, [9, 3]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.5;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([1, 2]).unwrap() = 0.5;
        *points.get_mut([2, 0]).unwrap() = 1.0;
        *points.get_mut([2, 1]).unwrap() = 0.0;
        *points.get_mut([2, 2]).unwrap() = 0.0;
        *points.get_mut([3, 0]).unwrap() = 0.0;
        *points.get_mut([3, 1]).unwrap() = 0.5;
        *points.get_mut([3, 2]).unwrap() = 0.0;
        *points.get_mut([4, 0]).unwrap() = 0.5;
        *points.get_mut([4, 1]).unwrap() = 0.5;
        *points.get_mut([4, 2]).unwrap() = 0.0;
        *points.get_mut([5, 0]).unwrap() = 1.0;
        *points.get_mut([5, 1]).unwrap() = 0.5;
        *points.get_mut([5, 2]).unwrap() = 0.0;
        *points.get_mut([6, 0]).unwrap() = 0.0;
        *points.get_mut([6, 1]).unwrap() = 1.0;
        *points.get_mut([6, 2]).unwrap() = 0.0;
        *points.get_mut([7, 0]).unwrap() = 0.5;
        *points.get_mut([7, 1]).unwrap() = 1.0;
        *points.get_mut([7, 2]).unwrap() = 0.0;
        *points.get_mut([8, 0]).unwrap() = 1.0;
        *points.get_mut([8, 1]).unwrap() = 1.0;
        *points.get_mut([8, 2]).unwrap() = 0.0;
        SerialSingleElementGeometry::new(
            points,
            &[0, 2, 8, 5, 4, 1, 0, 8, 6, 7, 3, 4],
            p2triangle,
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
            HashMap::from([
                (0, 0),
                (1, 1),
                (2, 2),
                (3, 3),
                (4, 4),
                (5, 5),
                (6, 6),
                (7, 7),
                (8, 8),
            ]),
            vec![0, 1],
            HashMap::from([(0, 0), (1, 1)]),
        )
    }

    fn example_geometry_quad() -> SerialSingleElementGeometry<f64> {
        //! A 3D quadrilateral geometry
        let p1quad = create_element(
            ElementFamily::Lagrange,
            ReferenceCellType::Quadrilateral,
            1,
            Continuity::Continuous,
        );
        let mut points = rlst_dynamic_array2!(f64, [6, 3]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([1, 2]).unwrap() = 0.0;
        *points.get_mut([2, 0]).unwrap() = 2.0;
        *points.get_mut([2, 1]).unwrap() = 0.0;
        *points.get_mut([2, 2]).unwrap() = 1.0;
        *points.get_mut([3, 0]).unwrap() = 0.0;
        *points.get_mut([3, 1]).unwrap() = 1.0;
        *points.get_mut([3, 2]).unwrap() = 0.0;
        *points.get_mut([4, 0]).unwrap() = 1.0;
        *points.get_mut([4, 1]).unwrap() = 1.0;
        *points.get_mut([4, 2]).unwrap() = 0.0;
        *points.get_mut([5, 0]).unwrap() = 2.0;
        *points.get_mut([5, 1]).unwrap() = 1.0;
        *points.get_mut([5, 2]).unwrap() = 1.0;
        SerialSingleElementGeometry::new(
            points,
            &[0, 1, 3, 4, 1, 2, 4, 5],
            p1quad,
            vec![0, 1, 2, 3, 4, 5],
            HashMap::from([(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]),
            vec![0, 1],
            HashMap::from([(0, 0), (1, 1)]),
        )
    }

    fn triangle_points() -> Array<f64, BaseArray<f64, VectorContainer<f64>, 2>, 2> {
        //! Create a set of points in the reference triangle
        let mut points = rlst_dynamic_array2!(f64, [2, 2]);
        *points.get_mut([0, 0]).unwrap() = 0.2;
        *points.get_mut([0, 1]).unwrap() = 0.5;
        *points.get_mut([1, 0]).unwrap() = 0.6;
        *points.get_mut([1, 1]).unwrap() = 0.1;
        points
    }

    #[test]
    fn test_counts() {
        //! Test the point and cell counts
        let g = example_geometry_2d();
        assert_eq!(g.point_count(), 4);
        assert_eq!(g.cell_count(), 2);
    }

    #[test]
    fn test_cell_points() {
        //! Test the cell points
        let g = example_geometry_2d();
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
    fn test_compute_point_2d() {
        //! Test the compute_point function of an evaluator
        let g = example_geometry_2d();
        let points = triangle_points();

        let evaluator = g.get_evaluator(points.data());
        let mut mapped_point = vec![0.0; 2];
        for (cell_i, points) in [
            vec![vec![0.7, 0.5], vec![0.7, 0.1]],
            vec![vec![0.2, 0.7], vec![0.6, 0.7]],
        ]
        .iter()
        .enumerate()
        {
            for (point_i, point) in points.iter().enumerate() {
                evaluator.compute_point(cell_i, point_i, &mut mapped_point);
                for (i, j) in mapped_point.iter().zip(point) {
                    assert_relative_eq!(*i, *j, epsilon = 1e-12);
                }
            }
        }
    }

    #[test]
    fn test_compute_point_3d() {
        //! Test the compute_point function of an evaluator
        let g = example_geometry_3d();
        let points = triangle_points();
        let evaluator = g.get_evaluator(points.data());

        let mut mapped_point = vec![0.0; 3];
        for (cell_i, points) in [
            vec![vec![0.7, 0.5, 0.12], vec![0.7, 0.1, 0.36]],
            vec![vec![0.2, 0.7, 0.0], vec![0.6, 0.7, 0.0]],
        ]
        .iter()
        .enumerate()
        {
            for (point_i, point) in points.iter().enumerate() {
                evaluator.compute_point(cell_i, point_i, &mut mapped_point);
                for (i, j) in mapped_point.iter().zip(point) {
                    assert_relative_eq!(*i, *j, epsilon = 1e-12);
                }
            }
        }
    }

    #[test]
    fn test_compute_jacobian_3d() {
        //! Test the compute_jacobian function of an evaluator
        let g = example_geometry_3d();
        let points = triangle_points();
        let evaluator = g.get_evaluator(points.data());

        let mut computed_jacobian = rlst_dynamic_array2!(f64, [3, 2]);
        for (cell_i, jacobians) in [
            vec![
                vec![vec![1.0, 1.0], vec![0.0, 1.0], vec![0.2, -0.4]],
                vec![vec![1.0, 1.0], vec![0.0, 1.0], vec![-0.6, -1.2]],
            ],
            vec![
                vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![0.0, 0.0]],
                vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![0.0, 0.0]],
            ],
        ]
        .iter()
        .enumerate()
        {
            for (point_i, jacobian) in jacobians.iter().enumerate() {
                evaluator.compute_jacobian(cell_i, point_i, computed_jacobian.data_mut());
                for (i, row) in jacobian.iter().enumerate() {
                    for (j, entry) in row.iter().enumerate() {
                        assert_relative_eq!(
                            *entry,
                            *computed_jacobian.get([i, j]).unwrap(),
                            epsilon = 1e-12
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_compute_normal_3d() {
        //! Test the compute_normal function of an evaluator
        let g = example_geometry_3d();
        let points = triangle_points();
        let evaluator = g.get_evaluator(points.data());

        let mut computed_normal = vec![0.0; 3];
        for (cell_i, normals) in [
            vec![
                vec![
                    -0.2 / f64::sqrt(1.4),
                    0.6 / f64::sqrt(1.4),
                    1.0 / f64::sqrt(1.4),
                ],
                vec![
                    0.6 / f64::sqrt(1.72),
                    0.6 / f64::sqrt(1.72),
                    1.0 / f64::sqrt(1.72),
                ],
            ],
            vec![vec![0.0, 0.0, 1.0], vec![0.0, 0.0, 1.0]],
        ]
        .iter()
        .enumerate()
        {
            for (point_i, normal) in normals.iter().enumerate() {
                evaluator.compute_normal(cell_i, point_i, &mut computed_normal);
                assert_relative_eq!(
                    computed_normal[0] * computed_normal[0]
                        + computed_normal[1] * computed_normal[1]
                        + computed_normal[2] * computed_normal[2],
                    1.0,
                    epsilon = 1e-12
                );
                for (i, j) in computed_normal.iter().zip(normal) {
                    assert_relative_eq!(*i, *j, epsilon = 1e-12);
                }
            }
        }
    }

    #[test]
    fn test_midpoint_2d() {
        //! Test midpoints
        let g = example_geometry_2d();

        let mut midpoint = vec![0.0; 2];
        for (cell_i, point) in [vec![2.0 / 3.0, 1.0 / 3.0], vec![1.0 / 3.0, 2.0 / 3.0]]
            .iter()
            .enumerate()
        {
            g.midpoint(cell_i, &mut midpoint);
            for (i, j) in midpoint.iter().zip(point) {
                assert_relative_eq!(*i, *j, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn test_midpoint_3d() {
        //! Test midpoints
        let g = example_geometry_3d();

        let mut midpoint = vec![0.0; 3];
        for (cell_i, point) in [
            vec![2.0 / 3.0, 1.0 / 3.0, 2.0 / 9.0],
            vec![1.0 / 3.0, 2.0 / 3.0, 0.0],
        ]
        .iter()
        .enumerate()
        {
            g.midpoint(cell_i, &mut midpoint);
            for (i, j) in midpoint.iter().zip(point) {
                assert_relative_eq!(*i, *j, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn test_diameter() {
        //! Test diameters
        let g = example_geometry_2d();

        for cell_i in 0..2 {
            assert_relative_eq!(
                g.diameter(cell_i),
                2.0 * f64::sqrt(1.5 - f64::sqrt(2.0)),
                epsilon = 1e-12
            );
        }

        let g = example_geometry_3d();

        for cell_i in 0..2 {
            assert_relative_eq!(
                g.diameter(cell_i),
                2.0 * f64::sqrt(1.5 - f64::sqrt(2.0)),
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn test_volume() {
        //! Test cell volumes
        let g = example_geometry_2d();

        for cell_i in 0..2 {
            assert_relative_eq!(g.volume(cell_i), 0.5, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_volume_3d() {
        //! Test cell volumes
        let g = example_geometry_3d();

        for (cell_i, d) in [0.7390096708393067, 0.5].iter().enumerate() {
            assert_relative_eq!(g.volume(cell_i), d, epsilon = 1e-5);
        }
    }

    #[test]
    fn test_volume_quad() {
        //! Test cell volumes
        let g = example_geometry_quad();

        for (cell_i, d) in [1.0, f64::sqrt(2.0)].iter().enumerate() {
            assert_relative_eq!(g.volume(cell_i), d, epsilon = 1e-5);
        }
    }
}
