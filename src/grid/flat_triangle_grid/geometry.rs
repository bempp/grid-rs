//! Implementation of grid geometry

use crate::grid::traits::{Geometry, GeometryEvaluator};
use crate::reference_cell;
use crate::types::ReferenceCellType;
use bempp_element::element::{create_element, CiarletElement, ElementFamily};
use bempp_traits::element::{Continuity, FiniteElement};
use num::Float;
use rlst::rlst_static_array;
use rlst_common::types::Scalar;
use rlst_dense::{
    array::{views::ArrayViewMut, Array, SliceArray},
    base_array::BaseArray,
    data_container::VectorContainer,
    rlst_array_from_slice2,
    traits::{
        DefaultIterator, DefaultIteratorMut, MatrixInverse, RandomAccessByRef, RawAccess, Shape,
        UnsafeRandomAccessByRef,
    },
};
use rlst_proc_macro::rlst_static_type;

/// Geometry of a serial grid
pub struct SerialFlatTriangleGeometry<T: Float + Scalar> {
    dim: usize,
    index_map: Vec<usize>,
    pub(crate) coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
    pub(crate) cells: Vec<Vec<usize>>,
    pub(crate) element: CiarletElement<T>,
    midpoints: Vec<rlst_static_type!(T, 3)>,
    diameters: Vec<T>,
    volumes: Vec<T>,
    pub(crate) normals: Vec<rlst_static_type!(T, 3)>,
    pub(crate) jacobians: Vec<rlst_static_type!(T, 3, 2)>,
    cell_indices: Vec<usize>,
}

unsafe impl<T: Float + Scalar> Sync for SerialFlatTriangleGeometry<T> {}

impl<T: Float + Scalar<Real = T>> SerialFlatTriangleGeometry<T>
where
    for<'a> Array<T, ArrayViewMut<'a, T, BaseArray<T, VectorContainer<T>, 2>, 2>, 2>: MatrixInverse,
{
    pub fn new(
        coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
        cells_input: &[usize],
    ) -> Self {
        let dim = coordinates.shape()[1];
        assert_eq!(dim, 3);
        let ncells = cells_input.len() / dim;

        let mut cells = vec![vec![]; ncells];
        for (cell_i, cell) in cells.iter_mut().enumerate() {
            cell.extend_from_slice(&cells_input[cell_i * 3..(cell_i + 1) * 3])
        }

        let mut index_map = vec![0; ncells];
        let mut midpoints = vec![];
        let mut diameters = vec![T::from(0.0).unwrap(); ncells];
        let mut volumes = vec![T::from(0.0).unwrap(); ncells];
        let mut normals = vec![];
        let mut jacobians = vec![];

        let mut a = rlst_static_array!(T, 3);
        let mut b = rlst_static_array!(T, 3);
        let mut c = rlst_static_array!(T, 3);

        let mut v0 = rlst_static_array!(T, 3);
        let mut v1 = rlst_static_array!(T, 3);
        let mut v2 = rlst_static_array!(T, 3);

        for (cell_i, index) in index_map.iter_mut().enumerate() {
            *index = cell_i;

            midpoints.push(rlst_static_array!(T, 3));
            normals.push(rlst_static_array!(T, 3));
            jacobians.push(rlst_static_array!(T, 3, 2));

            for (i, c) in v0.iter_mut().enumerate() {
                *c = unsafe { *coordinates.get_unchecked([cells[cell_i][0], i]) };
            }
            for (i, c) in v1.iter_mut().enumerate() {
                *c = unsafe { *coordinates.get_unchecked([cells[cell_i][1], i]) };
            }
            for (i, c) in v2.iter_mut().enumerate() {
                *c = unsafe { *coordinates.get_unchecked([cells[cell_i][2], i]) };
            }

            midpoints[cell_i].fill_from(
                (v0.view() + v1.view() + v2.view()).scalar_mul(T::from(1.0 / 3.0).unwrap()),
            );

            a.fill_from(v1.view() - v0.view());
            b.fill_from(v2.view() - v0.view());
            c.fill_from(v2.view() - v1.view());
            jacobians[cell_i].view_mut().slice(1, 0).fill_from(a.view());
            jacobians[cell_i].view_mut().slice(1, 1).fill_from(b.view());

            let a_norm = a.view().norm_2();
            let b_norm = b.view().norm_2();
            let c_norm = c.view().norm_2();

            a.cross(b.view(), normals[cell_i].view_mut());

            let normal_length = normals[cell_i].view().norm_2();
            normals[cell_i].scale_in_place(T::one() / normal_length);

            volumes[cell_i] = normal_length / T::from(2.0).unwrap();

            let s = (a_norm + b_norm + c_norm) / T::from(2.0).unwrap();

            diameters[cell_i] = T::from(2.0).unwrap()
                * Scalar::sqrt(((s - a_norm) * (s - b_norm) * (s - c_norm)) / s);
        }

        let element = create_element(
            ElementFamily::Lagrange,
            ReferenceCellType::Triangle,
            1,
            Continuity::Continuous,
        );
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
            normals,
            jacobians,
            cell_indices,
        }
    }
}

impl<T: Float + Scalar> Geometry for SerialFlatTriangleGeometry<T> {
    type IndexType = usize;
    type T = T;
    type Element = CiarletElement<T>;
    type Evaluator<'a> = GeometryEvaluatorFlatTriangle<'a, T> where Self: 'a;

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
        if index < self.cells.len() {
            Some(&self.cells[index])
        } else {
            None
        }
    }

    fn cell_count(&self) -> usize {
        self.cells.len()
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
        point.copy_from_slice(self.midpoints[index].data());
    }

    fn diameter(&self, index: usize) -> Self::T {
        self.diameters[index]
    }
    fn volume(&self, index: usize) -> Self::T {
        self.volumes[index]
    }

    fn get_evaluator<'a>(&'a self, points: &'a [Self::T]) -> GeometryEvaluatorFlatTriangle<'a, T> {
        GeometryEvaluatorFlatTriangle::<T>::new(self, points)
    }
}

pub struct GeometryEvaluatorFlatTriangle<'a, T: Float + Scalar> {
    geometry: &'a SerialFlatTriangleGeometry<T>,
    points: SliceArray<'a, T, 2>,
}

impl<'a, T: Float + Scalar> GeometryEvaluatorFlatTriangle<'a, T> {
    fn new(geometry: &'a SerialFlatTriangleGeometry<T>, points: &'a [T]) -> Self {
        let tdim = reference_cell::dim(geometry.element.cell_type());
        assert_eq!(points.len() % tdim, 0);
        let npoints = points.len() / tdim;
        Self {
            geometry,
            points: rlst_array_from_slice2!(T, points, [npoints, tdim]),
        }
    }
}

impl<'a, T: Float + Scalar> GeometryEvaluator for GeometryEvaluatorFlatTriangle<'a, T> {
    type T = T;

    fn point_count(&self) -> usize {
        self.points.shape()[0]
    }

    fn compute_point(&self, cell_index: usize, point_index: usize, point: &mut [T]) {
        let jacobian = &self.geometry.jacobians[cell_index];
        for (index, val_out) in point.iter_mut().enumerate() {
            *val_out = self.geometry.coordinates[[self.geometry.cells[cell_index][0], index]]
                + jacobian[[index, 0]] * self.points[[point_index, 0]]
                + jacobian[[index, 1]] * self.points[[point_index, 1]];
        }
    }

    fn compute_jacobian(&self, cell_index: usize, _point_index: usize, jacobian: &mut [T]) {
        for (i, j) in jacobian
            .iter_mut()
            .zip(self.geometry.jacobians[cell_index].iter())
        {
            *i = j;
        }
    }

    fn compute_normal(&self, cell_index: usize, _point_index: usize, normal: &mut [T]) {
        for (i, j) in normal
            .iter_mut()
            .zip(self.geometry.normals[cell_index].iter())
        {
            *i = j;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::*;
    use rlst_dense::{
        rlst_dynamic_array2,
        traits::{RandomAccessMut, RawAccessMut},
    };

    fn example_geometry_flat() -> SerialFlatTriangleGeometry<f64> {
        let mut points = rlst_dynamic_array2!(f64, [4, 3]);
        points[[0, 0]] = 0.0;
        points[[0, 1]] = 0.0;
        points[[0, 2]] = 0.0;
        points[[1, 0]] = 1.0;
        points[[1, 1]] = 0.0;
        points[[1, 2]] = 0.0;
        points[[2, 0]] = 1.0;
        points[[2, 1]] = 1.0;
        points[[2, 2]] = 0.0;
        points[[3, 0]] = 0.0;
        points[[3, 1]] = 1.0;
        points[[3, 2]] = 0.0;
        let cells = vec![0, 1, 2, 0, 2, 3];
        SerialFlatTriangleGeometry::new(points, &cells)
    }

    fn example_geometry_3d() -> SerialFlatTriangleGeometry<f64> {
        let mut points = rlst_dynamic_array2!(f64, [4, 3]);
        points[[0, 0]] = 0.0;
        points[[0, 1]] = 0.0;
        points[[0, 2]] = 0.0;
        points[[1, 0]] = 1.0;
        points[[1, 1]] = 0.0;
        points[[1, 2]] = 1.0;
        points[[2, 0]] = 1.0;
        points[[2, 1]] = 1.0;
        points[[2, 2]] = 0.0;
        points[[3, 0]] = 0.0;
        points[[3, 1]] = 1.0;
        points[[3, 2]] = 0.0;
        let cells = vec![0, 1, 2, 0, 2, 3];
        SerialFlatTriangleGeometry::new(points, &cells)
    }

    #[test]
    fn test_counts() {
        let g = example_geometry_flat();
        assert_eq!(g.point_count(), 4);
        assert_eq!(g.cell_count(), 2);
    }

    #[test]
    fn test_cell_points() {
        let g = example_geometry_flat();
        for (cell_i, points) in [
            vec![
                vec![0.0, 0.0, 0.0],
                vec![1.0, 0.0, 0.0],
                vec![1.0, 1.0, 0.0],
            ],
            vec![
                vec![0.0, 0.0, 0.0],
                vec![1.0, 1.0, 0.0],
                vec![0.0, 1.0, 0.0],
            ],
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

    fn triangle_points() -> Array<f64, BaseArray<f64, VectorContainer<f64>, 2>, 2> {
        let mut points = rlst_dynamic_array2!(f64, [2, 2]);
        *points.get_mut([0, 0]).unwrap() = 0.2;
        *points.get_mut([0, 1]).unwrap() = 0.5;
        *points.get_mut([1, 0]).unwrap() = 0.6;
        *points.get_mut([1, 1]).unwrap() = 0.1;
        points
    }

    #[test]
    fn test_compute_point_flat() {
        let g = example_geometry_flat();
        let points = triangle_points();

        let evaluator = g.get_evaluator(points.data());
        let mut mapped_point = vec![0.0; 3];
        for (cell_i, pts) in [
            vec![vec![0.7, 0.5, 0.0], vec![0.7, 0.1, 0.0]],
            vec![vec![0.2, 0.7, 0.0], vec![0.6, 0.7, 0.0]],
        ]
        .iter()
        .enumerate()
        {
            for (point_i, point) in pts.iter().enumerate() {
                evaluator.compute_point(cell_i, point_i, &mut mapped_point);
                for (i, j) in mapped_point.iter().zip(point) {
                    assert_relative_eq!(*i, *j, epsilon = 1e-12);
                }
            }
        }
    }

    #[test]
    fn test_compute_point_3d() {
        let g = example_geometry_3d();
        let points = triangle_points();
        let evaluator = g.get_evaluator(points.data());

        let mut mapped_point = vec![0.0; 3];
        for (cell_i, pts) in [
            vec![vec![0.7, 0.5, 0.2], vec![0.7, 0.1, 0.6]],
            vec![vec![0.2, 0.7, 0.0], vec![0.6, 0.7, 0.0]],
        ]
        .iter()
        .enumerate()
        {
            for (point_i, point) in pts.iter().enumerate() {
                evaluator.compute_point(cell_i, point_i, &mut mapped_point);
                for (i, j) in mapped_point.iter().zip(point) {
                    assert_relative_eq!(*i, *j, epsilon = 1e-12);
                }
            }
        }
    }

    #[test]
    fn test_compute_jacobian_3d() {
        let g = example_geometry_3d();
        let points = triangle_points();
        let evaluator = g.get_evaluator(points.data());

        let mut computed_jacobian = rlst_dynamic_array2!(f64, [3, 2]);
        for (cell_i, jacobian) in [
            vec![vec![1.0, 1.0], vec![0.0, 1.0], vec![1.0, 0.0]],
            vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![0.0, 0.0]],
        ]
        .iter()
        .enumerate()
        {
            for point_i in 0..points.shape()[0] {
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
        let g = example_geometry_3d();
        let points = triangle_points();
        let evaluator = g.get_evaluator(points.data());

        let mut computed_normal = vec![0.0; 3];
        for (cell_i, normal) in [
            vec![
                -1.0 / f64::sqrt(3.0),
                1.0 / f64::sqrt(3.0),
                1.0 / f64::sqrt(3.0),
            ],
            vec![0.0, 0.0, 1.0],
        ]
        .iter()
        .enumerate()
        {
            for point_i in 0..points.shape()[0] {
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
    fn test_midpoint_flat() {
        let g = example_geometry_flat();

        let mut midpoint = vec![0.0; 3];
        for (cell_i, point) in [
            vec![2.0 / 3.0, 1.0 / 3.0, 0.0],
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
    fn test_midpoint_3d() {
        let g = example_geometry_3d();

        let mut midpoint = vec![0.0; 3];
        for (cell_i, point) in [
            vec![2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
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
}
