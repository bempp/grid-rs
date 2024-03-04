//! Implementation of grid geometry

use crate::grid_impl::traits::Geometry;
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
    traits::{RandomAccessByRef, Shape, UnsafeRandomAccessByRef},
};

/// Geometry of a serial grid
pub struct SerialSingleElementGeometry<T: Float + Scalar> {
    dim: usize,
    index_map: Vec<usize>,
    pub(crate) coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
    pub(crate) cells: Vec<usize>,
    pub(crate) element: CiarletElement<T>,
    midpoints: Vec<Vec<T>>,
    diameters: Vec<T>,
    volumes: Vec<T>,
}

unsafe impl<T: Float + Scalar> Sync for SerialSingleElementGeometry<T> {}

impl<T: Float + Scalar> SerialSingleElementGeometry<T> {
    pub fn new(
        coordinates: Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>,
        cells_input: &[usize],
        element: CiarletElement<T>,
    ) -> Self {
        let dim = coordinates.shape()[1];
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

impl<T: Float + Scalar> Geometry for SerialSingleElementGeometry<T> {
    type IndexType = usize;
    type T = T;
    type Element = CiarletElement<T>;

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

pub struct GeometryEvaluatorSingleElement<'a, T: Float + Scalar> {
    geometry: &'a SerialSingleElementGeometry<T>,
    tdim: usize,
    table: Array<T, BaseArray<T, VectorContainer<T>, 4>, 4>,
}

impl<'a, T: Float + Scalar> GeometryEvaluatorSingleElement<'a, T> {
    fn new<Points: RandomAccessByRef<2, Item = T> + Shape<2>>(
        geometry: &'a SerialSingleElementGeometry<T>,
        points: &Points,
    ) -> Self {
        let tdim = reference_cell::dim(
            // TODO: remove this match once bempp-rs and grid-rs use the same ReferenceCellType
            match geometry.element.cell_type() {
                bempp_element::cell::ReferenceCellType::Interval => ReferenceCellType::Interval,
                bempp_element::cell::ReferenceCellType::Triangle => ReferenceCellType::Triangle,
                bempp_element::cell::ReferenceCellType::Quadrilateral => {
                    ReferenceCellType::Quadrilateral
                }
                _ => {
                    panic!("Unsupported cell type: {:?}", geometry.element.cell_type());
                }
            },
        );
        assert_eq!(tdim, points.shape()[1]);
        let mut table = rlst_dynamic_array4!(
            T,
            geometry.element.tabulate_array_shape(1, points.shape()[0])
        );
        geometry.element.tabulate(points, 1, &mut table);
        Self {
            geometry,
            tdim,
            table,
        }
    }

    pub fn compute_point(&self, cell_index: usize, point_index: usize, point: &mut [T]) {
        let gdim = self.geometry.dim();
        let element_npts = self.geometry.element.dim();
        assert_eq!(point.len(), gdim);

        for component in point.iter_mut() {
            *component = T::from(0.0).unwrap();
        }
        for (i, v) in self.geometry.cells[cell_index * element_npts..]
            .iter()
            .take(element_npts)
            .enumerate()
        {
            let t = unsafe { *self.table.get_unchecked([0, point_index, i, 0]) };
            for (j, component) in point.iter_mut().enumerate() {
                *component += unsafe { *self.geometry.coordinates.get_unchecked([*v, j]) } * t;
            }
        }
    }

    pub fn compute_jacobian(&self, cell_index: usize, point_index: usize, jacobian: &mut [T]) {
        let gdim = self.geometry.dim();
        let tdim = self.tdim;
        let element_npts = self.geometry.element.dim();
        assert_eq!(jacobian.len(), gdim * tdim);

        for component in jacobian.iter_mut() {
            *component = T::from(0.0).unwrap();
        }
        for (i, v) in self.geometry.cells[cell_index * element_npts..]
            .iter()
            .take(element_npts)
            .enumerate()
        {
            for gd in 0..gdim {
                for td in 0..tdim {
                    jacobian[td * gdim + gd] +=
                        unsafe { *self.geometry.coordinates.get_unchecked([*v, gd]) }
                            * unsafe { *self.table.get_unchecked([1 + td, point_index, i, 0]) };
                }
            }
        }
    }

    pub fn compute_normal(&self, cell_index: usize, point_index: usize, normal: &mut [T]) {
        let gdim = self.geometry.dim();
        let tdim = self.tdim;
        assert_eq!(tdim, 2); // TODO: remove this
        assert_eq!(tdim, gdim - 1);
        let mut jacobian = [T::from(0.0).unwrap(); 6];

        self.compute_jacobian(cell_index, point_index, &mut jacobian[..]);

        for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
            normal[i] = jacobian[j] * jacobian[3 + k] - jacobian[k] * jacobian[3 + j];
        }
        let size = Scalar::sqrt(normal.iter().map(|&x| Scalar::powi(x, 2)).sum::<T>());
        for i in normal.iter_mut() {
            *i /= size;
        }
    }
}

impl<T: Float + Scalar> SerialSingleElementGeometry<T> {
    pub fn get_evaluator<Points: RandomAccessByRef<2, Item = T> + Shape<2>>(
        &self,
        points: &Points,
    ) -> GeometryEvaluatorSingleElement<T> {
        GeometryEvaluatorSingleElement::<T>::new(self, points)
    }
}

#[cfg(test)]
mod test {
    use crate::grid_impl::single_element_grid::geometry::*;
    use approx::*;
    use bempp_element::element::{create_element, ElementFamily};
    use bempp_traits::element::Continuity;
    use rlst_dense::{
        rlst_dynamic_array2,
        traits::{RandomAccessMut, RawAccessMut},
    };

    fn example_geometry_2d() -> SerialSingleElementGeometry<f64> {
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
        SerialSingleElementGeometry::new(points, &[0, 1, 2, 0, 2, 3], p1triangle)
    }

    fn example_geometry_3d() -> SerialSingleElementGeometry<f64> {
        let p2triangle = create_element(
            ElementFamily::Lagrange,
            bempp_element::cell::ReferenceCellType::Triangle,
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
        SerialSingleElementGeometry::new(points, &[0, 2, 8, 5, 4, 1, 0, 8, 6, 7, 3, 4], p2triangle)
    }

    #[test]
    fn test_counts() {
        let g = example_geometry_2d();
        assert_eq!(g.point_count(), 4);
        assert_eq!(g.cell_count(), 2);
    }

    #[test]
    fn test_cell_points() {
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

    fn triangle_points() -> Array<f64, BaseArray<f64, VectorContainer<f64>, 2>, 2> {
        let mut points = rlst_dynamic_array2!(f64, [2, 2]);
        *points.get_mut([0, 0]).unwrap() = 0.2;
        *points.get_mut([0, 1]).unwrap() = 0.5;
        *points.get_mut([1, 0]).unwrap() = 0.6;
        *points.get_mut([1, 1]).unwrap() = 0.1;
        points
    }

    #[test]
    fn test_compute_point_2d() {
        let g = example_geometry_2d();
        let points = triangle_points();

        let evaluator = g.get_evaluator(&points);
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
        let g = example_geometry_3d();
        let points = triangle_points();
        let evaluator = g.get_evaluator(&points);

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
        let g = example_geometry_3d();
        let points = triangle_points();
        let evaluator = g.get_evaluator(&points);

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
                evaluator.compute_jacobian(cell_i, point_i, &mut computed_jacobian.data_mut());
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
        let evaluator = g.get_evaluator(&points);

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
}
