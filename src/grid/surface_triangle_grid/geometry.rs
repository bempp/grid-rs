//! Definition of the Geometry

use std::{iter::Copied, marker::PhantomData};

use num::{Float, Zero};
use rlst::rlst_static_array;
use rlst_common::types::Scalar;
use rlst_dense::{rlst_array_from_slice1, rlst_array_from_slice_mut1};

use crate::{
    traits::{GeometryType, GridType},
    types::vertex_iterator::PointIterator,
};

use super::grid::TriangleSurfaceGrid;

pub struct TriangleGeometry<'a, T: Float + Scalar> {
    cell_index: usize,
    diameter: T,
    volume: T,
    midpoint: [T; 3],
    grid: &'a TriangleSurfaceGrid<T>,
}

impl<'a, T: Float + Scalar> TriangleGeometry<'a, T> {
    pub fn new(grid: &'a TriangleSurfaceGrid<T>, cell_index: usize) -> Self {
        let vertex_indices = &grid.cells[cell_index];
        let v0 = rlst_array_from_slice1!(T, grid.vertices[vertex_indices[0]].as_slice(), [3]);
        let v1 = rlst_array_from_slice1!(T, grid.vertices[vertex_indices[1]].as_slice(), [3]);
        let v2 = rlst_array_from_slice1!(T, grid.vertices[vertex_indices[2]].as_slice(), [3]);

        let mut a = rlst_static_array!(T, 3);
        let mut b = rlst_static_array!(T, 3);
        let mut c = rlst_static_array!(T, 3);
        let mut cross = rlst_static_array!(T, 3);

        let mut midpoint = [
            <T as Zero>::zero(),
            <T as Zero>::zero(),
            <T as Zero>::zero(),
        ];

        let mut midpoint_arr = rlst_array_from_slice_mut1!(T, midpoint.as_mut(), [3]);

        a.fill_from(v1.view() - v0.view());
        b.fill_from(v2.view() - v0.view());
        c.fill_from(v2.view() - v1.view());

        let a_norm = a.view().norm_2();
        let b_norm = b.view().norm_2();
        let c_norm = c.view().norm_2();

        midpoint_arr.fill_from(
            (v0.view() + v1.view() + v2.view()).scalar_mul(T::from_f64(1.0 / 3.0).unwrap()),
        );

        a.cross(b.view(), cross.view_mut());

        let volume = num::cast::<f64, T::Real>(0.5).unwrap() * cross.view().norm_2();

        let s = num::cast::<f64, T::Real>(0.5).unwrap() * (a_norm + b_norm + c_norm);

        let diameter = num::cast::<f64, T::Real>(2.0).unwrap()
            * num::Float::sqrt(((s - a_norm) * (s - b_norm) * (s - c_norm)) / s);

        Self {
            cell_index,
            diameter: num::cast(diameter).unwrap(),
            volume: num::cast(volume).unwrap(),
            midpoint,
            grid,
        }
    }
}

impl<'a, T: Float + Scalar> GeometryType for TriangleGeometry<'a, T> {
    type Grid = TriangleSurfaceGrid<T>;

    type VertexIterator<'iter> =
        PointIterator<'iter, Self::Grid, Copied<std::slice::Iter<'iter, usize>>> where Self: 'iter;

    type PointsIterator<'iter> = Self::VertexIterator<'iter> where Self: 'iter;

    fn physical_dimension(&self) -> usize {
        3
    }

    fn midpoint(&self, point: &mut [<Self::Grid as GridType>::T]) {
        for (ind_out, &ind_in) in point.iter_mut().zip(self.midpoint.iter()) {
            *ind_out = ind_in
        }
    }

    fn diameter(&self) -> <Self::Grid as GridType>::T {
        self.diameter
    }

    fn volume(&self) -> <Self::Grid as GridType>::T {
        self.volume
    }

    fn vertices(&self) -> Self::VertexIterator<'_> {
        self.grid
            .iter_points(self.grid.cells[self.cell_index].iter().copied())
    }

    fn points(&self) -> Self::PointsIterator<'_> {
        self.vertices()
    }
}
