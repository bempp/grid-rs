//! Implementation of grid geometry

use num::Float;
use crate::grid::traits::Geometry;
//use crate::grid::traits::{Ownership, Geometry};
//use crate::reference_cell;
//use crate::traits::cell::ReferenceCellType;
use bempp_element::element::CiarletElement;

//use std::collections::HashMap;

/// Geometry of a serial grid
pub struct SerialGeometry<T: Float> {
    pub id: T,
}

unsafe impl<T: Float> Sync for SerialGeometry<T> {}

impl<T: Float> Geometry for SerialGeometry<T> {

    type T = T;
    type Element = CiarletElement;

    fn dim(&self) -> usize { panic!(); }

    fn index_map(&self) -> &[usize] { panic!(); }

    fn coordinate(&self, _point_index: usize, _coord_index: usize) -> Option<&Self::T> { panic!(); }

    fn point_count(&self) -> usize { panic!(); }

    fn cell_vertices(&self, _index: usize) -> Option<&[usize]> { panic!(); }

    fn cell_count(&self) -> usize { panic!(); }

    fn cell_element(&self, _index: usize) -> Option<&Self::Element> { panic!(); }

    fn element_count(&self) -> usize { panic!(); }
    fn element(&self, _i: usize) -> Option<&Self::Element> { panic!(); }
    fn cells(&self, _i: usize) -> Option<&[usize]> { panic!(); }

    fn elements_and_cells(&self) -> std::slice::Iter<'_, (&Self::Element, &[usize])> { panic!(); }

    fn compute_points(&self, _table: &[Self::T], _cell: usize, _physical_points: &mut [Self::T]) { panic!(); }

    fn compute_normals_and_jacobian_determinants(
        &self,
        _table: &[Self::T],
        _cell: usize,
        _physical_points: &mut [Self::T],
    ) { panic!(); }

}
