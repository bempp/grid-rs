//! Definition of a grid

use crate::traits::cell::CellType;
use crate::traits::point::PointType;
use crate::types::cell_iterator::CellIterator;
use crate::types::point_iterator::PointIterator;
use crate::types::{CellLocalIndexPair, Float};

use super::ReferenceMapType;

pub trait GridType: std::marker::Sized {
    type T: Float;
    type IndexType: std::fmt::Debug + Eq + Copy;

    type Point<'a>: PointType
    where
        Self: 'a;

    type Cell<'a>: CellType
    where
        Self: 'a;

    type ReferenceMap<'a>: ReferenceMapType
    where
        Self: 'a;

    fn number_of_vertices(&self) -> usize;

    fn number_of_points(&self) -> usize;

    fn number_of_cells(&self) -> usize;

    fn point_index_from_id(&self, id: usize) -> usize;
    fn point_id_from_index(&self, index: usize) -> usize;

    fn cell_index_from_id(&self, id: usize) -> usize;
    fn cell_id_from_index(&self, index: usize) -> usize;

    fn point_from_index(&self, index: usize) -> Self::Point<'_>;

    fn cell_from_index(&self, index: usize) -> Self::Cell<'_>;

    fn iter_points<Iter: std::iter::Iterator<Item = usize>>(
        &self,
        index_iter: Iter,
    ) -> PointIterator<'_, Self, Iter> {
        PointIterator::new(index_iter, self)
    }

    fn iter_all_points(&self) -> PointIterator<'_, Self, std::ops::Range<usize>> {
        self.iter_points(0..self.number_of_points())
    }

    fn iter_cells<Iter: std::iter::Iterator<Item = usize>>(
        &self,
        index_iter: Iter,
    ) -> CellIterator<'_, Self, Iter> {
        CellIterator::new(index_iter, self)
    }

    fn iter_all_cells(&self) -> CellIterator<'_, Self, std::ops::Range<usize>> {
        self.iter_cells(0..self.number_of_cells())
    }

    fn reference_to_physical_map<'a>(
        &'a self,
        reference_points: &'a [Self::T],
    ) -> Self::ReferenceMap<'a>;

    fn vertex_to_cells(
        &self,
        vertex_index: Self::IndexType,
    ) -> &[CellLocalIndexPair<Self::IndexType>];

    fn edge_to_cells(&self, edge_index: Self::IndexType) -> &[CellLocalIndexPair<Self::IndexType>];

    fn face_to_cells(&self, face_index: Self::IndexType) -> &[CellLocalIndexPair<Self::IndexType>];
}
