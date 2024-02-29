//! Definition of a grid

use crate::traits::cell::CellType;
use crate::traits::point::PointType;
use crate::types::cell_iterator::CellIterator;
use crate::types::vertex_iterator::PointIterator;
use num::Float;

pub trait GridType: std::marker::Sized {
    type T: Float;

    type Point<'a>: PointType
    where
        Self: 'a;
    type Edge;
    type Face;
    type Cell<'a>: CellType
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
        self.iter_points(0..self.number_of_vertices())
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
}
