//! Definition of a grid

use crate::traits::cell::CellType;
use crate::traits::vertex::VertexType;
use crate::types::cell_iterator::CellIterator;
use crate::types::vertex_iterator::VertexIterator;
use crate::types::Float;

pub trait GridType: std::marker::Sized {
    type T: Float;

    type Vertex<'a>: VertexType
    where
        Self: 'a;
    type Edge;
    type Face;
    type Cell<'a>: CellType
    where
        Self: 'a;

    fn number_of_vertices(&self) -> usize;
    fn number_of_cells(&self) -> usize;

    fn vertex_index_from_id(&self, id: usize) -> usize;
    fn vertex_id_from_index(&self, index: usize) -> usize;

    fn cell_index_from_id(&self, id: usize) -> usize;
    fn cell_id_from_index(&self, index: usize) -> usize;

    fn vertex_from_index(&self, index: usize) -> Self::Vertex<'_>;

    fn cell_from_index(&self, index: usize) -> Self::Cell<'_>;

    fn iter_vertices<Iter: std::iter::Iterator<Item = usize>>(
        &self,
        index_iter: Iter,
    ) -> VertexIterator<'_, Self, Iter> {
        VertexIterator::new(index_iter, self)
    }

    fn iter_all_vertices(&self) -> VertexIterator<'_, Self, std::ops::Range<usize>> {
        self.iter_vertices(0..self.number_of_vertices())
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
