//! Trait definitions.

use num::traits::Float;

pub trait VertexType {
    type T: Float;

    fn coords(&self, data: &mut [Self::T]);
    fn index(&self) -> usize;
    fn id(&self) -> usize;
}

pub trait CellType {
    type Topology<'a>: TopologyType
    where
        Self: 'a;

    fn id(&self) -> usize;
    fn index(&self) -> usize;

    fn topology(&self) -> Self::Topology<'_>;
}

pub struct VertexIterator<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> {
    iter: Iter,
    grid: &'a Grid,
}

impl<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> std::iter::Iterator
    for VertexIterator<'a, Grid, Iter>
{
    type Item = Grid::Vertex<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.iter.next() {
            Some(self.grid.vertex_from_index(index))
        } else {
            None
        }
    }
}

pub struct CellIterator<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> {
    iter: Iter,
    grid: &'a Grid,
}

impl<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> std::iter::Iterator
    for CellIterator<'a, Grid, Iter>
{
    type Item = Grid::Cell<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.iter.next() {
            Some(self.grid.cell_from_index(index))
        } else {
            None
        }
    }
}

pub trait TopologyType {
    type Grid: GridType;
    type VertexIndexIter<'a>: std::iter::Iterator<Item = usize>
    where
        Self: 'a;

    fn vertex_indices(&self) -> Self::VertexIndexIter<'_>;
}

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
        indices: Iter,
    ) -> VertexIterator<'_, Self, Iter> {
        VertexIterator {
            iter: indices,
            grid: self,
        }
    }

    fn iter_all_vertices(&self) -> VertexIterator<'_, Self, std::ops::Range<usize>> {
        self.iter_vertices(0..self.number_of_vertices())
    }

    fn iter_cells<Iter: std::iter::Iterator<Item = usize>>(
        &self,
        indices: Iter,
    ) -> CellIterator<'_, Self, Iter> {
        CellIterator {
            iter: indices,
            grid: self,
        }
    }

    fn iter_all_cells(&self) -> CellIterator<'_, Self, std::ops::Range<usize>> {
        self.iter_cells(0..self.number_of_cells())
    }
}
