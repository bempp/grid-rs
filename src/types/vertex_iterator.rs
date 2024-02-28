//! Iterator over vertices

use crate::traits::grid::GridType;

pub struct VertexIterator<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> {
    iter: Iter,
    grid: &'a Grid,
}

impl<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> VertexIterator<'a, Grid, Iter> {
    pub fn new(iter: Iter, grid: &'a Grid) -> Self {
        Self { iter, grid }
    }
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
