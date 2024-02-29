//! Iterator over vertices

use crate::traits::grid::GridType;

pub struct PointIterator<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> {
    iter: Iter,
    grid: &'a Grid,
}

impl<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> PointIterator<'a, Grid, Iter> {
    pub fn new(iter: Iter, grid: &'a Grid) -> Self {
        Self { iter, grid }
    }
}

impl<'a, Grid: GridType, Iter: std::iter::Iterator<Item = usize>> std::iter::Iterator
    for PointIterator<'a, Grid, Iter>
{
    type Item = Grid::Point<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.iter.next() {
            Some(self.grid.point_from_index(index))
        } else {
            None
        }
    }
}
