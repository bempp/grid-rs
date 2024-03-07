//! Grid builder
use crate::traits::grid::GridType;
use num::Float;

pub trait Builder<const GDIM: usize> {
    const GDIM: usize = GDIM;
    type GridType: GridType;
    type T: Float;
    type CellData;
    type GridMetadata;

    /// Create a new grid builder
    fn new(data: Self::GridMetadata) -> Self;

    /// Create a new grid builder with capacaty for a given number of points and cells
    fn new_with_capacity(npoints: usize, ncells: usize, data: Self::GridMetadata) -> Self;

    /// Add a point to the grid
    fn add_point(&mut self, id: usize, data: [Self::T; GDIM]);

    /// Add a cell to the grid
    fn add_cell(&mut self, id: usize, cell_data: Self::CellData);

    /// Create the grid
    fn create_grid(&self) -> Self::GridType;
}
