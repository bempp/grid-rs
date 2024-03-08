//! Grid builder
use crate::traits::grid::GridType;
use num::Float;

pub trait Builder<const GDIM: usize> {
    //! Object that can be used to build a mesh

    /// The geometric/physical dimension
    const GDIM: usize = GDIM;
    /// The type of the grid that the builder created
    type GridType: GridType;
    /// The floating point type used for coordinates
    type T: Float;
    /// The type of the data that is input to add a cell
    type CellData;
    /// The type of the data that must be provided when initialising the builder
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
    fn create_grid(self) -> Self::GridType;
}
