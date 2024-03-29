//! Definition of a cell.

use crate::traits::topology::TopologyType;

use super::{GeometryType, GridType};

pub trait CellType {
    //! A cell

    /// The type of the grid that the cell is part of
    type Grid: GridType;
    /// The type of the cell topology
    type Topology<'a>: TopologyType
    where
        Self: 'a;
    /// The type of the cell geometry
    type Geometry<'a>: GeometryType
    where
        Self: 'a;

    /// The id of the cell
    fn id(&self) -> usize;

    /// The index of the cell
    fn index(&self) -> usize;

    /// Get the cell's topology
    fn topology(&self) -> Self::Topology<'_>;

    /// Get the grid that the cell is part of
    fn grid(&self) -> &Self::Grid;

    /// Get the cell's geometry
    fn geometry(&self) -> Self::Geometry<'_>;
}
