//! Definition of a cell.

use crate::traits::topology::TopologyType;

use super::GridType;

pub trait CellType {
    type Grid: GridType;

    type Topology<'a>: TopologyType
    where
        Self: 'a;

    fn id(&self) -> usize;
    fn index(&self) -> usize;

    fn topology(&self) -> Self::Topology<'_>;

    fn grid(&self) -> &Self::Grid;
}
