//! Definition of a cell.

use crate::traits::topology::TopologyType;

pub trait CellType {
    type Topology<'a>: TopologyType
    where
        Self: 'a;

    fn id(&self) -> usize;
    fn index(&self) -> usize;

    fn topology(&self) -> Self::Topology<'_>;
}
