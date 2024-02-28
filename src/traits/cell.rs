//! Definition of a cell.

use crate::traits::topology::TopologyType;

use super::{GeometryType, GridType};

pub trait CellType {
    type Grid: GridType;

    type Topology<'a>: TopologyType
    where
        Self: 'a;

    type Geometry<'a>: GeometryType
    where
        Self: 'a;

    fn id(&self) -> usize;
    fn index(&self) -> usize;

    fn topology(&self) -> Self::Topology<'_>;

    fn grid(&self) -> &Self::Grid;

    fn geometry(&self) -> Self::Geometry<'_>;
}
