//! Definition of a topology

use crate::traits::grid::GridType;

pub trait TopologyType {
    type Grid: GridType;
    type VertexIndexIter<'a>: std::iter::Iterator<Item = usize>
    where
        Self: 'a;

    fn vertex_indices(&self) -> Self::VertexIndexIter<'_>;
}
