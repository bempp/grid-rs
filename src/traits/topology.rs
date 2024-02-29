//! Definition of a topology

use crate::traits::grid::GridType;

pub trait TopologyType {
    type Grid: GridType;
    type IndexType: std::fmt::Debug + Eq + Copy;
    type VertexIndexIter<'a>: std::iter::Iterator<Item = Self::IndexType>
    where
        Self: 'a;

    type EdgeIndexIter<'a>: std::iter::Iterator<Item = Self::IndexType>
    where
        Self: 'a;

    fn vertex_indices(&self) -> Self::VertexIndexIter<'_>;

    fn edge_indices(&self) -> Self::EdgeIndexIter<'_>;
}
