//! Trait definitions.

use num::traits::Float;

pub trait Grid {
    type T: Float;

    type Vertex<'a>
    where
        Self: 'a;
    type Edge;
    type Face;
    type Cell<'a>
    where
        Self: 'a;

    type VertexIterator<'v>: std::iter::Iterator<Item = Self::Vertex<'v>>
    where
        Self: 'v;

    type CellIterator<'v>: std::iter::Iterator<Item = Self::Cell<'v>>
    where
        Self: 'v;

    fn iter_vertices(&self) -> Self::VertexIterator<'_>;
    fn iter_cells(&self) -> Self::CellIterator<'_>;
}
