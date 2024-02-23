//! Definition of a grid

use num::traits::Float;

pub trait Grid {
    type T: Float;

    type Vertex<'a>
    where
        Self: 'a;
    type Edge<'a>
    where
        Self: 'a;
    type Face<'a>
    where
        Self: 'a;

    type VertexIterator<'a>: std::iter::Iterator<Item = Self::Vertex>
    where
        Self: 'a;

    type EdgeIterator<'a>: std::iter::Iterator<Item = Self::Edge>
    where
        Self: 'a;

    type FaceIterator<'a>: std::iter::Iterator<Item = Self::Face>
    where
        Self: 'a;

    type CellIterator<'a>: std::iter::Iterator<Item = Self::CellType>
    where
        Self: 'a;

    fn iter_cell(&self) -> Self::CellIterator<'_>;

    fn iter_vertex(&self) -> Self::VertexIterator<'_>;

    fn iter_face(&self) -> Self::FaceIterator<'_>;

    fn iter_edge(&self) -> Self::EdgeIterator<'_>;
}
