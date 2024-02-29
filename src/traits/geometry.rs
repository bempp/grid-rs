//! Geometry of a physical cell

use crate::types::Float;

use super::{GridType, VertexType};

pub trait GeometryType {
    type Grid: GridType;

    type VertexIterator<'iter>: std::iter::Iterator<Item = <Self::Grid as GridType>::Vertex<'iter>>
    where
        Self: 'iter;

    type PointsIterator<'iter>: std::iter::Iterator<Item = <Self::Grid as GridType>::Vertex<'iter>>
    where
        Self: 'iter;

    fn physical_dimension(&self) -> usize;

    fn midpoint(&self, point: &mut [<Self::Grid as GridType>::T]);

    fn diameter(&self) -> <Self::Grid as GridType>::T;

    fn volume(&self) -> <Self::Grid as GridType>::T;

    fn vertices(&self) -> Self::VertexIterator<'_>;

    fn points(&self) -> Self::PointsIterator<'_>;
}
