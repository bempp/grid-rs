//! Geometry of a physical cell

use crate::types::Float;

pub trait GeometryType {
    type T: Float;

    type PointIterator<'a>: std::iter::Iterator<Item = &'a [Self::T]>
    where
        Self: 'a;

    fn physical_dimension(&self) -> usize;

    fn midpoint(&self, point: &mut [Self::T]);

    fn diameter(&self) -> Self::T;

    fn volume(&self) -> Self::T;

    fn corners(&self) -> Self::PointIterator<'_>;
}
