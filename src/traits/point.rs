//! Definition of a vertex

use crate::types::Float;

pub trait PointType {
    type T: Float;

    fn coords(&self, data: &mut [Self::T]);
    fn index(&self) -> usize;
    fn id(&self) -> usize;
}
