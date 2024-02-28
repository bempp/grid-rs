//! Definition of a vertex

use crate::types::Float;

pub trait VertexType {
    type T: Float;

    fn coords(&self, data: &mut [Self::T]);
    fn index(&self) -> usize;
    fn id(&self) -> usize;
}
