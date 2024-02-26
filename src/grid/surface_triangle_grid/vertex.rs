//! Vertex of a triangle grid

use crate::types::Float;
use rlst_common::types::Scalar;

use crate::traits::*;

pub struct TriangleVertex<'a, T: Float + Scalar> {
    id: usize,
    index: usize,
    data: &'a [T; 3],
}

impl<'a, T: Float + Scalar> TriangleVertex<'a, T> {
    pub fn new(id: usize, index: usize, data: &'a [T; 3]) -> Self {
        Self { id, index, data }
    }
}

impl<'a, T: Float + Scalar> VertexType for TriangleVertex<'a, T> {
    type T = T;

    fn coords(&self, data: &mut [Self::T]) {
        for (out_data, &in_data) in data.iter_mut().zip(self.data) {
            *out_data = in_data;
        }
    }

    fn index(&self) -> usize {
        self.index
    }

    fn id(&self) -> usize {
        self.id
    }
}
