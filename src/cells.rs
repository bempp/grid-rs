//! Cells

use num::traits::Float;
use crate::traits::cell::{ReferenceCellType, ReferenceCell, InvalidConnectivity};


struct Interval<T: Float> {
    vertices: Vec<T>
}

impl<T: Float> Interval<T> {
    fn new() -> Self {
        let zero = T::zero();
        let one = T::one();
        Interval {
            vertices: vec![zero, one],
        }
    }
}

impl<T: Float> ReferenceCell for Interval<T> {
    type T = T;

    fn dim(&self) -> usize {
        1
    }
    fn is_simplex(&self) -> bool {
        true
    }
    fn vertices<'a>(&'a self) -> std::slice::Iter<'a, &[Self::T]> {
        // TODO
        [].iter()
        //[self.vertices[0..1], self.vertices[1..2]].iter()
    }
    fn midpoint(&self, midpoint: &mut [Self::T]) {
        midpoint[0] = (self.vertices[0] + self.vertices[1]) / T::from(2.0).unwrap();
    }
    fn edges(&self) -> std::slice::Iter<'_, (usize, usize)> {
        [(0,1)].iter()
    }
    fn faces(&self) -> std::slice::Iter<'_, &[usize]> {
        [].iter()
    }
    fn faces_nvertices(&self) -> &[usize] {
        &[]
    }
    fn entity_types(&self, cell_index: usize, dim: usize) -> ReferenceCellType {
        // TODO
        ReferenceCellType::Interval
    }
    fn vertex_count(&self) -> usize {
        2
    }
    fn edge_count(&self) -> usize {
        1
    }
    fn face_count(&self) -> usize {
        0
    }
    fn volume_count(&self) -> usize {
        0
    }
    fn connected_vertices(
        &self,
        entity_dim: usize,
        entity_number: usize,
    ) -> Result<std::slice::Iter<'_, usize>, InvalidConnectivity> {
        // TODO
        Err(InvalidConnectivity)
    }
    fn connected_edges(
        &self,
        entity_dim: usize,
        entity_number: usize,
    ) -> Result<std::slice::Iter<'_, (usize, usize)>, InvalidConnectivity> {
        // TODO
        Err(InvalidConnectivity)
    }
    fn connected_faces(
        &self,
        entity_dim: usize,
        entity_number: usize,
    ) -> Result<std::slice::Iter<'_, &[usize]>, InvalidConnectivity> {
        // TODO
        Err(InvalidConnectivity)
    }
    fn cell_type(&self) -> ReferenceCellType {
        // TODO
        ReferenceCellType::Interval
    }
    fn label(&self) -> &'static str {
        "interval"
    }
}

#[cfg(test)]
mod test {
    use crate::cells::*;

    #[test]
    fn test_interval() {
        let cell = Interval::<f64>::new();
        assert_eq!(cell.vertex_count(), 2);
        assert_eq!(cell.edge_count(), 1);
        assert_eq!(cell.face_count(), 0);
        assert_eq!(cell.volume_count(), 0);
    }
}
