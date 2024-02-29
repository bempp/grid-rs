//! Type definitions

pub mod cell_iterator;
pub mod vertex_iterator;

pub use num::Float;

pub enum ReferenceCellType {
    Interval,
    Triangle,
    Quadrilateral,
    Tetrahedron,
    Hexahedron,
    Prism,
    Pyramid,
}

pub struct CellLocalIndexPair {
    pub cell: usize,
    pub local_index: usize,
}

impl CellLocalIndexPair {
    pub fn new(cell: usize, local_index: usize) -> Self {
        Self { cell, local_index }
    }
}
