//! Type definitions

pub mod cell_iterator;
pub mod point_iterator;

pub use num::Float;

pub use bempp_element::cell::ReferenceCellType;

pub struct CellLocalIndexPair {
    pub cell: usize,
    pub local_index: usize,
}

impl CellLocalIndexPair {
    pub fn new(cell: usize, local_index: usize) -> Self {
        Self { cell, local_index }
    }
}
