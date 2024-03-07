//! Type definitions

pub mod cell_iterator;
pub mod point_iterator;

pub use num::Float;

pub use bempp_element::cell::ReferenceCellType;

#[derive(Debug, Clone)]
pub struct CellLocalIndexPair<IndexType: std::fmt::Debug + Eq + Copy> {
    pub cell: IndexType,
    pub local_index: usize,
}

impl<IndexType: std::fmt::Debug + Eq + Copy> CellLocalIndexPair<IndexType> {
    pub fn new(cell: IndexType, local_index: usize) -> Self {
        Self { cell, local_index }
    }
}
