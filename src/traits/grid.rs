//! Definition of a grid.

use crate::traits::cell::ReferenceCellType;
use num::traits::Float;

/// The ownership of a mesh entity
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub enum Ownership {
    Owned,
    Ghost(usize, usize),
}

/// The toplogy of a grid.
///
/// This provides information about which mesh entities are connected to other mesh entities
pub trait Topology {
    /// The dimension of the topology (eg a triangle's dimension is 2, tetrahedron's dimension is 3)
    fn dim(&self) -> usize;

    /// Return the index map from the input cell numbers to the storage numbers
    fn index_map(&self) -> &[usize];

    /// The number of entities of dimension `dim`
    fn entity_count(&self, dim: usize) -> usize;

    /// The indices of the vertices of the cell with topological index `index`
    fn cell(&self, index: usize) -> Option<&[usize]>;

    /// The cell type of the cell with topological index `index`
    fn cell_type(&self, index: usize) -> Option<ReferenceCellType>;

    /// All entity types of the given dimension that are included in the grid
    fn entity_types(&self, dim: usize) -> &[ReferenceCellType];

    /// Get the connectivity of entities of type `type0` and entities of type `type1`
    fn connectivity(&self, type0: ReferenceCellType, type1: ReferenceCellType) -> &[usize];

    /// Get the ownership of a mesh entity
    fn entity_ownership(&self, dim: usize, index: usize) -> Ownership;
}
