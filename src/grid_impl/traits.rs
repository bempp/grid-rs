//! Traits used in the implementation of a grid

use crate::reference_cell::ReferenceCellType;
use bempp_traits::element::FiniteElement;
use num::Float;

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub enum Ownership {
    Owned,
    Ghost(usize, usize),
}

/// The toplogy of a grid.
///
/// This provides information about which mesh entities are connected to other mesh entities
pub trait Topology {
    type IndexType: std::fmt::Debug + Eq + Copy;

    /// The dimension of the topology (eg a triangle's dimension is 2, tetrahedron's dimension is 3)
    fn dim(&self) -> usize;

    /// Return the index map from the input cell numbers to the storage numbers
    fn index_map(&self) -> &[Self::IndexType];

    /// The number of entities of type `etype`
    fn entity_count(&self, etype: ReferenceCellType) -> usize;

    /// The number of entities of dimension `dim`
    fn entity_count_by_dim(&self, dim: usize) -> usize;

    /// The indices of the vertices of the cell with topological index `index`
    fn cell(&self, index: Self::IndexType) -> Option<&[Self::IndexType]>;

    /// The cell type of the cell with topological index `index`
    fn cell_type(&self, index: Self::IndexType) -> Option<ReferenceCellType>;

    /// All entity types of the given dimension that are included in the grid
    fn entity_types(&self, dim: usize) -> &[ReferenceCellType];

    /// Get the indices of entities and types of entities that are connected to each cell of type `cell_type`
    fn cell_entities(&self, cell_type: ReferenceCellType, dim: usize)
        -> Option<&[Self::IndexType]>;

    /// Get the indices and types of entities of dimension `dim` that are connected to the entity of type `etype` with index `index`
    fn connectivity(&self, index: Self::IndexType, dim: usize) -> Option<&[Self::IndexType]>;

    /// Get the ownership of a mesh entity
    fn entity_ownership(&self, dim: usize, index: Self::IndexType) -> Ownership;

    /// Extract a flat index from an IndexType
    fn extract_index(&self, index: Self::IndexType) -> usize;
}

/// The geometry of a grid
///
/// This provides information about the physical locations of mesh points in space
pub trait Geometry {
    type T: Float;
    type Element: FiniteElement;

    /// The geometric dimension
    fn dim(&self) -> usize;

    /// Return the index map from the input cell numbers to the storage numbers
    fn index_map(&self) -> &[usize];

    /// Get one of the coordinates of a point
    fn coordinate(&self, point_index: usize, coord_index: usize) -> Option<&Self::T>;

    /// The number of points stored in the geometry
    fn point_count(&self) -> usize;

    /// Get the vertex numbers of a cell
    fn cell_vertices(&self, index: usize) -> Option<&[usize]>;

    /// The number of cells
    fn cell_count(&self) -> usize;

    /// Get the element used to represent a cell
    fn cell_element(&self, index: usize) -> Option<&Self::Element>;

    /// Get the number of distinct geometry elements
    fn element_count(&self) -> usize;
    /// Get the `i`th element
    fn element(&self, i: usize) -> Option<&Self::Element>;
    /// Get the cells associated with the `i`th element
    fn cells(&self, i: usize) -> Option<&[usize]>;

    // ... or would it be better to replace the above 3 functions with an Iter?
    // fn elements_and_cells(&self) -> std::slice::Iter<'_, (&Self::Element, &[usize])>;

    /// Midpoint of a cell
    fn midpoint(&self, index: usize, point: &mut [Self::T]);

    /// Diameter of a cell
    fn diameter(&self, index: usize) -> Self::T;

    /// Volume of a cell
    fn volume(&self, index: usize) -> Self::T;
}

/// A grid
pub trait Grid {
    type T: Float;

    /// The type that implements [Topology]
    type Topology: Topology;

    /// The type that implements [Geometry]
    type Geometry: Geometry<T = Self::T>;

    /// Get the grid topology (See [Topology])
    fn topology(&self) -> &Self::Topology;

    /// Get the grid geometry (See [Geometry])
    fn geometry(&self) -> &Self::Geometry;

    // Check if the grid is stored in serial
    fn is_serial(&self) -> bool;
}
