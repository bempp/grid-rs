//! Traits used in the implementation of a grid

use crate::traits::cell::ReferenceCellType;
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
    /// The dimension of the topology (eg a triangle's dimension is 2, tetrahedron's dimension is 3)
    fn dim(&self) -> usize;

    /// Return the index map from the input cell numbers to the storage numbers
    fn index_map(&self) -> &[usize];

    /// The number of entities of dimension `dim`
    fn entity_count(&self, etype: ReferenceCellType) -> usize;

    /// The indices of the vertices of the cell with topological index `index`
    fn cell(&self, index: usize) -> Option<&[usize]>;

    /// The cell type of the cell with topological index `index`
    fn cell_type(&self, index: usize) -> Option<ReferenceCellType>;

    /// All entity types of the given dimension that are included in the grid
    fn entity_types(&self, dim: usize) -> &[ReferenceCellType];

    // NO: If type0 = Point and type1 = Triangle, there will not be the same number of indices for each `type0`
    /// Get the connectivity of entities of type `type0` and entities of type `type1`
    // fn connectivity(&self, type0: ReferenceCellType, type1: ReferenceCellType) -> &[usize];

    // Alternative not requiring an adjacency list. Possible issue: unable to get connectivity information for all non-cell entities at once
    // This may be ok, as in assembly we're only using the connectivity from cells to vertices, which can be obtained via cell_entities.
    /// Get the indices of entities of type `etype` that are connected to each cell of type `cell_type`
    fn cell_entities(
        &self,
        cell_type: ReferenceCellType,
        etype: ReferenceCellType,
    ) -> Option<&[usize]>;
    /// Get the indices of entities of dimension `dim` that are connected to the entity of type `etype` with index `index`
    fn connectivity(&self, etype: ReferenceCellType, index: usize, dim: usize) -> Option<&[usize]>;

    /// Get the ownership of a mesh entity
    fn entity_ownership(&self, dim: usize, index: usize) -> Ownership;
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
    fn elements_and_cells(&self) -> std::slice::Iter<'_, (&Self::Element, &[usize])>;

    // In progress: functions to compute geometry information
    // Need these to be able to compute something for every cell in turn without needing to re-tabulate the geometry element every time
    // Question: what's better: passing in 1D slices or rlst matrices? For now I've used 1D, but current bempp-rs uses rlst
    // I'm not a big fan of the user having to sort out these tables

    /// Compute the coordinates of a set of points on a physical cell
    ///
    /// `table` should be the geometry element tabulated at the points on the reference cell, ie:
    ///    table = rlst_dynamic_array4!(T, self.cell_element(cell).tabulate_array_shape(0, npts));
    ///    self.cell_element(cell).tabulate(points, 0, table.data_mut());
    /// For this function, `table` may include 0 derivatives.
    // We want to run this function on multiple cells with the same points, so don't want to retabuilate inside the function every time?
    fn compute_points(&self, table: &[Self::T], cell: usize, physical_points: &mut [Self::T]);

    /// Compute the normals and jacobian determinants at a set of points on a physical cell
    ///
    /// `table` should be the geometry element tabulated at the points on the reference cell, ie:
    ///    table = rlst_dynamic_array4!(T, self.cell_element(cell).tabulate_array_shape(0, npts));
    ///    self.cell_element(cell).tabulate(points, 0, table.data_mut());
    /// For this function, `table` must include at least 1 derivative.
    fn compute_normals_and_jacobian_determinants(
        &self,
        table: &[Self::T],
        cell: usize,
        physical_points: &mut [Self::T],
    );
}

/// A grid
pub trait Grid {
    /// The type that implements [Topology]
    type Topology: Topology;

    /// The type that implements [Geometry]
    type Geometry: Geometry;

    /// Get the grid topology (See [Topology])
    fn topology(&self) -> &Self::Topology;

    /// Get the grid geometry (See [Geometry])
    fn geometry(&self) -> &Self::Geometry;

    // Check if the grid is stored in serial
    fn is_serial(&self) -> bool;
}
