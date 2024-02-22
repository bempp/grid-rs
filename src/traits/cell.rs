//! Definition of a reference cell.

use num::traits::Float;
use std::fmt;

#[derive(Debug)]
pub struct InvalidConnectivity;

impl fmt::Display for InvalidConnectivity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Invalid connectivity")
    }
}

impl std::error::Error for InvalidConnectivity {}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
#[repr(u8)]
pub enum ReferenceCellType {
    Point = 0,
    Interval = 1,
    Triangle = 2,
    Quadrilateral = 3,
    Tetrahedron = 4,
    Hexahedron = 5,
    Prism = 6,
    Pyramid = 7,
}

// <start>
//To discuss with Timo:
// In Basix, all cell information is obtained via free functions, eg:
fn dim(cell: ReferenceCellType) {
    match(cell) {
        ReferenceCellType::Point => 0,
        ReferenceCellType::Interval => 1,
        ReferenceCellType::Triangle => 2,
        ReferenceCellType::Quadrilateral => 2,
        ReferenceCellType::Tetrahedron => 3,
        ReferenceCellType::Hexahedron => 3,
        ReferenceCellType::Prism => 3,
        ReferenceCellType::Pyramid => 3,
    }
}
fn is_simplex(cell: ReferenceCellType) {
    match(cell) {
        ReferenceCellType::Point => true,
        ReferenceCellType::Interval => true,
        ReferenceCellType::Triangle => true,
        ReferenceCellType::Quadrilateral => false,
        ReferenceCellType::Tetrahedron => true,
        ReferenceCellType::Hexahedron => false,
        ReferenceCellType::Prism => false,
        ReferenceCellType::Pyramid => false,
    }
}
// This is ok, as these functions should only be called rarely, and means we don't have to
// convert a ReferenceCellType to a struct whenever we need one of these functions.
// Another question: should we implement these in the same file as we define the enum?
// <end>

/// A 0- to 3- dimensional reference cell
pub trait ReferenceCell {
    type T: Float;

    /// The dimension of the reference cell (eg a triangle's dimension is 2, tetrahedron's dimension is 3)
    fn dim(&self) -> usize;

    /// Check if the cell is a simplex
    fn is_simplex(&self) -> bool;

    /// The vertices of the cell
    ///
    /// The first dim components represent the first vertex, the next dim the second vertex, and so on.
    fn vertices<'a>(&'a self) -> std::slice::Iter<'a, usize>;

    /// Ths midpoint of the cell
    fn midpoint(&self, midpoint: &mut [Self::T]);

    /// The edges of the cell
    ///
    /// The first 2 components are the vertex numbers of the endpoints of the first edge, the next 2 the second edge, and so on.
    fn edges(&self) -> std::slice::Iter<'_, (usize, usize)>;

    /// The faces of the cell
    ///
    /// The first `self.faces_nvertices()[0]` components are the vertex numbers of vertices of the first face, the next `self.faces_nvertices()[1]` the second face, and so on.
    fn faces(&self) -> std::slice::Iter<'_, &[usize]>;

    /// The number of vertices adjacent to each face
    fn faces_nvertices(&self) -> &[usize];

    /// The number of entities of dimension `dim`
    fn entity_count(&self, dim: usize) -> usize {
        match dim {
            0 => self.vertex_count(),
            1 => self.edge_count(),
            2 => self.face_count(),
            3 => self.volume_count(),
            _ => 0,
        }
    }

    /// The cell type of the entity with dimension `dim` and cell index `index`.
    fn entity_types(&self, cell_index: usize, dim: usize) -> ReferenceCellType;

    /// The number of vertices
    fn vertex_count(&self) -> usize;

    /// The number of edges
    fn edge_count(&self) -> usize;

    /// The number of faces
    fn face_count(&self) -> usize;

    /// The number of volumes
    fn volume_count(&self) -> usize;

    /// Get the entities connected to an entity
    ///
    /// This function returns a list of entity numbers of entities of dimension `connected_dim` that are attached to the entity numbered `entity_dim` of number entity_number.
    /// For example connectivity(1, 0, 2) will return a list of faces (2D entities) that are connected to edge (1D entity) 0.
    fn connected_vertices(
        &self,
        entity_dim: usize,
        entity_number: usize,
    ) -> Result<std::slice::Iter<'_, usize>, InvalidConnectivity>;

    fn connected_edges(
        &self,
        entity_dim: usize,
        entity_number: usize,
    ) -> Result<std::slice::Iter<'_, (usize, usize)>, InvalidConnectivity>;

    fn connected_faces(
        &self,
        entity_dim: usize,
        entity_number: usize,
    ) -> Result<std::slice::Iter<'_, &[usize]>, InvalidConnectivity>;

    /// The reference cell type
    fn cell_type(&self) -> ReferenceCellType;

    /// The reference cell label
    fn label(&self) -> &'static str;
}
