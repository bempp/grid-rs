// //! Definition of a reference cell.

// use num::traits::Float;
// use std::fmt;

// use super::map::Map;

// #[derive(Debug)]
// pub struct InvalidConnectivity;

// impl fmt::Display for InvalidConnectivity {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         write!(f, "Invalid connectivity")
//     }
// }

// impl std::error::Error for InvalidConnectivity {}

// #[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
// #[repr(u8)]
// pub enum CellType {
//     Point,
//     Interval,
//     Triangle,
//     Quadrilateral,
//     Tetrahedron,
//     Hexahedron,
//     Prism,
//     Pyramid,
// }

// /// A 0- to 3- dimensional reference cell
// pub trait Cell {
//     type T: Float;

//     /// The dimension of the reference cell (eg a triangle's dimension is 2, tetrahedron's dimension is 3)
//     fn dim(&self) -> usize;

//     /// Check if the cell is a simplex
//     fn is_simplex(&self) -> bool;

//     /// The vertices of the cell
//     ///
//     fn vertices(&self) -> std::slice::Iter<'_, Self::Vertex>;

//     /// The edges of the cell
//     ///
//     /// The first 2 components are the vertex numbers of the endpoints of the first edge, the next 2 the second edge, and so on.
//     fn edges(&self) -> std::slice::Iter<'_, Self::Edge>;

//     /// The faces of the cell
//     ///
//     /// The first `self.faces_nvertices()[0]` components are the vertex numbers of vertices of the first face, the next `self.faces_nvertices()[1]` the second face, and so on.
//     fn faces(&self) -> std::slice::Iter<'_, Self::Face>;

//     /// The number of vertices adjacent to each face
//     fn faces_nvertices(&self) -> &[usize];

//     /// The number of entities of dimension `dim`
//     fn entity_count(&self, dim: usize) -> usize {
//         match dim {
//             0 => self.vertex_count(),
//             1 => self.edge_count(),
//             2 => self.face_count(),
//             3 => self.volume_count(),
//             _ => 0,
//         }
//     }

//     /// The cell type of the entity with dimension `dim` and cell index `index`.
//     fn entity_types(&self, cell_index: usize, dim: usize) -> CellType;

//     /// The number of vertices
//     fn vertex_count(&self) -> usize;

//     /// The number of edges
//     fn edge_count(&self) -> usize;

//     /// The number of faces
//     fn face_count(&self) -> usize;

//     /// The number of volumes
//     fn volume_count(&self) -> usize;

//     fn connected_vertices(
//         &self,
//         entity_dim: usize,
//         entity_number: usize,
//     ) -> Result<std::slice::Iter<'_, Self::Vertex>, InvalidConnectivity>;

//     fn connected_edges(
//         &self,
//         entity_dim: usize,
//         entity_number: usize,
//     ) -> Result<std::slice::Iter<'_, Self::Edge>, InvalidConnectivity>;

//     fn connected_faces(
//         &self,
//         entity_dim: usize,
//         entity_number: usize,
//     ) -> Result<std::slice::Iter<'_, Self::Face>, InvalidConnectivity>;

//     /// The reference cell type
//     fn cell_type(&self) -> CellType;
// }
