//! Surface triangle grid

pub mod cell;
pub mod geometry;
pub mod grid;
pub mod reference_map;
pub mod topology;
pub mod vertex;

// use std::{collections::HashMap, iter::Copied, marker::PhantomData, ops::Add};

// use num::{traits::Float, One};
// use rlst::{rlst_dynamic_array2, rlst_static_array, Array};
// use rlst_common::types::Scalar;
// use rlst_dense::{
//     array::DynamicArray,
//     rlst_array_from_slice1, rlst_array_from_slice_mut1, rlst_array_from_slice_mut2,
//     traits::{DefaultIteratorMut, RandomAccessMut, RawAccess},
// };
// use rlst_proc_macro::rlst_static_type;

// use crate::traits::{CellType, GridType, TopologyType, VertexType};

// pub struct TriangleVertex<'a, T: Float + Scalar> {
//     id: usize,
//     index: usize,
//     data: &'a [T; 3],
// }

// impl<'a, T: Float + Scalar> VertexType for TriangleVertex<'a, T> {
//     type T = T;

//     fn coords(&self, data: &mut [Self::T]) {
//         for (out_data, &in_data) in data.iter_mut().zip(self.data) {
//             *out_data = in_data;
//         }
//     }

//     fn index(&self) -> usize {
//         self.index
//     }

//     fn id(&self) -> usize {
//         self.id
//     }
// }

// pub struct TriangleCell<'a, T: Float + Scalar> {
//     id: usize,
//     index: usize,
//     data: &'a [usize; 3],
//     _marker: PhantomData<T>,
// }

// impl<'a, T: Float + Scalar> CellType for TriangleCell<'a, T> {
//     type Topology<'top> = TriangleTopology<'top, T> where Self: 'top;

//     fn id(&self) -> usize {
//         self.id
//     }

//     fn index(&self) -> usize {
//         self.index
//     }

//     fn topology(&self) -> Self::Topology<'_> {
//         TriangleTopology::new(self)
//     }
// }

// pub struct TriangleTopology<'a, T: Float + Scalar> {
//     cell: &'a TriangleCell<'a, T>,
//     _marker: PhantomData<T>,
// }

// impl<'a, T: Float + Scalar> TriangleTopology<'a, T> {
//     pub fn new(cell: &'a TriangleCell<'a, T>) -> Self {
//         Self {
//             cell,
//             _marker: PhantomData,
//         }
//     }
// }

// impl<'a, T: Float + Scalar> TopologyType for TriangleTopology<'a, T> {
//     type Grid = TriangleSurfaceGrid<T>;

//     type VertexIndexIter<'iter> = Copied<std::slice::Iter<'iter, usize>> where Self: 'iter;

//     fn vertex_indices(&self) -> Self::VertexIndexIter<'_> {
//         self.cell.data.iter().copied()
//     }
// }

// // pub struct TriangleMap<T: Float + Scalar> {
// //     v0: rlst_static_type!(T, 3),
// //     gradient: rlst_static_type!(T, 3, 2),
// // }

// // impl<T: Float + Scalar> TriangleMap<T> {
// //     pub fn new(cell: &TriangleCell, grid: &TriangleSurfaceGrid<T>) -> Self {
// //         let mut v0 = rlst_static_array!(T, 3);
// //         v0.fill_from(grid.vertices[cell.vertices()[0]].data().view());

// //         let mut gradient = rlst_static_array!(T, 3, 2);
// //         let vertices = [
// //             &grid.vertices[cell.vertices()[1]],
// //             &grid.vertices[cell.vertices()[2]],
// //         ];

// //         for (mut gradient_col, vertex) in gradient.col_iter_mut().zip(vertices) {
// //             gradient_col.fill_from(
// //                 vertex.data().view().add(
// //                     v0.view()
// //                         .add(vertex.data().view().scalar_mul(-<T as One>::one())),
// //                 ),
// //             );
// //         }

// //         Self { v0, gradient }
// //     }
// // }

// pub struct TriangleSurfaceGrid<T: Float + Scalar> {
//     vertices: Vec<[T; 3]>,
//     cells: Vec<[usize; 3]>,
//     ids_from_vertex_indices: Vec<usize>,
//     ids_from_cell_indices: Vec<usize>,
//     vertex_ids: HashMap<usize, usize>,
//     cell_ids: HashMap<usize, usize>,
// }

// impl<T: Float + Scalar> Default for TriangleSurfaceGrid<T> {
//     fn default() -> Self {
//         Self::new()
//     }
// }

// impl<T: Float + Scalar> TriangleSurfaceGrid<T> {
//     pub fn new() -> Self {
//         Self {
//             vertices: Vec::new(),
//             cells: Vec::new(),
//             ids_from_vertex_indices: Vec::new(),
//             ids_from_cell_indices: Vec::new(),
//             vertex_ids: HashMap::new(),
//             cell_ids: HashMap::new(),
//         }
//     }

//     pub fn new_with_capacity(nvertices: usize, ncells: usize) -> Self {
//         Self {
//             vertices: Vec::with_capacity(nvertices),
//             cells: Vec::with_capacity(ncells),
//             ids_from_vertex_indices: Vec::with_capacity(nvertices),
//             ids_from_cell_indices: Vec::with_capacity(ncells),
//             vertex_ids: HashMap::new(),
//             cell_ids: HashMap::new(),
//         }
//     }

//     pub fn add_vertex(&mut self, id: usize, data: [T; 3]) {
//         let index = self.vertices.len();
//         self.vertices.push(data);
//         self.ids_from_vertex_indices.push(id);
//         self.vertex_ids.insert(id, index);
//     }

//     pub fn add_cell(&mut self, id: usize, vertex_ids: [usize; 3]) {
//         let index = self.cells.len();
//         self.cells.push([
//             self.vertex_ids[&vertex_ids[0]],
//             self.vertex_ids[&vertex_ids[1]],
//             self.vertex_ids[&vertex_ids[2]],
//         ]);
//         self.ids_from_cell_indices.push(id);
//         self.cell_ids.insert(id, index);
//     }
// }

// impl<T: Float + Scalar> GridType for TriangleSurfaceGrid<T> {
//     type T = T;

//     type Vertex<'a> = TriangleVertex<'a, T> where Self: 'a;
//     type Cell<'a> = TriangleCell<'a, T> where Self: 'a;

//     type Edge = [usize; 2];

//     type Face = ();

//     fn number_of_vertices(&self) -> usize {
//         self.vertices.len()
//     }
//     fn number_of_cells(&self) -> usize {
//         self.cells.len()
//     }

//     fn vertex_index_from_id(&self, id: usize) -> usize {
//         self.vertex_ids[&id]
//     }
//     fn vertex_id_from_index(&self, index: usize) -> usize {
//         self.ids_from_vertex_indices[index]
//     }

//     fn cell_index_from_id(&self, id: usize) -> usize {
//         self.cell_ids[&id]
//     }
//     fn cell_id_from_index(&self, index: usize) -> usize {
//         self.ids_from_cell_indices[index]
//     }

//     fn vertex_from_index(&self, index: usize) -> Self::Vertex<'_> {
//         TriangleVertex {
//             id: self.vertex_id_from_index(index),
//             index,
//             data: &self.vertices[index],
//         }
//     }

//     fn cell_from_index(&self, index: usize) -> Self::Cell<'_> {
//         TriangleCell {
//             id: self.cell_id_from_index(index),
//             index,
//             data: &self.cells[index],
//             _marker: PhantomData,
//         }
//     }
// }

// #[cfg(test)]
// mod test {
//     use crate::traits::{CellType, GridType, TopologyType};

//     use super::TriangleSurfaceGrid;

//     #[test]
//     fn test_grid() {
//         let mut grid = TriangleSurfaceGrid::<f64>::new();

//         grid.add_vertex(0, [0.0, 0.0, 0.0]);
//         grid.add_vertex(1, [1.0, 0.0, 0.0]);
//         grid.add_vertex(2, [0.0, 1.0, 0.0]);
//         grid.add_vertex(3, [1.0, 1.0, 0.0]);

//         grid.add_cell(0, [0, 1, 2]);
//         grid.add_cell(0, [1, 3, 2]);

//         for vertex in grid.iter_all_vertices() {
//             println!("{:#?}", vertex.data)
//         }

//         for cell in grid.iter_all_cells() {
//             for (local_index, vertex_index) in cell.topology().vertex_indices().enumerate() {
//                 println!(
//                     "Cell: {}, Vertex: {}, {}",
//                     cell.index(),
//                     local_index,
//                     vertex_index
//                 )
//             }
//         }
//     }
// }
