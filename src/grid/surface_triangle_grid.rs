//! Surface triangle grid

use std::{collections::HashMap, marker::PhantomData, ops::Add};

use num::{traits::Float, One};
use rlst::{rlst_dynamic_array2, rlst_static_array, Array};
use rlst_common::types::Scalar;
use rlst_dense::{
    array::DynamicArray,
    rlst_array_from_slice1, rlst_array_from_slice_mut1, rlst_array_from_slice_mut2,
    traits::{DefaultIteratorMut, RandomAccessMut, RawAccess},
};
use rlst_proc_macro::rlst_static_type;

use crate::traits::Grid;

pub struct TriangleVertex<'a, T: Float + Scalar> {
    id: usize,
    index: usize,
    data: &'a [T; 3],
}

pub struct TriangleCell<'a> {
    id: usize,
    index: usize,
    data: &'a [usize; 3],
}

pub struct VertexIterator<'a, T: Float + Scalar> {
    grid: &'a TriangleSurfaceGrid<T>,
    current: usize,
    finished: bool,
}

pub struct CellIterator<'a, T: Float + Scalar> {
    grid: &'a TriangleSurfaceGrid<T>,
    current: usize,
    finished: bool,
}

impl<'a, T: Float + Scalar> VertexIterator<'a, T> {
    pub fn new(grid: &'a TriangleSurfaceGrid<T>) -> Self {
        Self {
            grid,
            current: 0,
            finished: false,
        }
    }
}

impl<'a, T: Float + Scalar> CellIterator<'a, T> {
    pub fn new(grid: &'a TriangleSurfaceGrid<T>) -> Self {
        Self {
            grid,
            current: 0,
            finished: false,
        }
    }
}

impl<'a, T: Float + Scalar> std::iter::Iterator for VertexIterator<'a, T> {
    type Item = TriangleVertex<'a, T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished || self.grid.vertices.is_empty() {
            None
        } else {
            let out = Some(TriangleVertex {
                id: self.grid.ids_from_vertex_indices[self.current],
                index: self.current,
                data: &self.grid.vertices[self.current],
            });

            self.current += 1;
            self.finished = self.current == self.grid.vertices.len();

            out
        }
    }
}

impl<'a, T: Float + Scalar> std::iter::Iterator for CellIterator<'a, T> {
    type Item = TriangleCell<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished || self.grid.vertices.is_empty() {
            None
        } else {
            let out = Some(TriangleCell {
                id: self.grid.ids_from_cell_indices[self.current],
                index: self.current,
                data: &self.grid.cells[self.current],
            });

            self.current += 1;
            self.finished = self.current == self.grid.cells.len();

            out
        }
    }
}

// pub struct TriangleMap<T: Float + Scalar> {
//     v0: rlst_static_type!(T, 3),
//     gradient: rlst_static_type!(T, 3, 2),
// }

// impl<T: Float + Scalar> TriangleMap<T> {
//     pub fn new(cell: &TriangleCell, grid: &TriangleSurfaceGrid<T>) -> Self {
//         let mut v0 = rlst_static_array!(T, 3);
//         v0.fill_from(grid.vertices[cell.vertices()[0]].data().view());

//         let mut gradient = rlst_static_array!(T, 3, 2);
//         let vertices = [
//             &grid.vertices[cell.vertices()[1]],
//             &grid.vertices[cell.vertices()[2]],
//         ];

//         for (mut gradient_col, vertex) in gradient.col_iter_mut().zip(vertices) {
//             gradient_col.fill_from(
//                 vertex.data().view().add(
//                     v0.view()
//                         .add(vertex.data().view().scalar_mul(-<T as One>::one())),
//                 ),
//             );
//         }

//         Self { v0, gradient }
//     }
// }

pub struct TriangleSurfaceGrid<T: Float + Scalar> {
    vertices: Vec<[T; 3]>,
    cells: Vec<[usize; 3]>,
    ids_from_vertex_indices: Vec<usize>,
    ids_from_cell_indices: Vec<usize>,
    vertex_ids: HashMap<usize, usize>,
    cell_ids: HashMap<usize, usize>,
}

impl<T: Float + Scalar> Default for TriangleSurfaceGrid<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Float + Scalar> TriangleSurfaceGrid<T> {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            cells: Vec::new(),
            ids_from_vertex_indices: Vec::new(),
            ids_from_cell_indices: Vec::new(),
            vertex_ids: HashMap::new(),
            cell_ids: HashMap::new(),
        }
    }

    pub fn new_with_capacity(nvertices: usize, ncells: usize) -> Self {
        Self {
            vertices: Vec::with_capacity(nvertices),
            cells: Vec::with_capacity(ncells),
            ids_from_vertex_indices: Vec::with_capacity(nvertices),
            ids_from_cell_indices: Vec::with_capacity(ncells),
            vertex_ids: HashMap::new(),
            cell_ids: HashMap::new(),
        }
    }

    pub fn add_vertex(&mut self, id: usize, data: [T; 3]) {
        let index = self.vertices.len();
        self.vertices.push(data);
        self.ids_from_vertex_indices.push(id);
        self.vertex_ids.insert(id, index);
    }

    pub fn add_cell(&mut self, id: usize, vertex_ids: [usize; 3]) {
        let index = self.cells.len();
        self.cells.push([
            self.vertex_ids[&vertex_ids[0]],
            self.vertex_ids[&vertex_ids[1]],
            self.vertex_ids[&vertex_ids[2]],
        ]);
        self.ids_from_cell_indices.push(id);
        self.cell_ids.insert(id, index);
    }
}

impl<T: Float + Scalar> Grid for TriangleSurfaceGrid<T> {
    type T = T;

    type Vertex<'a> = TriangleVertex<'a, T> where Self: 'a;
    type Cell<'a> = TriangleCell<'a> where Self: 'a;

    type Edge = [usize; 2];

    type Face = ();

    type VertexIterator<'v> = VertexIterator<'v, T>
    where
        Self: 'v;

    fn iter_vertices(&self) -> Self::VertexIterator<'_> {
        VertexIterator::new(self)
    }

    type CellIterator<'v> = CellIterator<'v, T>
    where
        Self: 'v;

    fn iter_cells(&self) -> Self::CellIterator<'_> {
        CellIterator::new(self)
    }
}

#[cfg(test)]
mod test {
    use crate::traits::Grid;

    use super::TriangleSurfaceGrid;

    #[test]
    fn test_grid() {
        let mut grid = TriangleSurfaceGrid::<f64>::new();

        grid.add_vertex(0, [0.0, 0.0, 0.0]);
        grid.add_vertex(1, [1.0, 0.0, 0.0]);
        grid.add_vertex(2, [0.0, 1.0, 0.0]);
        grid.add_vertex(3, [1.0, 1.0, 0.0]);

        grid.add_cell(0, [0, 1, 2]);
        grid.add_cell(0, [1, 3, 2]);

        for vertex in grid.iter_cells() {
            println!("{:#?}", vertex.data)
        }
    }
}
