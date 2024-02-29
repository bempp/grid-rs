//! A surface triangle grid.

use std::collections::HashMap;

use crate::types::Float;
use rlst_common::types::Scalar;

use crate::traits::*;

use super::{cell::TriangleCell, vertex::TriangleVertex};

pub struct TriangleSurfaceGrid<T: Float + Scalar> {
    pub(crate) vertices: Vec<[T; 3]>,
    pub(crate) cells: Vec<[usize; 3]>,
    ids_from_vertex_indices: Vec<usize>,
    ids_from_cell_indices: Vec<usize>,
    vertex_ids: HashMap<usize, usize>,
    cell_ids: HashMap<usize, usize>,
    edge_connectivity: HashMap<(usize, usize), (usize, Vec<usize>)>,
    pub(crate) cell_to_edges: Vec<[usize; 3]>,
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
            edge_connectivity: HashMap::new(),
            cell_to_edges: Vec::new(),
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
            edge_connectivity: HashMap::new(),
            cell_to_edges: Vec::new(),
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

    pub fn finalize(&mut self) {
        self.create_edges();
    }

    fn create_edges(&mut self) {
        let mut nedges: usize = 0;
        self.cell_to_edges
            .resize_with(self.cells.len(), Default::default);
        let edge_local: [(usize, usize); 3] = [(1, 2), (0, 2), (0, 1)];
        for (cell_index, cell_vertices) in self.cells.iter().enumerate() {
            for (local_index, &(first, second)) in edge_local.iter().enumerate() {
                let mut first = cell_vertices[first];
                let mut second = cell_vertices[second];
                if first > second {
                    std::mem::swap(&mut first, &mut second);
                }
                if let Some((edge_index, adjacent_cells)) =
                    self.edge_connectivity.get_mut(&(first, second))
                {
                    adjacent_cells.push(cell_index);
                    self.cell_to_edges[cell_index][local_index] = *edge_index;
                } else {
                    let mut adjacent_cells = Vec::<usize>::with_capacity(2);
                    adjacent_cells.push(cell_index);
                    self.cell_to_edges[cell_index][local_index] = nedges;
                    self.edge_connectivity
                        .insert((first, second), (nedges, adjacent_cells));
                    nedges += 1;
                }
            }
        }
    }
}

impl<T: Float + Scalar> GridType for TriangleSurfaceGrid<T> {
    // type T = T;

    type Vertex<'a> = TriangleVertex<'a, T> where Self: 'a;
    type Cell<'a> = TriangleCell<'a, T> where Self: 'a;

    type Edge = [usize; 2];

    type Face = ();

    fn number_of_vertices(&self) -> usize {
        self.vertices.len()
    }
    fn number_of_cells(&self) -> usize {
        self.cells.len()
    }

    fn vertex_index_from_id(&self, id: usize) -> usize {
        self.vertex_ids[&id]
    }
    fn vertex_id_from_index(&self, index: usize) -> usize {
        self.ids_from_vertex_indices[index]
    }

    fn cell_index_from_id(&self, id: usize) -> usize {
        self.cell_ids[&id]
    }
    fn cell_id_from_index(&self, index: usize) -> usize {
        self.ids_from_cell_indices[index]
    }

    fn vertex_from_index(&self, index: usize) -> Self::Vertex<'_> {
        TriangleVertex::new(
            self.vertex_id_from_index(index),
            index,
            &self.vertices[index],
        )
    }

    fn cell_from_index(&self, index: usize) -> Self::Cell<'_> {
        TriangleCell::new(index, self)
    }
}

#[cfg(test)]
mod test {
    use crate::traits::*;

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

        grid.finalize();

        let mut coords = [0.0; 3];
        for vertex in grid.iter_all_vertices() {
            vertex.coords(coords.as_mut_slice());
            println!("{:#?}", coords);
        }

        for cell in grid.iter_all_cells() {
            for (local_index, (vertex_index, edge_index)) in cell
                .topology()
                .vertex_indices()
                .zip(cell.topology().edge_indices())
                .enumerate()
            {
                println!(
                    "Cell: {}, Vertex: {}, {}, Edge: {}, {}, Volume: {}",
                    cell.index(),
                    local_index,
                    vertex_index,
                    local_index,
                    edge_index,
                    cell.geometry().volume(),
                )
            }
        }
    }
}
