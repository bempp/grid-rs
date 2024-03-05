//! A surface triangle grid.

use std::collections::HashMap;

use crate::types::{CellLocalIndexPair, Float};
use rlst::rlst_static_array;
use rlst_common::types::Scalar;
use rlst_dense::rlst_array_from_slice1;
use rlst_proc_macro::rlst_static_type;

use crate::traits::*;

use super::{
    cell::TriangleCell,
    reference_map::TriangleReferenceMap,
    vertex::TriangleVertex,
};

pub struct TriangleSurfaceGrid<T: Float + Scalar> {
    pub(crate) vertices: Vec<[T; 3]>,
    pub(crate) cells: Vec<[usize; 3]>,
    ids_from_vertex_indices: Vec<usize>,
    ids_from_cell_indices: Vec<usize>,
    vertex_ids: HashMap<usize, usize>,
    cell_ids: HashMap<usize, usize>,
    edge_to_cells: Vec<Vec<CellLocalIndexPair>>,
    point_to_cells: Vec<Vec<CellLocalIndexPair>>,
    pub(crate) cell_to_edges: Vec<[usize; 3]>,
    pub(crate) jacobians: Vec<rlst_static_type!(T, 3, 2)>,
    pub(crate) volumes: Vec<T>,
    pub(crate) diameters: Vec<T>,
    pub(crate) normals: Vec<rlst_static_type!(T, 3)>,
    pub(crate) midpoints: Vec<rlst_static_type!(T, 3)>,
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
            edge_to_cells: Vec::new(),
            point_to_cells: Vec::new(),
            cell_to_edges: Vec::new(),
            jacobians: Vec::new(),
            volumes: Vec::new(),
            diameters: Vec::new(),
            normals: Vec::new(),
            midpoints: Vec::new(),
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
            edge_to_cells: Vec::new(),
            point_to_cells: Vec::new(),
            cell_to_edges: Vec::new(),
            jacobians: Vec::with_capacity(ncells),
            volumes: Vec::with_capacity(ncells),
            diameters: Vec::with_capacity(ncells),
            normals: Vec::with_capacity(ncells),
            midpoints: Vec::with_capacity(ncells),
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
        self.create_edge_connectivity();
        self.create_point_connectivity();
        self.compute_geometry_information();
    }

    fn create_edge_connectivity(&mut self) {
        let mut nedges: usize = 0;
        self.cell_to_edges
            .resize_with(self.cells.len(), Default::default);
        let edge_local: [(usize, usize); 3] = [(1, 2), (0, 2), (0, 1)];
        let mut edge_connectivity =
            HashMap::<(usize, usize), (usize, Vec<CellLocalIndexPair>)>::new();
        for (cell_index, cell_vertices) in self.cells.iter().enumerate() {
            for (local_index, &(first, second)) in edge_local.iter().enumerate() {
                let mut first = cell_vertices[first];
                let mut second = cell_vertices[second];
                if first > second {
                    std::mem::swap(&mut first, &mut second);
                }
                if let Some((edge_index, adjacent_cells)) =
                    edge_connectivity.get_mut(&(first, second))
                {
                    adjacent_cells.push(CellLocalIndexPair::new(cell_index, local_index));
                    self.cell_to_edges[cell_index][local_index] = *edge_index;
                } else {
                    let mut adjacent_cells = Vec::<CellLocalIndexPair>::with_capacity(2);
                    adjacent_cells.push(CellLocalIndexPair::new(cell_index, local_index));
                    self.cell_to_edges[cell_index][local_index] = nedges;
                    edge_connectivity.insert((first, second), (nedges, adjacent_cells));
                    nedges += 1;
                }
            }
        }
        self.edge_to_cells.resize_with(nedges, Default::default);
        for (edge_index, cell_connectivity_vec) in edge_connectivity.into_values() {
            self.edge_to_cells[edge_index] = cell_connectivity_vec;
        }
    }

    fn create_point_connectivity(&mut self) {
        self.point_to_cells
            .resize_with(self.vertices.len(), Default::default);

        for (cell_index, cell) in self.cells.iter().enumerate() {
            for (local_index, &vertex_index) in cell.iter().enumerate() {
                self.point_to_cells[vertex_index]
                    .push(CellLocalIndexPair::new(cell_index, local_index));
            }
        }
    }

    fn compute_geometry_information(&mut self) {
        for vertex_indices in &self.cells {
            let v0 = rlst_array_from_slice1!(T, self.vertices[vertex_indices[0]].as_slice(), [3]);
            let v1 = rlst_array_from_slice1!(T, self.vertices[vertex_indices[1]].as_slice(), [3]);
            let v2 = rlst_array_from_slice1!(T, self.vertices[vertex_indices[2]].as_slice(), [3]);

            let mut a = rlst_static_array!(T, 3);
            let mut b = rlst_static_array!(T, 3);
            let mut c = rlst_static_array!(T, 3);
            let mut cross = rlst_static_array!(T, 3);

            let mut midpoint = rlst_static_array!(T, 3);

            let mut jacobian = rlst_static_array!(T, 3, 2);

            a.fill_from(v1.view() - v0.view());
            b.fill_from(v2.view() - v0.view());
            c.fill_from(v2.view() - v1.view());

            jacobian.view_mut().slice(1, 0).fill_from(a.view());
            jacobian.view_mut().slice(1, 1).fill_from(b.view());

            let a_norm = a.view().norm_2();
            let b_norm = b.view().norm_2();
            let c_norm = c.view().norm_2();

            midpoint.fill_from(
                (v0.view() + v1.view() + v2.view()).scalar_mul(T::from_f64(1.0 / 3.0).unwrap()),
            );

            a.cross(b.view(), cross.view_mut());
            let normal_length = cross.view().norm_2();
            cross.scale_in_place(T::one() / num::cast(normal_length).unwrap());

            let volume = num::cast::<f64, T::Real>(0.5).unwrap() * normal_length;

            let s = num::cast::<f64, T::Real>(0.5).unwrap() * (a_norm + b_norm + c_norm);

            let diameter = num::cast::<f64, T::Real>(2.0).unwrap()
                * num::Float::sqrt(((s - a_norm) * (s - b_norm) * (s - c_norm)) / s);

            self.jacobians.push(jacobian);
            self.volumes.push(num::cast(volume).unwrap());
            self.diameters.push(num::cast(diameter).unwrap());
            self.normals.push(cross);
            self.midpoints.push(midpoint);
        }
    }
}

impl<T: Float + Scalar> GridType for TriangleSurfaceGrid<T> {
    type T = T;

    type Point<'a> = TriangleVertex<'a, T> where Self: 'a;
    type Cell<'a> = TriangleCell<'a, T> where Self: 'a;

    type ReferenceMap<'a> = TriangleReferenceMap<'a, T>
    where
        Self: 'a;

    fn number_of_vertices(&self) -> usize {
        self.vertices.len()
    }

    fn number_of_points(&self) -> usize {
        self.vertices.len()
    }

    fn number_of_cells(&self) -> usize {
        self.cells.len()
    }

    fn point_index_from_id(&self, id: usize) -> usize {
        self.vertex_ids[&id]
    }
    fn point_id_from_index(&self, index: usize) -> usize {
        self.ids_from_vertex_indices[index]
    }

    fn cell_index_from_id(&self, id: usize) -> usize {
        self.cell_ids[&id]
    }
    fn cell_id_from_index(&self, index: usize) -> usize {
        self.ids_from_cell_indices[index]
    }

    fn point_from_index(&self, index: usize) -> Self::Point<'_> {
        TriangleVertex::new(
            self.point_id_from_index(index),
            index,
            &self.vertices[index],
        )
    }

    fn cell_from_index(&self, index: usize) -> Self::Cell<'_> {
        TriangleCell::new(index, self)
    }

    fn reference_to_physical_map<'a>(
        &'a self,
        reference_points: &'a [Self::T],
    ) -> Self::ReferenceMap<'a> {
        TriangleReferenceMap::new(reference_points, self)
    }

    fn edge_to_cells(&self, edge_index: usize) -> &[CellLocalIndexPair] {
        self.edge_to_cells[edge_index].as_slice()
    }

    fn point_to_cells(&self, point_index: usize) -> &[CellLocalIndexPair] {
        self.point_to_cells[point_index].as_slice()
    }

    fn face_to_cells(&self, _face_index: usize) -> &[CellLocalIndexPair] {
        std::unimplemented!()
    }
}

#[cfg(test)]
mod test {
    use rlst::rlst_dynamic_array2;
    use rlst_dense::{
        tools::PrettyPrint,
        traits::{RawAccess, RawAccessMut},
    };

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
        for vertex in grid.iter_all_points() {
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

        let mut reference_points = rlst_dynamic_array2!(f64, [2, 2]);

        reference_points[[0, 0]] = 1.0 / 3.0;
        reference_points[[1, 0]] = 1.0 / 3.0;
        reference_points[[0, 1]] = 1.0;
        reference_points[[1, 1]] = 0.0;

        let map = grid.reference_to_physical_map(reference_points.data());
        for cell_index in 0..grid.number_of_cells()
        {
            let mut points = rlst_dynamic_array2!(f64, [3, map.number_of_reference_points()]);
            for point_index in 0..map.number_of_reference_points() {
                map.reference_to_physical(
                    cell_index, point_index,
                    points.view_mut().slice(1, point_index).data_mut(),
                );
            }
            points.pretty_print();
        }
    }
}
