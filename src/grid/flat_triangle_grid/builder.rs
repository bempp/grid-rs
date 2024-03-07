//! Grid builder

use crate::grid::flat_triangle_grid::grid::SerialFlatTriangleGrid;
use crate::traits::builder::Builder;
use bempp_element::element::Inverse;
use num::Float;
use rlst_common::types::Scalar;
use rlst_dense::{rlst_array_from_slice2, rlst_dynamic_array2};
use std::collections::HashMap;

pub struct SerialFlatTriangleGridBuilder<T: Float + Scalar<Real = T> + Inverse> {
    points: Vec<T>,
    cells: Vec<usize>,
    point_indices_to_ids: Vec<usize>,
    cell_indices_to_ids: Vec<usize>,
    point_ids_to_indices: HashMap<usize, usize>,
    cell_ids_to_indices: HashMap<usize, usize>,
}

impl<T: Float + Scalar<Real = T> + Inverse> Builder<3> for SerialFlatTriangleGridBuilder<T> {
    type GridType = SerialFlatTriangleGrid<T>;
    type T = T;
    type CellData = [usize; 3];
    type GridMetadata = ();

    fn new(_data: ()) -> Self {
        Self {
            points: vec![],
            cells: vec![],
            point_indices_to_ids: vec![],
            cell_indices_to_ids: vec![],
            point_ids_to_indices: HashMap::new(),
            cell_ids_to_indices: HashMap::new(),
        }
    }

    fn new_with_capacity(npoints: usize, ncells: usize, _data: ()) -> Self {
        Self {
            points: Vec::with_capacity(npoints * Self::GDIM),
            cells: Vec::with_capacity(ncells * 3),
            point_indices_to_ids: Vec::with_capacity(npoints),
            cell_indices_to_ids: Vec::with_capacity(ncells),
            point_ids_to_indices: HashMap::new(),
            cell_ids_to_indices: HashMap::new(),
        }
    }

    fn add_point(&mut self, id: usize, data: [T; 3]) {
        self.point_ids_to_indices
            .insert(id, self.point_indices_to_ids.len());
        self.point_indices_to_ids.push(id);
        self.points.extend_from_slice(&data);
    }

    fn add_cell(&mut self, id: usize, cell_data: [usize; 3]) {
        self.cell_ids_to_indices
            .insert(id, self.cell_indices_to_ids.len());
        self.cell_indices_to_ids.push(id);
        self.cells.extend_from_slice(&cell_data);
    }

    fn create_grid(&self) -> Self::GridType {
        // TODO: remove this transposing
        let npts = self.point_indices_to_ids.len();
        let mut points = rlst_dynamic_array2!(T, [npts, 3]);
        points.fill_from(rlst_array_from_slice2!(
            T,
            &self.points,
            [npts, 3],
            [1, npts]
        ));
        SerialFlatTriangleGrid::new(points, &self.cells)
    }
}
