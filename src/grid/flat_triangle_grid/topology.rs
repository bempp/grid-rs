//! Implementation of grid topology

use crate::grid::traits::{Ownership, Topology};
use crate::reference_cell;
use crate::reference_cell::ReferenceCellType;
use crate::types::CellLocalIndexPair;
use std::collections::HashMap;

/// Topology of a serial grid
pub struct SerialFlatTriangleTopology {
    index_map: Vec<usize>,
    entities_to_vertices: Vec<Vec<Vec<usize>>>,
    cells_to_entities: Vec<Vec<Vec<usize>>>,
    entities_to_cells: Vec<Vec<Vec<CellLocalIndexPair<usize>>>>,
    entity_types: Vec<ReferenceCellType>,
}

unsafe impl Sync for SerialFlatTriangleTopology {}

impl SerialFlatTriangleTopology {
    pub fn new(cells: &[usize]) -> Self {
        let ncells = cells.len() / 3;
        let nvertices = *cells.iter().max().unwrap() + 1;

        let mut index_map = vec![0; ncells];

        let entity_types = vec![
            ReferenceCellType::Point,
            ReferenceCellType::Interval,
            ReferenceCellType::Triangle,
        ];

        let mut cells_to_entities = vec![vec![vec![]; ncells]; 3];
        let mut entities_to_cells = vec![vec![]; 3];
        let mut entities_to_vertices = vec![vec![]; 2];

        let mut edge_indices = HashMap::new();

        entities_to_cells[2] = vec![vec![]; ncells];
        entities_to_vertices[0] = (0..nvertices).map(|i| vec![i]).collect::<Vec<_>>();
        entities_to_cells[0] = vec![vec![]; nvertices];

        for (cell_i, i) in index_map.iter_mut().enumerate() {
            let cell = &cells[cell_i * 3..(cell_i + 1) * 3];
            *i = cell_i;
            for (local_index, v) in cell.iter().enumerate() {
                entities_to_cells[0][*v].push(CellLocalIndexPair::new(cell_i, local_index));
            }
            entities_to_cells[2][cell_i] = vec![CellLocalIndexPair::new(cell_i, 0)];
            cells_to_entities[0][cell_i].extend_from_slice(cell);
            cells_to_entities[2][cell_i] = vec![cell_i];
        }

        let ref_conn = &reference_cell::connectivity(ReferenceCellType::Triangle)[1];
        for cell_i in 0..ncells {
            let cell = &cells[cell_i * 3..(cell_i + 1) * 3];
            for (local_index, rc) in ref_conn.iter().enumerate() {
                let mut first = cell[rc[0][0]];
                let mut second = cell[rc[0][1]];
                if first > second {
                    std::mem::swap(&mut first, &mut second);
                }
                if let Some(edge_index) = edge_indices.get(&(first, second)) {
                    cells_to_entities[1][cell_i].push(*edge_index);
                    entities_to_cells[1][*edge_index]
                        .push(CellLocalIndexPair::new(cell_i, local_index));
                } else {
                    edge_indices.insert((first, second), entities_to_vertices[1].len());
                    cells_to_entities[1][cell_i].push(entities_to_vertices[1].len());
                    entities_to_cells[1].push(vec![CellLocalIndexPair::new(cell_i, local_index)]);
                    entities_to_vertices[1].push(vec![first, second]);
                }
            }
        }

        Self {
            index_map,
            entities_to_vertices,
            cells_to_entities,
            entities_to_cells,
            entity_types,
        }
    }
}

impl Topology for SerialFlatTriangleTopology {
    type IndexType = usize;

    fn dim(&self) -> usize {
        2
    }
    fn index_map(&self) -> &[Self::IndexType] {
        &self.index_map
    }
    fn entity_count(&self, etype: ReferenceCellType) -> usize {
        if self.entity_types.contains(&etype) {
            self.entities_to_cells[reference_cell::dim(etype)].len()
        } else {
            0
        }
    }
    fn entity_count_by_dim(&self, dim: usize) -> usize {
        self.entity_count(self.entity_types[dim])
    }
    fn cell(&self, index: Self::IndexType) -> Option<&[usize]> {
        if index < self.cells_to_entities[2].len() {
            Some(&self.cells_to_entities[2][index])
        } else {
            None
        }
    }
    fn cell_type(&self, index: Self::IndexType) -> Option<ReferenceCellType> {
        if index < self.cells_to_entities[2].len() {
            Some(self.entity_types[2])
        } else {
            None
        }
    }

    fn entity_types(&self, dim: usize) -> &[ReferenceCellType] {
        &self.entity_types[dim..dim + 1]
    }

    fn entity_ownership(&self, _dim: usize, _index: Self::IndexType) -> Ownership {
        Ownership::Owned
    }
    fn cell_to_entities(&self, index: Self::IndexType, dim: usize) -> Option<&[Self::IndexType]> {
        if dim <= 2 && index < self.cells_to_entities[dim].len() {
            Some(&self.cells_to_entities[dim][index])
        } else {
            None
        }
    }
    fn entity_to_cells(
        &self,
        dim: usize,
        index: Self::IndexType,
    ) -> Option<&[CellLocalIndexPair<Self::IndexType>]> {
        if dim <= 2 && index < self.entities_to_cells[dim].len() {
            Some(&self.entities_to_cells[dim][index])
        } else {
            None
        }
    }

    fn entity_vertices(&self, dim: usize, index: Self::IndexType) -> Option<&[Self::IndexType]> {
        if dim == 2 {
            self.cell_to_entities(index, 0)
        } else if dim < 2 && index < self.entities_to_vertices[dim].len() {
            Some(&self.entities_to_vertices[dim][index])
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn example_topology() -> SerialFlatTriangleTopology {
        SerialFlatTriangleTopology::new(&[0, 1, 2, 2, 1, 3])
    }

    #[test]
    fn test_counts() {
        let t = example_topology();
        assert_eq!(t.dim(), 2);
        assert_eq!(t.entity_count(ReferenceCellType::Point), 4);
        assert_eq!(t.entity_count(ReferenceCellType::Interval), 5);
        assert_eq!(t.entity_count(ReferenceCellType::Triangle), 2);
    }

    #[test]
    fn test_cell_entities_points() {
        let t = example_topology();
        for (i, vertices) in [[0, 1, 2], [2, 1, 3]].iter().enumerate() {
            let c = t.cell_to_entities(i, 0).unwrap();
            assert_eq!(c.len(), 3);
            assert_eq!(c, vertices);
        }
    }

    #[test]
    fn test_cell_entities_intervals() {
        let t = example_topology();
        for (i, edges) in [[0, 1, 2], [3, 4, 0]].iter().enumerate() {
            let c = t.cell_to_entities(i, 1).unwrap();
            assert_eq!(c.len(), 3);
            assert_eq!(c, edges);
        }
    }

    #[test]
    fn test_cell_entities_triangles() {
        let t = example_topology();
        for i in 0..2 {
            let c = t.cell_to_entities(i, 2).unwrap();
            assert_eq!(c.len(), 1);
            assert_eq!(c[0], i);
        }
    }

    #[test]
    fn test_entities_to_cells_points() {
        let t = example_topology();
        let c_to_e = (0..t.entity_count(ReferenceCellType::Triangle))
            .map(|i| t.cell_to_entities(i, 0).unwrap())
            .collect::<Vec<_>>();
        let e_to_c = (0..t.entity_count(ReferenceCellType::Point))
            .map(|i| {
                t.entity_to_cells(0, i)
                    .unwrap()
                    .iter()
                    .map(|x| x.cell)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        for (i, cell) in c_to_e.iter().enumerate() {
            for v in *cell {
                assert!(e_to_c[*v].contains(&i));
            }
        }
        for (i, cells) in e_to_c.iter().enumerate() {
            for c in cells {
                assert!(c_to_e[*c].contains(&i));
            }
        }
    }

    #[test]
    fn test_entities_to_cells_edges() {
        let t = example_topology();
        let c_to_e = (0..t.entity_count(ReferenceCellType::Triangle))
            .map(|i| t.cell_to_entities(i, 1).unwrap())
            .collect::<Vec<_>>();
        let e_to_c = (0..t.entity_count(ReferenceCellType::Interval))
            .map(|i| {
                t.entity_to_cells(1, i)
                    .unwrap()
                    .iter()
                    .map(|x| x.cell)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        for (i, cell) in c_to_e.iter().enumerate() {
            for v in *cell {
                assert!(e_to_c[*v].contains(&i));
            }
        }
        for (i, cells) in e_to_c.iter().enumerate() {
            for c in cells {
                assert!(c_to_e[*c].contains(&i));
            }
        }
    }
}
