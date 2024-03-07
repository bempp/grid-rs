//! Implementation of grid topology

use crate::grid::traits::{Ownership, Topology};
use crate::reference_cell;
use crate::reference_cell::ReferenceCellType;
use crate::types::CellLocalIndexPair;

use std::collections::HashMap;

fn all_equal<T: Eq>(a: &[T], b: &[T]) -> bool {
    if a.len() != b.len() {
        false
    } else {
        all_in(a, b)
    }
}

fn all_in<T: Eq>(a: &[T], b: &[T]) -> bool {
    for i in a {
        if !b.contains(i) {
            return false;
        }
    }
    true
}

type IndexType = (ReferenceCellType, usize);

/// Topology of a serial grid
pub struct SerialMixedTopology {
    dim: usize,
    index_map: Vec<IndexType>,
    entities_to_vertices: Vec<HashMap<ReferenceCellType, Vec<Vec<IndexType>>>>,
    cells_to_entities: Vec<HashMap<ReferenceCellType, Vec<Vec<IndexType>>>>,
    entities_to_cells: Vec<HashMap<ReferenceCellType, Vec<Vec<CellLocalIndexPair<IndexType>>>>>,
    entity_types: Vec<Vec<ReferenceCellType>>,
}

unsafe impl Sync for SerialMixedTopology {}

impl SerialMixedTopology {
    pub fn new(cells_input: &[usize], cell_types: &[ReferenceCellType]) -> Self {
        let mut index_map = vec![(ReferenceCellType::Point, 0); cell_types.len()];
        let mut vertices = vec![];
        let dim = reference_cell::dim(cell_types[0]);

        let mut entity_types = vec![vec![]; 4];

        let mut entities_to_vertices = vec![HashMap::new(); dim];
        let mut cells_to_entities = vec![HashMap::new(); dim + 1];
        let mut entities_to_cells = vec![HashMap::new(); dim + 1];

        for c in cell_types {
            if dim != reference_cell::dim(*c) {
                panic!("Grids with cells of mixed topological dimension not supported.");
            }
            for (dim0, etypes) in reference_cell::entity_types(*c).iter().enumerate() {
                for e in etypes {
                    if !entity_types[dim0].contains(e) {
                        entity_types[dim0].push(*e);

                        entities_to_cells[dim0].insert(*e, vec![]);
                        if dim0 == dim {
                            for ce in cells_to_entities.iter_mut() {
                                ce.insert(*e, vec![]);
                            }
                        } else {
                            entities_to_vertices[dim0].insert(*e, vec![]);
                        }
                    }
                }
            }
        }

        // dim0 = dim, dim1 = 0
        for c in &entity_types[dim] {
            let n = reference_cell::entity_counts(*c)[0];
            let mut start = 0;
            for (i, ct) in cell_types.iter().enumerate() {
                if *ct == *c {
                    let cell = &cells_input[start..start + n];
                    let cell_i = (*c, cells_to_entities[0][c].len());
                    index_map[i] = cell_i;
                    let mut row = vec![];
                    for v in cell {
                        if !vertices.contains(v) {
                            for (_, ec) in entities_to_cells[0].iter_mut() {
                                ec.push(vec![]);
                            }
                            vertices.push(*v);
                        }
                        row.push((
                            ReferenceCellType::Point,
                            vertices.iter().position(|&r| r == *v).unwrap(),
                        ));
                    }

                    for (local_index, v) in row.iter().enumerate() {
                        entities_to_cells[0].get_mut(&v.0).unwrap()[v.1]
                            .push(CellLocalIndexPair::new(cell_i, local_index));
                    }
                    entities_to_cells[dim]
                        .get_mut(c)
                        .unwrap()
                        .push(vec![CellLocalIndexPair::new(cell_i, 0)]);

                    cells_to_entities[0].get_mut(c).unwrap().push(row);
                    cells_to_entities[dim]
                        .get_mut(c)
                        .unwrap()
                        .push(vec![cell_i]);
                }
                start += reference_cell::entity_counts(*ct)[0];
            }
        }
        for i in 0..vertices.len() {
            entities_to_vertices[0]
                .get_mut(&ReferenceCellType::Point)
                .unwrap()
                .push(vec![(ReferenceCellType::Point, i)]);
        }

        for (d, etypes0) in entity_types.iter().enumerate().take(dim).skip(1) {
            for etype in etypes0 {
                for cell_type in &entity_types[dim] {
                    let mut c_to_e = vec![];
                    let ref_conn = &reference_cell::connectivity(*cell_type)[d];
                    for (cell_i, cell) in cells_to_entities[0][&cell_type].iter().enumerate() {
                        let mut entity_ids = vec![];

                        for (local_index, rc) in ref_conn.iter().enumerate() {
                            let vertices = rc[0].iter().map(|x| cell[*x]).collect::<Vec<_>>();
                            let mut found = false;
                            for (entity_index, entity) in
                                entities_to_vertices[d][etype].iter().enumerate()
                            {
                                if all_equal(entity, &vertices) {
                                    entity_ids.push((*etype, entity_index));
                                    entities_to_cells[d].get_mut(etype).unwrap()[entity_index]
                                        .push(CellLocalIndexPair::new(
                                            (*cell_type, cell_i),
                                            local_index,
                                        ));

                                    found = true;
                                    break;
                                }
                            }
                            if !found {
                                entity_ids.push((*etype, entities_to_vertices[d][etype].len()));
                                entities_to_cells[d].get_mut(etype).unwrap().push(vec![
                                    CellLocalIndexPair::new((*cell_type, cell_i), local_index),
                                ]);
                                entities_to_vertices[d]
                                    .get_mut(etype)
                                    .unwrap()
                                    .push(vertices);
                            }
                        }
                        c_to_e.push(entity_ids);
                    }
                    *cells_to_entities[d].get_mut(cell_type).unwrap() = c_to_e;
                }
            }
        }

        Self {
            dim,
            index_map,
            entities_to_vertices,
            cells_to_entities,
            entities_to_cells,
            entity_types,
        }
    }
}

impl Topology for SerialMixedTopology {
    type IndexType = IndexType;

    fn dim(&self) -> usize {
        self.dim
    }
    fn index_map(&self) -> &[Self::IndexType] {
        &self.index_map
    }
    fn entity_count(&self, etype: ReferenceCellType) -> usize {
        let dim = reference_cell::dim(etype);
        if self.entity_types[dim].contains(&etype) {
            self.entities_to_cells[dim][&etype].len()
        } else {
            0
        }
    }
    fn entity_count_by_dim(&self, dim: usize) -> usize {
        self.entity_types[dim]
            .iter()
            .map(|e| self.entity_count(*e))
            .sum()
    }
    fn cell(&self, index: Self::IndexType) -> Option<&[IndexType]> {
        if self.cells_to_entities[0].contains_key(&index.0)
            && index.1 < self.cells_to_entities[0][&index.0].len()
        {
            Some(&self.cells_to_entities[0][&index.0][index.1])
        } else {
            None
        }
    }
    fn cell_type(&self, index: Self::IndexType) -> Option<ReferenceCellType> {
        if self.cells_to_entities[0].contains_key(&index.0)
            && index.1 < self.cells_to_entities[0][&index.0].len()
        {
            Some(index.0)
        } else {
            None
        }
    }

    fn entity_types(&self, dim: usize) -> &[ReferenceCellType] {
        &self.entity_types[dim]
    }

    fn entity_ownership(&self, _dim: usize, _index: Self::IndexType) -> Ownership {
        Ownership::Owned
    }

    fn entity_to_cells(
        &self,
        dim: usize,
        index: Self::IndexType,
    ) -> Option<&[CellLocalIndexPair<Self::IndexType>]> {
        if dim <= self.dim
            && self.entities_to_cells[dim].contains_key(&index.0)
            && index.1 < self.entities_to_cells[dim][&index.0].len()
        {
            Some(&self.entities_to_cells[dim][&index.0][index.1])
        } else {
            None
        }
    }
    fn cell_to_entities(&self, index: Self::IndexType, dim: usize) -> Option<&[Self::IndexType]> {
        if dim <= self.dim
            && self.cells_to_entities[dim].contains_key(&index.0)
            && index.1 < self.cells_to_entities[dim][&index.0].len()
        {
            Some(&self.cells_to_entities[dim][&index.0][index.1])
        } else {
            None
        }
    }

    fn entity_vertices(&self, dim: usize, index: Self::IndexType) -> Option<&[Self::IndexType]> {
        if dim == self.dim {
            self.cell_to_entities(index, 0)
        } else if dim < self.dim
            && self.entities_to_vertices[dim].contains_key(&index.0)
            && index.1 < self.entities_to_vertices[dim][&index.0].len()
        {
            Some(&self.entities_to_vertices[dim][&index.0][index.1])
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn example_topology() -> SerialMixedTopology {
        SerialMixedTopology::new(&[0, 1, 2, 2, 1, 3], &[ReferenceCellType::Triangle; 2])
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
    fn test_cell_to_entities_vertices() {
        let t = example_topology();
        for (i, vertices) in [[0, 1, 2], [2, 1, 3]].iter().enumerate() {
            let c = t
                .cell_to_entities((ReferenceCellType::Triangle, i), 0)
                .unwrap();
            assert_eq!(c.len(), 3);
            for i in c {
                assert_eq!(i.0, ReferenceCellType::Point);
            }
            assert_eq!(c[0].1, vertices[0]);
            assert_eq!(c[1].1, vertices[1]);
            assert_eq!(c[2].1, vertices[2]);
        }
    }
    #[test]
    fn test_cell_to_entities_intervals() {
        let t = example_topology();
        for (i, edges) in [[0, 1, 2], [3, 4, 0]].iter().enumerate() {
            let c = t
                .cell_to_entities((ReferenceCellType::Triangle, i), 1)
                .unwrap();
            assert_eq!(c.len(), 3);
            for i in c {
                assert_eq!(i.0, ReferenceCellType::Interval);
            }
            assert_eq!(c[0].1, edges[0]);
            assert_eq!(c[1].1, edges[1]);
            assert_eq!(c[2].1, edges[2]);
        }
    }
    #[test]
    fn test_cell_to_entities_triangles() {
        let t = example_topology();
        for i in 0..2 {
            let c = t
                .cell_to_entities((ReferenceCellType::Triangle, i), 2)
                .unwrap();
            assert_eq!(c.len(), 1);
            assert_eq!(c[0].0, ReferenceCellType::Triangle);
            assert_eq!(c[0].1, i);
        }
    }

    fn example_topology_mixed() -> SerialMixedTopology {
        SerialMixedTopology::new(
            &[0, 1, 2, 3, 1, 4, 3],
            &[
                ReferenceCellType::Quadrilateral,
                ReferenceCellType::Triangle,
            ],
        )
    }

    #[test]
    fn test_mixed_counts() {
        let t = example_topology_mixed();
        assert_eq!(t.dim(), 2);
        assert_eq!(t.entity_count(ReferenceCellType::Point), 5);
        assert_eq!(t.entity_count(ReferenceCellType::Interval), 6);
        assert_eq!(t.entity_count(ReferenceCellType::Triangle), 1);
        assert_eq!(t.entity_count(ReferenceCellType::Quadrilateral), 1);
    }

    #[test]
    fn test_mixed_cell_entities_points() {
        let t = example_topology_mixed();
        let c = t
            .cell_to_entities((ReferenceCellType::Quadrilateral, 0), 0)
            .unwrap();
        assert_eq!(c.len(), 4);
        for i in c {
            assert_eq!(i.0, ReferenceCellType::Point);
        }
        // cell 0
        assert_eq!(c[0].1, 0);
        assert_eq!(c[1].1, 1);
        assert_eq!(c[2].1, 2);
        assert_eq!(c[3].1, 3);

        let c = t
            .cell_to_entities((ReferenceCellType::Triangle, 0), 0)
            .unwrap();
        assert_eq!(c.len(), 3);
        // cell 1
        assert_eq!(c[0].1, 1);
        assert_eq!(c[1].1, 4);
        assert_eq!(c[2].1, 3);
    }

    #[test]
    fn test_mixed_cell_entities_intervals() {
        let t = example_topology_mixed();
        let c = t
            .cell_to_entities((ReferenceCellType::Quadrilateral, 0), 1)
            .unwrap();

        assert_eq!(c.len(), 4);
        for i in c {
            assert_eq!(i.0, ReferenceCellType::Interval);
        }
        // cell 0
        assert_eq!(c[0].1, 0);
        assert_eq!(c[1].1, 1);
        assert_eq!(c[2].1, 2);
        assert_eq!(c[3].1, 3);

        let c = t
            .cell_to_entities((ReferenceCellType::Triangle, 0), 1)
            .unwrap();
        assert_eq!(c.len(), 3);
        // cell 1
        assert_eq!(c[0].1, 4);
        assert_eq!(c[1].1, 2);
        assert_eq!(c[2].1, 5);
    }
    #[test]
    fn test_mixed_cell_entities_triangles() {
        let t = example_topology_mixed();
        let c = t
            .cell_to_entities((ReferenceCellType::Quadrilateral, 0), 2)
            .unwrap();
        assert_eq!(c.len(), 1);
        // cell 0
        assert_eq!(c[0].0, ReferenceCellType::Quadrilateral);
        assert_eq!(c[0].1, 0);

        let c = t
            .cell_to_entities((ReferenceCellType::Triangle, 0), 2)
            .unwrap();
        assert_eq!(c.len(), 1);
        // cell 1
        assert_eq!(c[0].0, ReferenceCellType::Triangle);
        assert_eq!(c[0].1, 0);
    }
}
