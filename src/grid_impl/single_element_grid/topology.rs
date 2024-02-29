//! Implementation of grid topology

use crate::grid_impl::traits::{Ownership, Topology};
use crate::reference_cell;
use crate::reference_cell::ReferenceCellType;

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

/// Topology of a serial grid
pub struct SerialSingleElementTopology {
    dim: usize,
    index_map: Vec<usize>,
    cells: Vec<Vec<usize>>, // TODO: use 2D array
    connectivity: Vec<Vec<Vec<Vec<usize>>>>,
    cell_connectivity: Vec<Vec<Vec<usize>>>,
    entity_types: Vec<ReferenceCellType>,
}

unsafe impl Sync for SerialSingleElementTopology {}

impl SerialSingleElementTopology {
    pub fn new(cells_input: &[usize], cell_type: ReferenceCellType) -> Self {
        let size = reference_cell::entity_counts(cell_type)[0];
        let ncells = cells_input.len() / size;

        let mut index_map = vec![0; ncells];
        let mut vertices = vec![];
        let dim = reference_cell::dim(cell_type);

        let entity_types = reference_cell::entity_types(cell_type)
            .iter()
            .filter(|t| !t.is_empty())
            .map(|t| t[0])
            .collect::<Vec<_>>();
        let mut cells = vec![];
        let mut connectivity = vec![vec![vec![]; dim + 1]; dim + 1];
        let mut cell_connectivity = vec![vec![]; dim + 1];

        // dim0 = dim, dim1 = 0
        let mut cty = vec![];
        let mut cell_cty = vec![vec![]; dim + 1];
        let mut start = 0;
        for (cell_i, i) in index_map.iter_mut().enumerate() {
            let cell = &cells_input[start..start + size];
            *i = cell_i;
            let mut row = vec![];
            for v in cell {
                if !vertices.contains(v) {
                    vertices.push(*v);
                }
                row.push(vertices.iter().position(|&r| r == *v).unwrap());
            }
            cell_cty[0].extend_from_slice(&row);
            cells.push(row.clone());
            cty.push(row);

            start += size;
        }
        connectivity[dim][0] = cty;
        cell_connectivity[dim] = cell_cty;
        // dim1 == 0
        for dim0 in 0..dim {
            let mut cty: Vec<Vec<usize>> = vec![];
            let ref_conn = &reference_cell::connectivity(cell_type)[dim0];
            for cell in &connectivity[dim][0] {
                for rc in ref_conn {
                    let vertices = rc[0].iter().map(|x| cell[*x]).collect::<Vec<_>>();
                    let mut found = false;
                    for entity in &cty {
                        if all_equal(entity, &vertices) {
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        cty.push(vertices);
                    }
                }
            }
            connectivity[dim0][0] = cty;
        }
        // dim0 == dim1 == 0
        let mut nvertices = 0;
        let mut cty = vec![];
        for cell in &connectivity[dim][0] {
            nvertices = std::cmp::max(nvertices, 1 + cell.iter().max().unwrap());
        }
        for i in 0..nvertices {
            cty.push(vec![i]);
        }
        connectivity[0][0] = cty;
        // dim0 == dim1
        for (dim0, c) in connectivity.iter_mut().enumerate() {
            let mut cty = vec![];
            for i in 0..c[0].len() {
                cty.push(vec![i]);
            }
            c[dim0] = cty;
        }
        // dim0 == dim1 == dim
        let mut cty = vec![];
        for i in 0..connectivity[dim][0].len() {
            cty.push(i);
        }
        cell_connectivity[dim][dim] = cty;
        // dim0 == dim
        for dim1 in 1..dim {
            let mut cty = vec![];
            let mut cell_cty = vec![];
            let entities0 = &connectivity[dim][0];
            let ref_conn = &reference_cell::connectivity(cell_type)[dim1];

            let entities1 = &connectivity[dim1][0];

            for entity0 in entities0 {
                let mut row = vec![];
                for rc in ref_conn {
                    let vertices = rc[0].iter().map(|x| entity0[*x]).collect::<Vec<_>>();
                    for (j, entity1) in entities1.iter().enumerate() {
                        if all_equal(&vertices, entity1) {
                            row.push(j);
                            break;
                        }
                    }
                }
                cell_cty.extend_from_slice(&row);
                cty.push(row);
            }
            cell_connectivity[dim][dim1] = cell_cty;
            connectivity[dim][dim1] = cty;
        }
        // dim1 < dim0
        for (dim0, etype0) in entity_types.iter().enumerate().take(dim + 1).skip(2) {
            for dim1 in 1..dim0 {
                let mut cty = vec![];
                let entities0 = &connectivity[dim0][0];
                let ref_conn = &reference_cell::connectivity(*etype0)[dim1];
                let entities1 = &connectivity[dim1][0];
                for entity0 in entities0 {
                    let mut row = vec![];
                    for rc in ref_conn {
                        let vertices = rc[0].iter().map(|x| entity0[*x]).collect::<Vec<_>>();
                        for (j, entity1) in entities1.iter().enumerate() {
                            if all_equal(&vertices, entity1) {
                                row.push(j);
                                break;
                            }
                        }
                    }
                    cty.push(row);
                }
                connectivity[dim0][dim1] = cty;
            }
        }
        // dim1 > dim0
        for dim0 in 0..dim {
            for dim1 in 1..dim + 1 {
                let mut data = vec![vec![]; connectivity[dim0][0].len()];
                for (i, row) in connectivity[dim1][dim0].iter().enumerate() {
                    for v in row {
                        data[*v].push(i);
                    }
                }
                connectivity[dim0][dim1] = data;
            }
        }
        Self {
            dim,
            index_map,
            cells,
            connectivity,
            cell_connectivity,
            entity_types,
        }
    }
}

impl Topology for SerialSingleElementTopology {
    type IndexType = usize;

    fn dim(&self) -> usize {
        self.dim
    }
    fn index_map(&self) -> &[Self::IndexType] {
        &self.index_map
    }
    fn entity_count(&self, etype: ReferenceCellType) -> usize {
        if self.entity_types.contains(&etype) {
            self.connectivity[reference_cell::dim(etype)][0].len()
        } else {
            0
        }
    }
    fn entity_count_by_dim(&self, dim: usize) -> usize {
        self.entity_count(self.entity_types[dim])
    }
    fn cell(&self, index: Self::IndexType) -> Option<&[usize]> {
        if index < self.cells.len() {
            Some(&self.cells[index])
        } else {
            None
        }
    }
    fn cell_type(&self, index: Self::IndexType) -> Option<ReferenceCellType> {
        if index < self.cells.len() {
            Some(self.entity_types[self.dim])
        } else {
            None
        }
    }

    fn entity_types(&self, dim: usize) -> &[ReferenceCellType] {
        &self.entity_types[dim..dim + 1]
    }

    fn cell_entities(&self, cell_type: ReferenceCellType, dim: usize) -> Option<&[usize]> {
        if self.entity_types.contains(&cell_type)
            && dim < self.cell_connectivity[reference_cell::dim(cell_type)].len()
        {
            Some(&self.cell_connectivity[reference_cell::dim(cell_type)][dim])
        } else {
            None
        }
    }

    fn connectivity(&self, dim0: usize, index: usize, dim1: usize) -> Option<&[Self::IndexType]> {
        if dim0 < self.connectivity.len()
            && dim1 < self.connectivity[dim0].len()
            && index < self.connectivity[dim0][dim1].len()
        {
            Some(&self.connectivity[dim0][dim1][index])
        } else {
            None
        }
    }

    fn entity_ownership(&self, _dim: usize, _index: Self::IndexType) -> Ownership {
        Ownership::Owned
    }

    fn extract_index(&self, index: usize) -> usize {
        index
    }
}

#[cfg(test)]
mod test {
    use crate::grid_impl::single_element_grid::topology::*;

    fn example_topology() -> SerialSingleElementTopology {
        SerialSingleElementTopology::new(&[0, 1, 2, 2, 1, 3], ReferenceCellType::Triangle)
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
        let c = t.cell_entities(ReferenceCellType::Triangle, 0).unwrap();
        assert_eq!(c.len(), 6);
        // cell 0
        assert_eq!(c[0], 0);
        assert_eq!(c[1], 1);
        assert_eq!(c[2], 2);
        // cell 1
        assert_eq!(c[3], 2);
        assert_eq!(c[4], 1);
        assert_eq!(c[5], 3);
    }

    #[test]
    fn test_cell_entities_intervals() {
        let t = example_topology();
        let c = t.cell_entities(ReferenceCellType::Triangle, 1).unwrap();
        assert_eq!(c.len(), 6);
        // cell 0
        assert_eq!(c[0], 0);
        assert_eq!(c[1], 1);
        assert_eq!(c[2], 2);
        // cell 1
        assert_eq!(c[3], 3);
        assert_eq!(c[4], 4);
        assert_eq!(c[5], 0);
    }
    #[test]
    fn test_cell_entities_triangles() {
        let t = example_topology();
        let c = t.cell_entities(ReferenceCellType::Triangle, 2).unwrap();
        assert_eq!(c.len(), 2);
        // cell 0
        assert_eq!(c[0], 0);
        // cell 1
        assert_eq!(c[1], 1);
    }
    #[test]
    fn test_vertex_connectivity() {
        let t = example_topology();

        for (id, vertices) in [vec![0], vec![1], vec![2], vec![3]].iter().enumerate() {
            let c = t.connectivity(0, id, 0).unwrap();
            for (i, j) in c.iter().zip(vertices) {
                assert_eq!(*i, *j);
            }
        }

        for (id, edges) in [vec![1, 2], vec![0, 2, 3], vec![0, 1, 4], vec![3, 4]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(0, id, 1).unwrap();
            for (i, j) in c.iter().zip(edges) {
                assert_eq!(*i, *j);
            }
        }

        for (id, faces) in [vec![0], vec![0, 1], vec![0, 1], vec![1]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(0, id, 2).unwrap();
            for (i, j) in c.iter().zip(faces) {
                assert_eq!(*i, *j);
            }
        }
    }
    #[test]
    fn test_edge_connectivity() {
        let t = example_topology();

        for (id, vertices) in [vec![1, 2], vec![0, 2], vec![0, 1], vec![1, 3], vec![2, 3]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(1, id, 0).unwrap();
            for (i, j) in c.iter().zip(vertices) {
                assert_eq!(*i, *j);
            }
        }

        for (id, edges) in [vec![0], vec![1], vec![2], vec![3], vec![4]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(1, id, 1).unwrap();
            for (i, j) in c.iter().zip(edges) {
                assert_eq!(*i, *j);
            }
        }

        for (id, faces) in [vec![0, 1], vec![0], vec![0], vec![1], vec![1]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(1, id, 2).unwrap();
            for (i, j) in c.iter().zip(faces) {
                assert_eq!(*i, *j);
            }
        }
    }

    #[test]
    fn test_cell_entities_vs_connectivity() {
        let t = example_topology();
        for cell_type in t.entity_types(t.dim()) {
            for dim in 0..t.dim() + 1 {
                let ce = t.cell_entities(*cell_type, dim).unwrap();
                let n = reference_cell::entity_counts(*cell_type)[dim];
                for i in 0..ce.len() / n {
                    // TODO: entity_count instead?
                    println!("{cell_type:?} {dim} {i}");
                    let con = t.connectivity(t.dim(), i, dim).unwrap();
                    println!(
                        "{cell_type:?} {dim} {i} -> {:?} {:?}",
                        con,
                        &ce[n * i..n * (i + 1)]
                    );
                }
            }
        }
        for cell_type in t.entity_types(t.dim()) {
            for dim in 0..t.dim() + 1 {
                let ce = t.cell_entities(*cell_type, dim).unwrap();
                let n = reference_cell::entity_counts(*cell_type)[dim];
                for i in 0..ce.len() / n {
                    // TODO: entity_count instead?
                    let con = t.connectivity(t.dim(), i, dim).unwrap();
                    assert_eq!(con, &ce[n * i..n * (i + 1)]);
                }
            }
        }
    }
}
