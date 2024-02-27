//! Implementation of grid topology

use crate::grid::traits::{Ownership, Topology};
use crate::reference_cell;
use crate::traits::cell::ReferenceCellType;

use std::collections::HashMap;

/// Topology of a serial grid
pub struct SerialTopology {
    dim: usize,
    connectivity: HashMap<ReferenceCellType, Vec<Vec<Vec<usize>>>>,
    cell_connectivity: HashMap<ReferenceCellType, HashMap<ReferenceCellType, Vec<usize>>>,
    index_map: Vec<usize>,
    entity_types: Vec<Vec<ReferenceCellType>>,
}

unsafe impl Sync for SerialTopology {}

impl SerialTopology {
    pub fn new(cells: &[usize], cell_types: &[ReferenceCellType]) -> Self {
        let mut index_map = vec![];
        let mut vertices = vec![];
        let dim = reference_cell::dim(cell_types[0]);

        let mut entity_types = vec![vec![]; 4];
        let mut connectivity = HashMap::new();
        for c in cell_types {
            if dim != reference_cell::dim(*c) {
                panic!("Grids with cells of mixed topological dimension not supported.");
            }
            for (dim0, etypes) in reference_cell::entity_types(*c).iter().enumerate() {
                for e in etypes {
                    if !entity_types[dim0].contains(e) {
                        entity_types[dim0].push(*e);
                        connectivity.insert(*e, vec![vec![]; dim + 1]);
                    }
                }
            }
        }

        // dim0 = dim, dim1 = 0
        for c in &entity_types[dim] {
            let mut cty = vec![];
            let n = reference_cell::entity_counts(*c)[0];
            let mut start = 0;
            for (i, ct) in cell_types.iter().enumerate() {
                if *ct == *c {
                    let cell = &cells[start..start + n];
                    index_map.push(i);
                    let mut row = vec![];
                    for v in cell {
                        if !vertices.contains(v) {
                            vertices.push(*v);
                        }
                        row.push(vertices.iter().position(|&r| r == *v).unwrap());
                    }
                    cty.push(row);
                }
                start += reference_cell::entity_counts(*ct)[0];
            }
            connectivity.get_mut(c).unwrap()[0] = cty;
        }

        // dim1 == 0
        for (dim0, etypes0) in entity_types.iter().enumerate().take(dim) {
            for etype in etypes0 {
                let mut cty: Vec<Vec<usize>> = vec![];
                for cell_type in &entity_types[dim] {
                    let ref_conn = reference_cell::connectivity(*cell_type);
                    let ref_entities = (0..reference_cell::entity_counts(*cell_type)[dim0])
                        .map(|x| ref_conn[dim0][x][0].clone())
                        .collect::<Vec<Vec<usize>>>();

                    let cells = &connectivity[&cell_type][0];
                    for cell in cells {
                        for e in &ref_entities {
                            let vertices = e.iter().map(|x| cell[*x]).collect::<Vec<usize>>();
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
                }
                connectivity.get_mut(etype).unwrap()[0] = cty;
            }
        }

        // dim0 == dim1 == 0
        let mut nvertices = 0;
        let mut cty = vec![];
        for cell_type in &entity_types[dim] {
            let cells = &connectivity[&cell_type][0];
            for cell in cells {
                for j in cell {
                    if *j >= nvertices {
                        nvertices = *j + 1;
                    }
                }
            }
        }
        for i in 0..nvertices {
            cty.push(vec![i]);
        }
        connectivity.get_mut(&ReferenceCellType::Point).unwrap()[0] = cty;

        // dim0 == dim1
        for (dim0, etypes) in entity_types.iter().enumerate().skip(1) {
            let mut start = 0;
            for etype in etypes {
                let mut cty = vec![];
                for i in 0..connectivity[etype][0].len() {
                    cty.push(vec![start + i]);
                }
                start += connectivity[etype][0].len();
                connectivity.get_mut(etype).unwrap()[dim0] = cty;
            }
        }

        // dim0 == dim
        for cell_type in &entity_types[dim] {
            for (dim1, etypes0) in entity_types.iter().enumerate().take(dim).skip(1) {
                let mut cty = vec![];
                let entities0 = &connectivity[cell_type][0];
                let mut start = 0;
                for etype in etypes0 {
                    let entities1 = &connectivity[etype][0];

                    for entity0 in entities0 {
                        let mut row = vec![];
                        for i in 0..reference_cell::entity_counts(*cell_type)[dim1] {
                            let vertices = reference_cell::connectivity(*cell_type)[dim1][i][0]
                                .iter()
                                .map(|x| entity0[*x])
                                .collect::<Vec<usize>>();
                            for (j, entity1) in entities1.iter().enumerate() {
                                if all_equal(&vertices, entity1) {
                                    row.push(start + j);
                                    break;
                                }
                            }
                        }
                        cty.push(row);
                    }
                    start += entities1.len();
                }
                connectivity.get_mut(cell_type).unwrap()[dim1] = cty;
            }
        }

        // dim1 < dim0
        for (dim0, etypes0) in entity_types.iter().enumerate().take(dim + 1).skip(2) {
            for etype0 in etypes0 {
                for (dim1, etypes1) in entity_types.iter().enumerate().take(dim0).skip(1) {
                    let mut cty = vec![];
                    let entities0 = &connectivity[etype0][0];
                    let mut start1 = 0;
                    for etype1 in etypes1 {
                        let entities1 = &connectivity[etype1][0];
                        for entity0 in entities0 {
                            let mut row = vec![];
                            for i in 0..reference_cell::entity_counts(*etype0)[dim1] {
                                let vertices = reference_cell::connectivity(*etype0)[dim1][i][0]
                                    .iter()
                                    .map(|x| entity0[*x])
                                    .collect::<Vec<usize>>();
                                for (j, entity1) in entities1.iter().enumerate() {
                                    if all_equal(&vertices, entity1) {
                                        row.push(start1 + j);
                                        break;
                                    }
                                }
                            }
                            cty.push(row);
                        }
                        start1 += entities1.len();
                    }
                    connectivity.get_mut(etype0).unwrap()[dim1] = cty;
                }
            }
        }

        // dim1 > dim0
        for (dim0, etypes0) in entity_types.iter().enumerate().take(dim) {
            for etype0 in etypes0 {
                for (dim1, etypes1) in entity_types.iter().enumerate().take(dim + 1).skip(1) {
                    let mut data = vec![vec![]; connectivity[etype0][0].len()];
                    let mut start = 0;
                    for etype1 in etypes1 {
                        for (i, row) in connectivity[etype1][dim0].iter().enumerate() {
                            for v in row {
                                data[*v].push(start + i);
                            }
                        }
                        start += connectivity[etype1][dim0].len();
                    }
                    connectivity.get_mut(etype0).unwrap()[dim1] = data;
                }
            }
        }

        let mut cell_connectivity = HashMap::new();
        for cell_type in &entity_types[dim] {
            let mut all = HashMap::new();
            for (edim, etypes) in entity_types.iter().enumerate().take(dim + 1) {
                for etype in etypes {
                    all.insert(*etype, vec![]);
                }

                for c in &connectivity[cell_type][edim] {
                    for i in c {
                        let mut start = 0;
                        for etype in &entity_types[edim] {
                            let count = connectivity[etype][0].len();
                            if *i < start + count {
                                all.get_mut(etype).unwrap().push(*i)
                            }
                            start += count;
                        }
                    }
                }
            }
            cell_connectivity.insert(*cell_type, all);
        }

        Self {
            dim,
            connectivity,
            cell_connectivity,
            index_map,
            entity_types,
        }
    }
}

fn all_equal(a: &[usize], b: &[usize]) -> bool {
    if a.len() != b.len() {
        false
    } else {
        all_in(a, b)
    }
}

fn all_in(a: &[usize], b: &[usize]) -> bool {
    for i in a {
        if !b.contains(i) {
            return false;
        }
    }
    true
}

impl Topology for SerialTopology {
    fn dim(&self) -> usize {
        self.dim
    }
    fn index_map(&self) -> &[usize] {
        &self.index_map
    }
    fn entity_count(&self, etype: ReferenceCellType) -> usize {
        self.connectivity[&etype][0].len()
    }
    fn cell(&self, index: usize) -> Option<&[usize]> {
        let mut start = 0;
        for etype in &self.entity_types[self.dim] {
            let count = self.connectivity[etype][0].len();
            if index < start + count {
                return Some(&self.connectivity[etype][0][index - start]);
            }
            start += count;
        }
        None
    }
    fn cell_type(&self, index: usize) -> Option<ReferenceCellType> {
        let mut start = 0;
        for etype in &self.entity_types[self.dim] {
            let count = self.connectivity[etype][0].len();
            if index < start + count {
                return Some(*etype);
            }
            start += count;
        }
        None
    }

    fn entity_types(&self, dim: usize) -> &[ReferenceCellType] {
        &self.entity_types[dim]
    }

    /// Get the indices of entities of type `etype` that are connected to each cell of type `cell_type`
    fn cell_entities(
        &self,
        cell_type: ReferenceCellType,
        etype: ReferenceCellType,
    ) -> Option<&[usize]> {
        if self.cell_connectivity.contains_key(&cell_type)
            && self.cell_connectivity[&cell_type].contains_key(&etype)
        {
            Some(&self.cell_connectivity[&cell_type][&etype])
        } else {
            None
        }
    }

    /// Get the indices of entities of dimension `dim` that are connected to the entity of type `etype` with index `index`
    fn connectivity(&self, etype: ReferenceCellType, index: usize, dim: usize) -> Option<&[usize]> {
        if self.connectivity.contains_key(&etype)
            && dim < self.connectivity[&etype].len()
            && index < self.connectivity[&etype][dim].len()
        {
            Some(&self.connectivity[&etype][dim][index])
        } else {
            None
        }
    }

    fn entity_ownership(&self, _dim: usize, _index: usize) -> Ownership {
        Ownership::Owned
    }
}

#[cfg(test)]
mod test {
    use crate::grid::topology::*;

    fn example_topology() -> SerialTopology {
        SerialTopology::new(&[0, 1, 2, 2, 1, 3], &[ReferenceCellType::Triangle; 2])
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
        let c = t
            .cell_entities(ReferenceCellType::Triangle, ReferenceCellType::Point)
            .unwrap();
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
        let c = t
            .cell_entities(ReferenceCellType::Triangle, ReferenceCellType::Interval)
            .unwrap();
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
        let c = t
            .cell_entities(ReferenceCellType::Triangle, ReferenceCellType::Triangle)
            .unwrap();
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
            let c = t.connectivity(ReferenceCellType::Point, id, 0).unwrap();
            assert_eq!(c, vertices);
        }

        for (id, edges) in [vec![1, 2], vec![0, 2, 3], vec![0, 1, 4], vec![3, 4]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Point, id, 1).unwrap();
            assert_eq!(c, edges);
        }

        for (id, faces) in [vec![0], vec![0, 1], vec![0, 1], vec![1]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Point, id, 2).unwrap();
            assert_eq!(c, faces);
        }
    }
    #[test]
    fn test_edge_connectivity() {
        let t = example_topology();

        for (id, vertices) in [vec![1, 2], vec![0, 2], vec![0, 1], vec![1, 3], vec![2, 3]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Interval, id, 0).unwrap();
            assert_eq!(c, vertices);
        }

        for (id, edges) in [vec![0], vec![1], vec![2], vec![3], vec![4]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Interval, id, 1).unwrap();
            assert_eq!(c, edges);
        }

        for (id, faces) in [vec![0, 1], vec![0], vec![0], vec![1], vec![1]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Interval, id, 2).unwrap();
            assert_eq!(c, faces);
        }
    }

    fn example_topology_mixed() -> SerialTopology {
        SerialTopology::new(
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
            .cell_entities(ReferenceCellType::Quadrilateral, ReferenceCellType::Point)
            .unwrap();
        assert_eq!(c.len(), 4);
        // cell 0
        assert_eq!(c[0], 0);
        assert_eq!(c[1], 1);
        assert_eq!(c[2], 2);
        assert_eq!(c[3], 3);

        let c = t
            .cell_entities(ReferenceCellType::Triangle, ReferenceCellType::Point)
            .unwrap();
        assert_eq!(c.len(), 3);
        // cell 1
        assert_eq!(c[0], 1);
        assert_eq!(c[1], 4);
        assert_eq!(c[2], 3);
    }

    #[test]
    fn test_mixed_cell_entities_intervals() {
        let t = example_topology_mixed();
        let c = t
            .cell_entities(
                ReferenceCellType::Quadrilateral,
                ReferenceCellType::Interval,
            )
            .unwrap();

        assert_eq!(c.len(), 4);
        // cell 0
        assert_eq!(c[0], 0);
        assert_eq!(c[1], 1);
        assert_eq!(c[2], 2);
        assert_eq!(c[3], 3);

        let c = t
            .cell_entities(ReferenceCellType::Triangle, ReferenceCellType::Interval)
            .unwrap();
        assert_eq!(c.len(), 3);
        // cell 1
        assert_eq!(c[0], 4);
        assert_eq!(c[1], 2);
        assert_eq!(c[2], 5);
    }
    #[test]
    fn test_mixed_cell_entities_triangles() {
        let t = example_topology_mixed();
        let c = t
            .cell_entities(
                ReferenceCellType::Quadrilateral,
                ReferenceCellType::Triangle,
            )
            .unwrap();
        assert_eq!(c.len(), 1);
        // cell 0
        assert_eq!(c[0], 0);

        let c = t
            .cell_entities(ReferenceCellType::Triangle, ReferenceCellType::Triangle)
            .unwrap();
        assert_eq!(c.len(), 1);
        // cell 1
        assert_eq!(c[0], 1);
    }
    #[test]
    fn test_mixed_vertex_connectivity() {
        let t = example_topology_mixed();

        for (id, vertices) in [vec![0], vec![1], vec![2], vec![3], vec![4]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Point, id, 0).unwrap();
            assert_eq!(c, vertices);
        }

        for (id, edges) in [
            vec![0, 1],
            vec![0, 2, 5],
            vec![1, 3],
            vec![2, 3, 4],
            vec![4, 5],
        ]
        .iter()
        .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Point, id, 1).unwrap();
            assert_eq!(c, edges);
        }

        for (id, faces) in [vec![0], vec![0, 1], vec![0], vec![0, 1], vec![1]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Point, id, 2).unwrap();
            assert_eq!(c, faces);
        }
    }
    #[test]
    fn test_mixed_edge_connectivity() {
        let t = example_topology_mixed();

        for (id, vertices) in [
            vec![0, 1],
            vec![0, 2],
            vec![1, 3],
            vec![2, 3],
            vec![4, 3],
            vec![1, 4],
        ]
        .iter()
        .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Interval, id, 0).unwrap();
            assert_eq!(c, vertices);
        }

        for (id, edges) in [vec![0], vec![1], vec![2], vec![3], vec![4], vec![5]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Interval, id, 1).unwrap();
            assert_eq!(c, edges);
        }

        for (id, faces) in [vec![0], vec![0], vec![0, 1], vec![0], vec![1], vec![1]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Interval, id, 2).unwrap();
            assert_eq!(c, faces);
        }
    }

    fn cell_entities_vs_connectivity(t: &impl Topology) {
        for cell_type in t.entity_types(t.dim()) {
            for dim in 0..t.dim() + 1 {
                for entity_type in t.entity_types(dim) {
                    let ce = t.cell_entities(*cell_type, *entity_type).unwrap();
                    if ce.len() > 0 {
                        let n = reference_cell::entity_counts(*cell_type)[dim];
                        for i in 0..t.entity_count(*cell_type) {
                            let con = t.connectivity(*cell_type, i, dim).unwrap();
                            assert_eq!(con, &ce[n * i..n * (i + 1)]);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_cell_entities_vs_connectivity() {
        let t = example_topology();
        cell_entities_vs_connectivity(&t);
    }

    #[test]
    fn test_cell_entities_vs_connectivity_mixes() {
        let t = example_topology_mixed();
        cell_entities_vs_connectivity(&t);
    }
}
