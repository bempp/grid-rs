//! Implementation of grid topology

use crate::grid::traits::{Ownership, Topology};
use crate::reference_cell;
use crate::traits::cell::ReferenceCellType;

use std::collections::HashMap;

type Connectivity = Vec<Vec<(ReferenceCellType, usize)>>;

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
pub struct SerialTopology {
    dim: usize,
    cells: Vec<(ReferenceCellType, Vec<usize>)>,
    connectivity: HashMap<ReferenceCellType, Vec<Connectivity>>,
    cell_connectivity: HashMap<ReferenceCellType, Connectivity>,
    index_map: Vec<usize>,
    entity_types: Vec<Vec<ReferenceCellType>>,
}

unsafe impl Sync for SerialTopology {}

impl SerialTopology {
    pub fn new(cells_input: &[usize], cell_types: &[ReferenceCellType]) -> Self {
        let mut index_map = vec![];
        let mut vertices = vec![];
        let dim = reference_cell::dim(cell_types[0]);

        let mut entity_types = vec![vec![]; 4];
        let mut cells = vec![];
        let mut connectivity = HashMap::new();
        let mut cell_connectivity = HashMap::new();

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
            let mut cty_neww = vec![vec![]; dim + 1];
            let n = reference_cell::entity_counts(*c)[0];
            let mut start = 0;
            for (i, ct) in cell_types.iter().enumerate() {
                if *ct == *c {
                    let cell = &cells_input[start..start + n];
                    index_map.push(i);
                    let mut row = vec![];
                    for v in cell {
                        if !vertices.contains(v) {
                            vertices.push(*v);
                        }
                        row.push((
                            ReferenceCellType::Point,
                            vertices.iter().position(|&r| r == *v).unwrap(),
                        ));
                    }
                    cty_neww[0].extend_from_slice(&row);
                    cells.push((*c, row.iter().map(|x| x.1).collect()));
                    cty.push(row);
                }
                start += reference_cell::entity_counts(*ct)[0];
            }
            connectivity.get_mut(c).unwrap()[0] = cty;
            cell_connectivity.insert(*c, cty_neww);
        }

        // dim1 == 0
        for (dim0, etypes0) in entity_types.iter().enumerate().take(dim) {
            for etype in etypes0 {
                let mut cty: Vec<Vec<(ReferenceCellType, usize)>> = vec![];
                for cell_type in &entity_types[dim] {
                    let ref_conn = reference_cell::connectivity(*cell_type);
                    let ref_entities = (0..reference_cell::entity_counts(*cell_type)[dim0])
                        .map(|x| ref_conn[dim0][x][0].clone())
                        .collect::<Vec<Vec<usize>>>();

                    for cell in &connectivity[&cell_type][0] {
                        for e in &ref_entities {
                            let vertices =
                                e.iter()
                                    .map(|x| cell[*x])
                                    .collect::<Vec<(ReferenceCellType, usize)>>();
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
            for cell in &connectivity[&cell_type][0] {
                for j in cell {
                    if j.1 >= nvertices {
                        nvertices = j.1 + 1;
                    }
                }
            }
        }
        for i in 0..nvertices {
            cty.push(vec![(ReferenceCellType::Point, i)]);
        }
        connectivity.get_mut(&ReferenceCellType::Point).unwrap()[0] = cty;

        // dim0 == dim1
        for (dim0, etypes) in entity_types.iter().enumerate().skip(1) {
            for etype in etypes {
                let mut cty = vec![];
                for i in 0..connectivity[etype][0].len() {
                    cty.push(vec![(*etype, i)]);
                }
                connectivity.get_mut(etype).unwrap()[dim0] = cty;
            }
        }
        // dim0 == dim1 == dim
        for etype in &entity_types[dim] {
            let mut cty = vec![];
            for i in 0..connectivity[etype][0].len() {
                cty.push((*etype, i));
            }
            cell_connectivity.get_mut(etype).unwrap()[dim] = cty;
        }
        // dim0 == dim
        for cell_type in &entity_types[dim] {
            for (dim1, etypes0) in entity_types.iter().enumerate().take(dim).skip(1) {
                let mut cty = vec![];
                let mut cty_neww = vec![];
                let entities0 = &connectivity[cell_type][0];
                for etype in etypes0 {
                    let entities1 = &connectivity[etype][0];

                    for entity0 in entities0 {
                        let mut row = vec![];
                        for i in 0..reference_cell::entity_counts(*cell_type)[dim1] {
                            let vertices = reference_cell::connectivity(*cell_type)[dim1][i][0]
                                .iter()
                                .map(|x| entity0[*x])
                                .collect::<Vec<(ReferenceCellType, usize)>>();
                            for (j, entity1) in entities1.iter().enumerate() {
                                if all_equal(&vertices, entity1) {
                                    row.push((*etype, j));
                                    break;
                                }
                            }
                        }
                        cty_neww.extend_from_slice(&row);
                        cty.push(row);
                    }
                }
                cell_connectivity.get_mut(cell_type).unwrap()[dim1] = cty_neww;
                connectivity.get_mut(cell_type).unwrap()[dim1] = cty;
            }
        }
        // dim1 < dim0
        for (dim0, etypes0) in entity_types.iter().enumerate().take(dim + 1).skip(2) {
            for etype0 in etypes0 {
                for (dim1, etypes1) in entity_types.iter().enumerate().take(dim0).skip(1) {
                    let mut cty = vec![];
                    let entities0 = &connectivity[etype0][0];
                    for etype1 in etypes1 {
                        let entities1 = &connectivity[etype1][0];
                        for entity0 in entities0 {
                            let mut row = vec![];
                            for i in 0..reference_cell::entity_counts(*etype0)[dim1] {
                                let vertices = reference_cell::connectivity(*etype0)[dim1][i][0]
                                    .iter()
                                    .map(|x| entity0[*x])
                                    .collect::<Vec<(ReferenceCellType, usize)>>();
                                for (j, entity1) in entities1.iter().enumerate() {
                                    if all_equal(&vertices, entity1) {
                                        row.push((*etype1, j));
                                        break;
                                    }
                                }
                            }
                            cty.push(row);
                        }
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
                    for etype1 in etypes1 {
                        for (i, row) in connectivity[etype1][dim0].iter().enumerate() {
                            for v in row {
                                data[v.1].push((*etype1, i));
                            }
                        }
                    }
                    connectivity.get_mut(etype0).unwrap()[dim1] = data;
                }
            }
        }

        Self {
            dim,
            cells,
            connectivity,
            cell_connectivity,
            index_map,
            entity_types,
        }
    }
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
        if index < self.cells.len() {
            Some(&self.cells[index].1)
        } else {
            None
        }
    }
    fn cell_type(&self, index: usize) -> Option<ReferenceCellType> {
        if index < self.cells.len() {
            Some(self.cells[index].0)
        } else {
            None
        }
    }

    fn entity_types(&self, dim: usize) -> &[ReferenceCellType] {
        &self.entity_types[dim]
    }

    fn cell_entities(
        &self,
        cell_type: ReferenceCellType,
        dim: usize,
    ) -> Option<&[(ReferenceCellType, usize)]> {
        if self.cell_connectivity.contains_key(&cell_type)
            && dim < self.cell_connectivity[&cell_type].len()
        {
            Some(&self.cell_connectivity[&cell_type][dim])
        } else {
            None
        }
    }

    fn connectivity(
        &self,
        etype: ReferenceCellType,
        index: usize,
        dim: usize,
    ) -> Option<&[(ReferenceCellType, usize)]> {
        if self.cell_connectivity.contains_key(&etype) {
            let n = reference_cell::entity_counts(etype)[dim];
            Some(&self.cell_connectivity[&etype][dim][n*index..n*(index + 1)])
        } else if self.connectivity.contains_key(&etype)
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
        let c = t.cell_entities(ReferenceCellType::Triangle, 0).unwrap();
        assert_eq!(c.len(), 6);
        for i in c {
            assert_eq!(i.0, ReferenceCellType::Point);
        }
        // cell 0
        assert_eq!(c[0].1, 0);
        assert_eq!(c[1].1, 1);
        assert_eq!(c[2].1, 2);
        // cell 1
        assert_eq!(c[3].1, 2);
        assert_eq!(c[4].1, 1);
        assert_eq!(c[5].1, 3);
    }

    #[test]
    fn test_cell_entities_intervals() {
        let t = example_topology();
        let c = t.cell_entities(ReferenceCellType::Triangle, 1).unwrap();
        assert_eq!(c.len(), 6);
        for i in c {
            assert_eq!(i.0, ReferenceCellType::Interval);
        }
        // cell 0
        assert_eq!(c[0].1, 0);
        assert_eq!(c[1].1, 1);
        assert_eq!(c[2].1, 2);
        // cell 1
        assert_eq!(c[3].1, 3);
        assert_eq!(c[4].1, 4);
        assert_eq!(c[5].1, 0);
    }
    #[test]
    fn test_cell_entities_triangles() {
        let t = example_topology();
        let c = t.cell_entities(ReferenceCellType::Triangle, 2).unwrap();
        assert_eq!(c.len(), 2);
        // cell 0
        assert_eq!(c[0].0, ReferenceCellType::Triangle);
        assert_eq!(c[0].1, 0);
        // cell 1
        assert_eq!(c[1].0, ReferenceCellType::Triangle);
        assert_eq!(c[1].1, 1);
    }
    #[test]
    fn test_vertex_connectivity() {
        let t = example_topology();

        for (id, vertices) in [vec![0], vec![1], vec![2], vec![3]].iter().enumerate() {
            let c = t.connectivity(ReferenceCellType::Point, id, 0).unwrap();
            for (i, j) in c.iter().zip(vertices) {
                assert_eq!(i.0, ReferenceCellType::Point);
                assert_eq!(i.1, *j);
            }
        }

        for (id, edges) in [vec![1, 2], vec![0, 2, 3], vec![0, 1, 4], vec![3, 4]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Point, id, 1).unwrap();
            for (i, j) in c.iter().zip(edges) {
                assert_eq!(i.0, ReferenceCellType::Interval);
                assert_eq!(i.1, *j);
            }
        }

        for (id, faces) in [vec![0], vec![0, 1], vec![0, 1], vec![1]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Point, id, 2).unwrap();
            for (i, j) in c.iter().zip(faces) {
                assert_eq!(i.0, ReferenceCellType::Triangle);
                assert_eq!(i.1, *j);
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
            let c = t.connectivity(ReferenceCellType::Interval, id, 0).unwrap();
            for (i, j) in c.iter().zip(vertices) {
                assert_eq!(i.0, ReferenceCellType::Point);
                assert_eq!(i.1, *j);
            }
        }

        for (id, edges) in [vec![0], vec![1], vec![2], vec![3], vec![4]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Interval, id, 1).unwrap();
            for (i, j) in c.iter().zip(edges) {
                assert_eq!(i.0, ReferenceCellType::Interval);
                assert_eq!(i.1, *j);
            }
        }

        for (id, faces) in [vec![0, 1], vec![0], vec![0], vec![1], vec![1]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Interval, id, 2).unwrap();
            for (i, j) in c.iter().zip(faces) {
                assert_eq!(i.0, ReferenceCellType::Triangle);
                assert_eq!(i.1, *j);
            }
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
            .cell_entities(ReferenceCellType::Quadrilateral, 0)
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

        let c = t.cell_entities(ReferenceCellType::Triangle, 0).unwrap();
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
            .cell_entities(ReferenceCellType::Quadrilateral, 1)
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

        let c = t.cell_entities(ReferenceCellType::Triangle, 1).unwrap();
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
            .cell_entities(ReferenceCellType::Quadrilateral, 2)
            .unwrap();
        assert_eq!(c.len(), 1);
        // cell 0
        assert_eq!(c[0].0, ReferenceCellType::Quadrilateral);
        assert_eq!(c[0].1, 0);

        let c = t.cell_entities(ReferenceCellType::Triangle, 2).unwrap();
        assert_eq!(c.len(), 1);
        // cell 1
        assert_eq!(c[0].0, ReferenceCellType::Triangle);
        assert_eq!(c[0].1, 0);
    }
    #[test]
    fn test_mixed_vertex_connectivity() {
        let t = example_topology_mixed();

        for (id, vertices) in [vec![0], vec![1], vec![2], vec![3], vec![4]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Point, id, 0).unwrap();
            for (i, j) in c.iter().zip(vertices) {
                assert_eq!(i.0, ReferenceCellType::Point);
                assert_eq!(i.1, *j);
            }
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
            for (i, j) in c.iter().zip(edges) {
                assert_eq!(i.0, ReferenceCellType::Interval);
                assert_eq!(i.1, *j);
            }
        }

        for (id, faces) in [
            vec![(ReferenceCellType::Quadrilateral, 0)],
            vec![
                (ReferenceCellType::Quadrilateral, 0),
                (ReferenceCellType::Triangle, 0),
            ],
            vec![(ReferenceCellType::Quadrilateral, 0)],
            vec![
                (ReferenceCellType::Quadrilateral, 0),
                (ReferenceCellType::Triangle, 0),
            ],
            vec![(ReferenceCellType::Triangle, 0)],
        ]
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
            for (i, j) in c.iter().zip(vertices) {
                assert_eq!(i.0, ReferenceCellType::Point);
                assert_eq!(i.1, *j);
            }
        }

        for (id, edges) in [vec![0], vec![1], vec![2], vec![3], vec![4], vec![5]]
            .iter()
            .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Interval, id, 1).unwrap();
            for (i, j) in c.iter().zip(edges) {
                assert_eq!(i.0, ReferenceCellType::Interval);
                assert_eq!(i.1, *j);
            }
        }

        for (id, faces) in [
            vec![(ReferenceCellType::Quadrilateral, 0)],
            vec![(ReferenceCellType::Quadrilateral, 0)],
            vec![
                (ReferenceCellType::Quadrilateral, 0),
                (ReferenceCellType::Triangle, 0),
            ],
            vec![(ReferenceCellType::Quadrilateral, 0)],
            vec![(ReferenceCellType::Triangle, 0)],
            vec![(ReferenceCellType::Triangle, 0)],
        ]
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
                let ce = t.cell_entities(*cell_type, dim).unwrap();
                let n = reference_cell::entity_counts(*cell_type)[dim];
                for i in 0..ce.len() / n {
                    // TODO: entity_count instead?
                    let con = t.connectivity(*cell_type, i, dim).unwrap();
                    assert_eq!(con, &ce[n * i..n * (i + 1)]);
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
    fn test_cell_entities_vs_connectivity_mixed() {
        let t = example_topology_mixed();
        cell_entities_vs_connectivity(&t);
    }
}
