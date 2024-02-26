//! Implementation of grid topology

use crate::grid::traits::{Ownership, Topology};
use crate::reference_cell;
use crate::traits::cell::ReferenceCellType;
use bempp_element::cell;
use bempp_traits::cell::ReferenceCell;

use std::collections::HashMap;

/// Topology of a serial grid
pub struct SerialTopology {
    dim: usize,
    connectivity: Vec<Vec<Vec<Vec<usize>>>>,
    cell_connectivity: Vec<Vec<usize>>, // TODO: get rid of this input
    connectivity_neww: HashMap<ReferenceCellType, Vec<Vec<Vec<usize>>>>,
    index_map: Vec<usize>,
    starts: Vec<usize>,
    cell_types: Vec<ReferenceCellType>,
    entity_types: Vec<Vec<ReferenceCellType>>,
}

fn get_reference_cell(cell_type: ReferenceCellType) -> Box<dyn ReferenceCell> {
    match cell_type {
        ReferenceCellType::Interval => Box::new(cell::Interval),
        ReferenceCellType::Triangle => Box::new(cell::Triangle),
        ReferenceCellType::Quadrilateral => Box::new(cell::Quadrilateral),
        _ => {
            panic!("Unsupported cell type (for now)");
        }
    }
}

unsafe impl Sync for SerialTopology {}

impl SerialTopology {
    pub fn new(cells: &[usize], cell_types: &[ReferenceCellType]) -> Self {
        let mut index_map = vec![];
        let mut vertices = vec![];
        let mut starts = vec![];
        let mut cell_types_new = vec![];
        let dim = reference_cell::dim(cell_types[0]);

        let mut connectivity = vec![];
        for i in 0..dim + 1 {
            connectivity.push(vec![]);
            for _j in 0..dim + 1 {
                connectivity[i].push(vec![]);
            }
        }

        let mut entity_types = vec![vec![]; 4];
        let mut connectivity_neww = HashMap::new();
        for c in cell_types {
            for (dim0, etypes) in reference_cell::entity_types(*c).iter().enumerate() {
                for e in etypes {
                    if !entity_types[dim0].contains(e) {
                        entity_types[dim0].push(*e);
                        connectivity_neww.insert(*e, vec![vec![]; dim + 1]);
                    }
                }
            }
        }

        // dim0 = dim, dim1 = 0
        for c in cell_types {
            if dim != reference_cell::dim(*c) {
                panic!("Grids with cells of mixed topological dimension not supported.");
            }
            if !cell_types_new.contains(c) {
                let mut cty = vec![];
                starts.push(connectivity[dim][0].len());
                cell_types_new.push(*c);
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
                        connectivity[dim][0].push(row.clone());
                        cty.push(row);
                    }
                    start += reference_cell::entity_counts(*ct)[0];
                }
                connectivity_neww.get_mut(c).unwrap()[0] = cty;
            }
        }

        // dim1 == 0
        for dim0 in 1..dim {
            let mut cty: Vec<Vec<usize>> = vec![];
            for etype in &entity_types[dim0] {
                let mut cty_neww = vec![];
                let cells = &connectivity[dim][0];
                for (i, cell_type) in cell_types_new.iter().enumerate() {
                    let ref_cell = get_reference_cell(*cell_type);
                    let ref_entities = (0..reference_cell::entity_counts(*cell_type)[dim0])
                        .map(|x| ref_cell.connectivity(dim0, x, 0).unwrap())
                        .collect::<Vec<Vec<usize>>>();

                    let cstart = starts[i];
                    let cend = if i == starts.len() - 1 {
                        connectivity[2][0].len()
                    } else {
                        starts[i + 1]
                    };
                    for cell in cells.iter().take(cend).skip(cstart) {
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
                                cty.push(vertices.clone());
                                cty_neww.push(vertices);
                            }
                        }
                    }
                }
                connectivity_neww.get_mut(etype).unwrap()[0] = cty_neww;
            }
            connectivity[dim0][0] = cty;
        }

        // dim0 == dim1 == 0
        let mut nvertices = 0;
        let mut cty = vec![];
        let cells = &connectivity[dim][0];
        for cell in cells {
            for j in cell {
                if *j >= nvertices {
                    nvertices = *j + 1;
                }
            }
        }
        for i in 0..nvertices {
            cty.push(vec![i]);
        }
        connectivity[0][0] = cty.clone();
        connectivity_neww.get_mut(&ReferenceCellType::Point).unwrap()[0] = cty;

        // dim0 == dim1
        for (dim0, c) in connectivity.iter_mut().enumerate().skip(1) {
            for i in 0..c[0].len() {
                c[dim0].push(vec![i]);
            }
        }
        for (dim0, etypes) in entity_types.iter().enumerate().skip(1) {
            for etype in etypes {
                let mut cty = vec![];
                for i in 0..connectivity_neww[etype][0].len() {
                    cty.push(vec![i]);
                }
                connectivity_neww.get_mut(etype).unwrap()[dim0] = cty;
            }
        }

        // dim0 == dim
        for dim1 in 1..dim {
            let mut cty = vec![];
            let entities0 = &connectivity[dim][0];
            let entities1 = &connectivity[dim1][0];

            let mut sub_cell_types = vec![ReferenceCellType::Point; entities0.len()];
            for (i, cell_type) in cell_types_new.iter().enumerate() {
                let etypes = &reference_cell::entity_types(*cell_type)[dim];

                let cstart = starts[i];
                let cend = if i == starts.len() - 1 {
                    connectivity[2][0].len()
                } else {
                    starts[i + 1]
                };
                for t in sub_cell_types.iter_mut().skip(cstart).take(cend) {
                    *t = etypes[0];
                }
            }
            for (ei, entity0) in entities0.iter().enumerate() {
                let entity = get_reference_cell(sub_cell_types[ei]);
                let mut row = vec![];
                for i in 0..entity.entity_count(dim1) {
                    let vertices = entity
                        .connectivity(dim1, i, 0)
                        .unwrap()
                        .iter()
                        .map(|x| entity0[*x])
                        .collect::<Vec<usize>>();
                    for (j, entity1) in entities1.iter().enumerate() {
                        if all_equal(&vertices, entity1) {
                            row.push(j);
                            break;
                        }
                    }
                }
                cty.push(row);
            }
            connectivity[dim][dim1] = cty;
        }
        for cell_type in &entity_types[dim] {
            for dim1 in 1..dim {
                let mut cty = vec![];
                let entities0 = &connectivity_neww[cell_type][0];
                let mut start = 0;
                for etype in &entity_types[dim] {
                    let entities1 = &connectivity_neww[etype][0];

                    let mut sub_cell_types = vec![ReferenceCellType::Point; entities0.len()];
                    for (i, cell_type) in cell_types_new.iter().enumerate() {
                        let etypes = &reference_cell::entity_types(*cell_type)[dim];

                        let cstart = starts[i];
                        let cend = if i == starts.len() - 1 {
                            connectivity[2][0].len()
                        } else {
                            starts[i + 1]
                        };
                        for t in sub_cell_types.iter_mut().skip(cstart).take(cend) {
                            *t = etypes[0];
                        }
                    }
                    for (ei, entity0) in entities0.iter().enumerate() {
                        let entity = get_reference_cell(sub_cell_types[ei]);
                        let mut row = vec![];
                        for i in 0..entity.entity_count(dim1) {
                            let vertices = entity
                                .connectivity(dim1, i, 0)
                                .unwrap()
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
                connectivity_neww.get_mut(cell_type).unwrap()[dim1] = cty;
            }
        }

        // dim1 < dim0
        // TODO: continue from here
        for dim1 in 1..dim + 1 {
            for dim0 in dim1 + 1..dim {
                let mut cty = vec![];
                let entities0 = &connectivity[dim0][0];
                let entities1 = &connectivity[dim1][0];
                let cell_to_entities0 = &connectivity[dim][dim0]; // TODO

                let mut sub_cell_types = vec![ReferenceCellType::Point; entities0.len()];
                for (i, cell_type) in cell_types_new.iter().enumerate() {
                    let etypes = &reference_cell::entity_types(*cell_type)[dim0];

                    let cstart = starts[i];
                    let cend = if i == starts.len() - 1 {
                        connectivity[2][0].len()
                    } else {
                        starts[i + 1]
                    };
                    for ces in cell_to_entities0.iter().take(cend).skip(cstart) {
                        for (e, t) in ces.iter().zip(etypes) {
                            sub_cell_types[*e] = *t;
                        }
                    }
                }
                for (ei, entity0) in entities0.iter().enumerate() {
                    let entity = get_reference_cell(sub_cell_types[ei]);
                    let mut row = vec![];
                    for i in 0..entity.entity_count(dim1) {
                        let vertices = entity
                            .connectivity(dim1, i, 0)
                            .unwrap()
                            .iter()
                            .map(|x| entity0[*x])
                            .collect::<Vec<usize>>();
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
        for dim0 in 2..dim + 1 {
            println!("{dim0}");
            for etype0 in &entity_types[dim0] {
                for dim1 in 1..dim0 {
                    let mut cty = vec![];
                    let entities0 = &connectivity_neww[etype0][0];
                    let mut start = 0;
                    for etype1 in &entity_types[dim1] {
                        let entities1 = &connectivity_neww[etype1][0];
                        let mut sub_cell_types = vec![ReferenceCellType::Point; entities0.len()];
                        for cell_type in &entity_types[dim] {
                            let cell_to_entities0 = &connectivity_neww[cell_type][dim0];
                            let etypes = &reference_cell::entity_types(*cell_type)[dim0];

                            for ces in cell_to_entities0 {
                                for (e, t) in ces.iter().zip(etypes) {
                                    sub_cell_types[*e] = *t;
                                }
                            }
                        }
                        for (ei, entity0) in entities0.iter().enumerate() {
                            let entity = get_reference_cell(sub_cell_types[ei]);
                            let mut row = vec![];
                            for i in 0..entity.entity_count(dim1) {
                                let vertices = entity
                                    .connectivity(dim1, i, 0)
                                    .unwrap()
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
                    println!("{dim0} {dim1}");
                    connectivity_neww.get_mut(etype0).unwrap()[dim1] = cty;
                }
            }
        }

        // dim1 > dim0
        for dim1 in 1..dim + 1 {
            for dim0 in 0..dim1 {
                let mut data = vec![vec![]; connectivity[dim0][0].len()];
                for (i, row) in connectivity[dim1][dim0].iter().enumerate() {
                    for v in row {
                        data[*v].push(i);
                    }
                }
                for row in data {
                    connectivity[dim0][dim1].push(row);
                }
            }
        }
        for dim0 in 0..dim + 1 {
            for etype0 in &entity_types[dim0] {
                for dim1 in dim0 + 1..dim + 1 {
                    let mut data = vec![vec![]; connectivity_neww[etype0][0].len()];
                    let mut start = 0;
                    for etype1 in &entity_types[dim1] {
                        for (i, row) in connectivity_neww[etype1][dim0].iter().enumerate() {
                            for v in row {
                                data[*v].push(start + i);
                            }
                        }
                        start += connectivity_neww[etype1][dim0].len();
                    }
                    for row in data {
                        connectivity_neww.get_mut(etype0).unwrap()[dim1].push(row);
                    }
                }
            }
        }


        // Check that old and new connectivity match
        for (dim0, etypes) in entity_types.iter().enumerate().take(dim + 1) {
            for dim1 in 0..dim + 1 {
                println!("{dim0} {dim1}");
                println!("{:?}", connectivity[dim0][dim1]);
                println!("{:?}", connectivity_neww[&etypes[0]][dim1]);
                println!();
            }
        }
        for (dim0, etypes) in entity_types.iter().enumerate().take(dim + 1) {
            for dim1 in 0..dim + 1 {
                assert_eq!(connectivity[dim0][dim1], connectivity_neww[&etypes[0]][dim1]);
            }
        }

        let mut cell_connectivity = vec![];
        for i in 0..dim + 1 {
            let mut all = vec![];
            for c in &connectivity[dim][i] {
                all.extend_from_slice(c)
            }
            cell_connectivity.push(all);
        }

        Self {
            dim,
            connectivity,
            cell_connectivity,
            connectivity_neww,
            index_map,
            starts,
            cell_types: cell_types_new,
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
    fn entity_count(&self, dim: usize) -> usize {
        let mut count = 0;
        for etype in &self.entity_types[dim] {
            count += self.connectivity_neww[etype][0].len();
        }
        count
    }
    fn cell(&self, index: usize) -> Option<&[usize]> {
        let mut start = 0;
        for etype in &self.entity_types[self.dim] {
            let count = self.connectivity_neww[etype][0].len();
            if index < start + count {
                return Some(&self.connectivity_neww[etype][0][index - start])
            }
            start += count;
        }
        None
    }
    fn cell_type(&self, index: usize) -> Option<ReferenceCellType> {
        let mut start = 0;
        for etype in &self.entity_types[self.dim] {
            let count = self.connectivity_neww[etype][0].len();
            if index < start + count {
                return Some(*etype)
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
        _cell_type: ReferenceCellType,
        etype: ReferenceCellType,
    ) -> Option<&[usize]> {
        // TODO
        Some(&self.cell_connectivity[reference_cell::dim(etype)])
    }

    /// Get the indices of entities of dimension `dim` that are connected to the entity of type `etype` with index `index`
    fn connectivity(&self, etype: ReferenceCellType, index: usize, dim: usize) -> Option<&[usize]> {
        Some(&self.connectivity_neww[&etype][dim][index])
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
        assert_eq!(t.entity_count(0), 4);
        assert_eq!(t.entity_count(1), 5);
        assert_eq!(t.entity_count(2), 2);
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

    /*
    fn example_topology_mixed() -> SerialTopology {
        SerialTopology::new(
            &[0, 1, 3, 4, 1, 2, 4],
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
        assert_eq!(t.entity_count(0), 5);
        assert_eq!(t.entity_count(1), 6);
        assert_eq!(t.entity_count(2), 2);
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
        assert_eq!(c[2], 3);
        assert_eq!(c[3], 4);

        let c = t
            .cell_entities(ReferenceCellType::Triangle, ReferenceCellType::Point)
            .unwrap();
        assert_eq!(c.len(), 3);
        // cell 1
        assert_eq!(c[0], 1);
        assert_eq!(c[1], 2);
        assert_eq!(c[2], 4);
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
        assert_eq!(c[3], 4);
        assert_eq!(c[4], 2);
        assert_eq!(c[5], 5);
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
        assert_eq!(c[0], 0);
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
            vec![4, 5],
            vec![1, 3],
            vec![2, 3, 4],
        ]
        .iter()
        .enumerate()
        {
            let c = t.connectivity(ReferenceCellType::Point, id, 1).unwrap();
            assert_eq!(c, edges);
        }

        for (id, faces) in [vec![0], vec![0, 1], vec![1], vec![0], vec![0, 1]]
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
            vec![0, 3],
            vec![1, 4],
            vec![3, 4],
            vec![2, 4],
            vec![1, 2],
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
    */
}
