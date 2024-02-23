//! Implementation of grid topology

use crate::grid::traits::{Ownership, Topology};
use crate::reference_cell::{dim, entity_counts, entity_types};
use crate::traits::cell::ReferenceCellType;
use bempp_element::cell;
use bempp_tools::arrays::AdjacencyList;
use bempp_traits::arrays::AdjacencyListAccess;
use bempp_traits::cell::ReferenceCell;

/// Topology of a serial grid
pub struct SerialTopology {
    dim: usize,
    connectivity: Vec<Vec<AdjacencyList<usize>>>,
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
    pub fn new(cells: &AdjacencyList<usize>, cell_types: &[ReferenceCellType]) -> Self {
        let mut index_map = vec![];
        let mut vertices = vec![];
        let mut starts = vec![];
        let mut cell_types_new = vec![];
        let dim = get_reference_cell(cell_types[0]).dim();

        let mut connectivity = vec![];
        for i in 0..dim + 1 {
            connectivity.push(vec![]);
            for _j in 0..dim + 1 {
                connectivity[i].push(AdjacencyList::<usize>::new());
            }
        }

        // dim0 = dim, dim1 = 0
        for c in cell_types {
            if dim != get_reference_cell(*c).dim() {
                panic!("Grids with cells of mixed topological dimension not supported.");
            }
            if !cell_types_new.contains(c) {
                starts.push(connectivity[dim][0].num_rows());
                cell_types_new.push(*c);
                let n = get_reference_cell(*c).vertex_count();
                for (i, cell) in cells.iter_rows().enumerate() {
                    if cell_types[i] == *c {
                        index_map.push(i);
                        // Note: this hard codes that the first n points are at the vertices
                        let mut row = vec![];
                        for v in &cell[..n] {
                            if !vertices.contains(v) {
                                vertices.push(*v);
                            }
                            row.push(vertices.iter().position(|&r| r == *v).unwrap());
                        }
                        connectivity[dim][0].add_row(&row);
                    }
                }
            }
        }

        // dim1 == 0
        for dim0 in 1..dim {
            let mut cty = AdjacencyList::<usize>::new();
            let cells = &connectivity[dim][0];
            for (i, cell_type) in cell_types_new.iter().enumerate() {
                let ref_cell = get_reference_cell(*cell_type);
                let ref_entities = (0..entity_counts(*cell_type)[dim0])
                    .map(|x| ref_cell.connectivity(dim0, x, 0).unwrap())
                    .collect::<Vec<Vec<usize>>>();

                let cstart = starts[i];
                let cend = if i == starts.len() - 1 {
                    connectivity[2][0].num_rows()
                } else {
                    starts[i + 1]
                };
                for c in cstart..cend {
                    let cell = cells.row(c).unwrap();
                    for e in &ref_entities {
                        let vertices = e.iter().map(|x| cell[*x]).collect::<Vec<usize>>();
                        let mut found = false;
                        for entity in cty.iter_rows() {
                            if all_equal(entity, &vertices) {
                                found = true;
                                break;
                            }
                        }
                        if !found {
                            cty.add_row(&vertices);
                        }
                    }
                }
            }
            connectivity[dim0][0] = cty;
        }

        // dim0 == dim1 == 0
        let mut nvertices = 0;
        let mut cty = AdjacencyList::<usize>::new();
        let cells = &connectivity[dim][0];
        for cell in cells.iter_rows() {
            for j in cell {
                if *j >= nvertices {
                    nvertices = *j + 1;
                }
            }
        }
        for i in 0..nvertices {
            cty.add_row(&[i]);
        }
        connectivity[0][0] = cty;

        // dim0 == dim1
        for (dim0, c) in connectivity.iter_mut().enumerate().skip(1) {
            for i in 0..c[0].num_rows() {
                c[dim0].add_row(&[i]);
            }
        }

        // dim0 == dim
        for dim1 in 1..dim + 1 {
            let mut cty = AdjacencyList::<usize>::new();
            let entities0 = &connectivity[dim][0];
            let entities1 = &connectivity[dim1][0];

            let mut sub_cell_types = vec![ReferenceCellType::Point; entities0.num_rows()];
            for (i, cell_type) in cell_types_new.iter().enumerate() {
                let etypes = &entity_types(*cell_type)[dim];

                let cstart = starts[i];
                let cend = if i == starts.len() - 1 {
                    connectivity[2][0].num_rows()
                } else {
                    starts[i + 1]
                };
                for t in sub_cell_types.iter_mut().skip(cstart).take(cend) {
                    *t = etypes[0];
                }
            }
            for (ei, entity0) in entities0.iter_rows().enumerate() {
                let entity = get_reference_cell(sub_cell_types[ei]);
                let mut row = vec![];
                for i in 0..entity.entity_count(dim1) {
                    let vertices = entity
                        .connectivity(dim1, i, 0)
                        .unwrap()
                        .iter()
                        .map(|x| entity0[*x])
                        .collect::<Vec<usize>>();
                    for (j, entity1) in entities1.iter_rows().enumerate() {
                        if all_equal(&vertices, entity1) {
                            row.push(j);
                            break;
                        }
                    }
                }
                cty.add_row(&row);
            }
            connectivity[dim][dim1] = cty
        }

        // dim1 < dim0
        for dim1 in 1..dim + 1 {
            for dim0 in dim1 + 1..dim {
                let mut cty = AdjacencyList::<usize>::new();
                let entities0 = &connectivity[dim0][0];
                let entities1 = &connectivity[dim1][0];
                let cell_to_entities0 = &connectivity[dim][dim0];

                let mut sub_cell_types = vec![ReferenceCellType::Point; entities0.num_rows()];
                for (i, cell_type) in cell_types_new.iter().enumerate() {
                    let etypes = &entity_types(*cell_type)[dim0];

                    let cstart = starts[i];
                    let cend = if i == starts.len() - 1 {
                        connectivity[2][0].num_rows()
                    } else {
                        starts[i + 1]
                    };
                    for c in cstart..cend {
                        for (e, t) in cell_to_entities0.row(c).unwrap().iter().zip(etypes) {
                            sub_cell_types[*e] = *t;
                        }
                    }
                }
                for (ei, entity0) in entities0.iter_rows().enumerate() {
                    let entity = get_reference_cell(sub_cell_types[ei]);
                    let mut row = vec![];
                    for i in 0..entity.entity_count(dim1) {
                        let vertices = entity
                            .connectivity(dim1, i, 0)
                            .unwrap()
                            .iter()
                            .map(|x| entity0[*x])
                            .collect::<Vec<usize>>();
                        for (j, entity1) in entities1.iter_rows().enumerate() {
                            if all_equal(&vertices, entity1) {
                                row.push(j);
                                break;
                            }
                        }
                    }
                    cty.add_row(&row);
                }
                connectivity[dim0][dim1] = cty;
            }
        }

        // dim1 > dim0
        for dim1 in 1..dim + 1 {
            for dim0 in 0..dim1 {
                let mut data = vec![vec![]; connectivity[dim0][0].num_rows()];
                for (i, row) in connectivity[dim1][dim0].iter_rows().enumerate() {
                    for v in row {
                        data[*v].push(i);
                    }
                }
                for row in data {
                    connectivity[dim0][dim1].add_row(&row);
                }
            }
        }

        let mut all_entity_types = vec![vec![], vec![], vec![], vec![]];
        for c in &cell_types_new {
            let et = entity_types(*c);
            for (dim, t) in et.iter().enumerate() {
                for e in t {
                    if !all_entity_types[dim].contains(e) {
                        all_entity_types[dim].push(*e);
                    }
                }
            }
        }

        Self {
            dim,
            connectivity,
            index_map,
            starts,
            cell_types: cell_types_new,
            entity_types: all_entity_types,
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
    //type Connectivity = AdjacencyList<usize>;

    fn dim(&self) -> usize {
        self.dim
    }
    fn index_map(&self) -> &[usize] {
        &self.index_map
    }
    fn entity_count(&self, dim: usize) -> usize {
        self.connectivity[dim][0].num_rows()
    }
    fn cell(&self, index: usize) -> Option<&[usize]> {
        if index < self.entity_count(self.dim) {
            Some(self.connectivity[self.dim][0].row(index).unwrap())
        } else {
            None
        }
    }
    fn cell_type(&self, index: usize) -> Option<ReferenceCellType> {
        for (i, start) in self.starts.iter().enumerate() {
            let end = if i == self.starts.len() - 1 {
                self.connectivity[2][0].num_rows()
            } else {
                self.starts[i + 1]
            };
            if *start <= index && index < end {
                return Some(self.cell_types[i]);
            }
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
        Some(&self.connectivity[self.dim][dim(etype)].data)
    }

    /// Get the indices of entities of dimension `dim` that are connected to the entity of type `etype` with index `index`
    fn connectivity(
        &self,
        _etype: ReferenceCellType,
        _index: usize,
        _dim: usize,
    ) -> Option<&[usize]> {
        None
    }

    fn entity_ownership(&self, _dim: usize, _index: usize) -> Ownership {
        Ownership::Owned
    }
}

#[cfg(test)]
mod test {
    use crate::grid::topology::*;

    fn example_topology() -> SerialTopology {
        SerialTopology::new(
            &AdjacencyList::from_data(vec![0, 1, 2, 2, 1, 3], vec![0, 3, 6]),
            &vec![ReferenceCellType::Triangle; 2],
        )
    }

    /*fn example_topology_mixed() -> SerialTopology {
        SerialTopology::new(
            &AdjacencyList::from_data(vec![0, 1, 3, 4, 1, 2, 4], vec![0, 4, 7]),
            &vec![ReferenceCellType::Quadrilateral, ReferenceCellType::Triangle],
        )
    }*/

    #[test]
    fn test_counts() {
        let t = example_topology();
        assert_eq!(t.dim(), 2);
        assert_eq!(t.entity_count(0), 4);
        assert_eq!(t.entity_count(1), 5);
        assert_eq!(t.entity_count(2), 2);
    }

    #[test]
    fn test_cell_connectivity_points() {
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
    fn test_cell_connectivity_intervals() {
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
    fn test_cell_connectivity_triangles() {
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
    /*
        fn test_cell_connectivity() {
            let t = example_topology();
            assert_eq!(t.connectivity(0, 0).row(0).unwrap().len(), 1);
            assert_eq!(t.connectivity(0, 0).row(0).unwrap()[0], 0);
            assert_eq!(t.connectivity(0, 0).row(1).unwrap().len(), 1);
            assert_eq!(t.connectivity(0, 0).row(1).unwrap()[0], 1);
            assert_eq!(t.connectivity(0, 0).row(2).unwrap().len(), 1);
            assert_eq!(t.connectivity(0, 0).row(2).unwrap()[0], 2);
            assert_eq!(t.connectivity(0, 0).row(3).unwrap().len(), 1);
            assert_eq!(t.connectivity(0, 0).row(3).unwrap()[0], 3);

            assert_eq!(t.connectivity(1, 0).row(0).unwrap().len(), 2);
            assert_eq!(t.connectivity(1, 0).row(0).unwrap()[0], 1);
            assert_eq!(t.connectivity(1, 0).row(0).unwrap()[1], 2);
            assert_eq!(t.connectivity(1, 0).row(1).unwrap().len(), 2);
            assert_eq!(t.connectivity(1, 0).row(1).unwrap()[0], 0);
            assert_eq!(t.connectivity(1, 0).row(1).unwrap()[1], 2);
            assert_eq!(t.connectivity(1, 0).row(2).unwrap().len(), 2);
            assert_eq!(t.connectivity(1, 0).row(2).unwrap()[0], 0);
            assert_eq!(t.connectivity(1, 0).row(2).unwrap()[1], 1);
            assert_eq!(t.connectivity(1, 0).row(3).unwrap().len(), 2);
            assert_eq!(t.connectivity(1, 0).row(3).unwrap()[0], 1);
            assert_eq!(t.connectivity(1, 0).row(3).unwrap()[1], 3);
            assert_eq!(t.connectivity(1, 0).row(4).unwrap().len(), 2);
            assert_eq!(t.connectivity(1, 0).row(4).unwrap()[0], 2);
            assert_eq!(t.connectivity(1, 0).row(4).unwrap()[1], 3);

            assert_eq!(t.connectivity(0, 1).row(0).unwrap().len(), 2);
            assert_eq!(t.connectivity(0, 1).row(0).unwrap()[0], 1);
            assert_eq!(t.connectivity(0, 1).row(0).unwrap()[1], 2);
            assert_eq!(t.connectivity(0, 1).row(1).unwrap().len(), 3);
            assert_eq!(t.connectivity(0, 1).row(1).unwrap()[0], 0);
            assert_eq!(t.connectivity(0, 1).row(1).unwrap()[1], 2);
            assert_eq!(t.connectivity(0, 1).row(1).unwrap()[2], 3);
            assert_eq!(t.connectivity(0, 1).row(2).unwrap().len(), 3);
            assert_eq!(t.connectivity(0, 1).row(2).unwrap()[0], 0);
            assert_eq!(t.connectivity(0, 1).row(2).unwrap()[1], 1);
            assert_eq!(t.connectivity(0, 1).row(2).unwrap()[2], 4);
            assert_eq!(t.connectivity(0, 1).row(3).unwrap().len(), 2);
            assert_eq!(t.connectivity(0, 1).row(3).unwrap()[0], 3);
            assert_eq!(t.connectivity(0, 1).row(3).unwrap()[1], 4);

            assert_eq!(t.connectivity(0, 2).row(0).unwrap().len(), 1);
            assert_eq!(t.connectivity(0, 2).row(0).unwrap()[0], 0);
            assert_eq!(t.connectivity(0, 2).row(1).unwrap().len(), 2);
            assert_eq!(t.connectivity(0, 2).row(1).unwrap()[0], 0);
            assert_eq!(t.connectivity(0, 2).row(1).unwrap()[1], 1);
            assert_eq!(t.connectivity(0, 2).row(2).unwrap().len(), 2);
            assert_eq!(t.connectivity(0, 2).row(2).unwrap()[0], 0);
            assert_eq!(t.connectivity(0, 2).row(2).unwrap()[1], 1);
            assert_eq!(t.connectivity(0, 2).row(3).unwrap().len(), 1);
            assert_eq!(t.connectivity(0, 2).row(3).unwrap()[0], 1);

            assert_eq!(t.connectivity(1, 2).row(0).unwrap().len(), 2);
            assert_eq!(t.connectivity(1, 2).row(0).unwrap()[0], 0);
            assert_eq!(t.connectivity(1, 2).row(0).unwrap()[1], 1);
            assert_eq!(t.connectivity(1, 2).row(1).unwrap().len(), 1);
            assert_eq!(t.connectivity(1, 2).row(1).unwrap()[0], 0);
            assert_eq!(t.connectivity(1, 2).row(2).unwrap().len(), 1);
            assert_eq!(t.connectivity(1, 2).row(2).unwrap()[0], 0);
            assert_eq!(t.connectivity(1, 2).row(3).unwrap().len(), 1);
            assert_eq!(t.connectivity(1, 2).row(3).unwrap()[0], 1);
            assert_eq!(t.connectivity(1, 2).row(4).unwrap().len(), 1);
            assert_eq!(t.connectivity(1, 2).row(4).unwrap()[0], 1);
        }
    */
}
