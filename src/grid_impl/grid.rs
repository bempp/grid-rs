use crate::grid_impl::traits::{Geometry, GeometryEvaluator, Grid, Topology};
use crate::reference_cell;
use crate::reference_cell::ReferenceCellType;
use crate::traits::{
    cell::CellType, geometry::GeometryType, grid::GridType, point::PointType,
    reference_map::ReferenceMapType, topology::TopologyType,
};
use crate::types::CellLocalIndexPair;
use bempp_traits::element::FiniteElement;
use num::Float;
use rlst_common::types::Scalar;
use std::iter::Copied;

pub struct Point<'a, T: Float, G: Geometry> {
    geometry: &'a G,
    index: usize,
    _t: std::marker::PhantomData<T>,
}
pub struct Cell<'a, T: Float, GridImpl: Grid> {
    grid: &'a GridImpl,
    index: usize,
    _t: std::marker::PhantomData<T>,
}
pub struct CellTopology<'a, GridImpl: Grid> {
    topology: &'a <GridImpl as Grid>::Topology,
    index: <<GridImpl as Grid>::Topology as Topology>::IndexType,
}
pub struct CellGeometry<'a, T: Float, GridImpl: Grid> {
    geometry: &'a <GridImpl as Grid>::Geometry,
    index: <<GridImpl as Grid>::Geometry as Geometry>::IndexType,
    _t: std::marker::PhantomData<T>,
}
pub struct ReferenceMap<'a, GridImpl: Grid> {
    grid: &'a GridImpl,
    evaluator: <<GridImpl as Grid>::Geometry as Geometry>::Evaluator<'a>,
}
pub struct PointIterator<'a, GridImpl: Grid, Iter: std::iter::Iterator<Item = usize>> {
    iter: Iter,
    geometry: &'a <GridImpl as Grid>::Geometry,
}

impl<'a, GridImpl: Grid, Iter: std::iter::Iterator<Item = usize>> std::iter::Iterator
    for PointIterator<'a, GridImpl, Iter>
{
    type Item = Point<'a, <GridImpl as Grid>::T, <GridImpl as Grid>::Geometry>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.iter.next() {
            Some(Point {
                geometry: self.geometry,
                index,
                _t: std::marker::PhantomData,
            })
        } else {
            None
        }
    }
}

impl<'a, T: Float, G: Geometry<T = T>> PointType for Point<'a, T, G> {
    type T = T;
    fn coords(&self, data: &mut [Self::T]) {
        assert_eq!(data.len(), self.geometry.dim());
        for (dim, d) in data.iter_mut().enumerate() {
            *d = *self.geometry.coordinate(self.index, dim).unwrap();
        }
    }
    fn index(&self) -> usize {
        self.index
    }
    fn id(&self) -> usize {
        self.index
    }
}

impl<'grid, T: Float + Scalar, GridImpl: Grid<T = T>> CellType for Cell<'grid, T, GridImpl>
where
    GridImpl: 'grid,
{
    type Grid = GridImpl;

    type Topology<'a> = CellTopology<'a, GridImpl> where Self: 'a;
    type Geometry<'a> = CellGeometry<'a, T, GridImpl> where Self: 'a;

    fn id(&self) -> usize {
        self.index
    }
    fn index(&self) -> usize {
        self.index
    }

    fn topology(&self) -> Self::Topology<'_> {
        CellTopology::<'_, GridImpl> {
            topology: self.grid.topology(),
            index: self.grid.topology().index_map()[self.index],
        }
    }

    fn grid(&self) -> &Self::Grid {
        self.grid
    }

    fn geometry(&self) -> Self::Geometry<'_> {
        CellGeometry::<'_, T, GridImpl> {
            geometry: self.grid.geometry(),
            index: self.grid.geometry().index_map()[self.index],
            _t: std::marker::PhantomData,
        }
    }
}

impl<'grid, T: Float + Scalar, GridImpl: Grid<T = T>> TopologyType for CellTopology<'grid, GridImpl>
where
    GridImpl: 'grid,
{
    type Grid = GridImpl;
    type IndexType = <<GridImpl as Grid>::Topology as Topology>::IndexType;
    type VertexIndexIter<'a> = Copied<std::slice::Iter<'a, Self::IndexType>>
    where
        Self: 'a;
    type EdgeIndexIter<'a> = Self::VertexIndexIter<'a>
    where
        Self: 'a;
    type FaceIndexIter<'a> = Self::VertexIndexIter<'a>
    where
        Self: 'a;

    fn vertex_indices(&self) -> Self::VertexIndexIter<'_> {
        self.topology
            .cell_to_entities(self.index, 0)
            .unwrap()
            .iter()
            .copied()
    }

    fn edge_indices(&self) -> Self::EdgeIndexIter<'_> {
        self.topology
            .cell_to_entities(self.index, 1)
            .unwrap()
            .iter()
            .copied()
    }

    fn face_indices(&self) -> Self::FaceIndexIter<'_> {
        self.topology
            .cell_to_entities(self.index, 2)
            .unwrap()
            .iter()
            .copied()
    }

    fn cell_type(&self) -> ReferenceCellType {
        self.topology.cell_type(self.index).unwrap()
    }
}

impl<'grid, T: Float + Scalar, GridImpl: Grid<T = T>> GeometryType
    for CellGeometry<'grid, T, GridImpl>
where
    GridImpl: 'grid,
{
    type Grid = GridImpl;

    type VertexIterator<'iter> =
        PointIterator<'iter, Self::Grid, Copied<std::slice::Iter<'iter, usize>>> where Self: 'iter;

    type PointIterator<'iter> = Self::VertexIterator<'iter> where Self: 'iter;

    fn physical_dimension(&self) -> usize {
        self.geometry.dim()
    }

    fn midpoint(&self, point: &mut [T]) {
        self.geometry.midpoint(self.index, point)
    }

    fn diameter(&self) -> T {
        self.geometry.diameter(self.index)
    }

    fn volume(&self) -> T {
        self.geometry.volume(self.index)
    }

    fn points(&self) -> Self::PointIterator<'_> {
        PointIterator {
            iter: self
                .geometry
                .cell_points(self.index)
                .unwrap()
                .iter()
                .copied(),
            geometry: self.geometry,
        }
    }
    fn vertices(&self) -> Self::VertexIterator<'_> {
        let cell_type = self.geometry.cell_element(self.index).unwrap().cell_type();
        let nvertices = reference_cell::entity_counts(cell_type)[0];
        PointIterator {
            iter: self.geometry.cell_points(self.index).unwrap()[..nvertices]
                .iter()
                .copied(),
            geometry: self.geometry,
        }
    }
}

impl<'a, T: Float + Scalar, GridImpl: Grid<T = T>> ReferenceMapType for ReferenceMap<'a, GridImpl> {
    type Grid = GridImpl;

    fn domain_dimension(&self) -> usize {
        self.grid.topology().dim()
    }

    fn physical_dimension(&self) -> usize {
        self.grid.geometry().dim()
    }

    fn number_of_reference_points(&self) -> usize {
        self.evaluator.point_count()
    }

    fn reference_to_physical(
        &self,
        cell_index: usize,
        point_index: usize,
        value: &mut [<Self::Grid as GridType>::T],
    ) {
        self.evaluator.compute_point(cell_index, point_index, value);
    }

    fn jacobian(&self, cell_index: usize, point_index: usize, value: &mut [T]) {
        self.evaluator
            .compute_jacobian(cell_index, point_index, value);
    }

    fn normal(&self, cell_index: usize, point_index: usize, value: &mut [T]) {
        self.evaluator
            .compute_normal(cell_index, point_index, value);
    }
}

impl<'grid, T: Float + Scalar, GridImpl: Grid<T = T>> GridType for GridImpl
where
    GridImpl: 'grid,
{
    type T = T;
    type IndexType = <<GridImpl as Grid>::Topology as Topology>::IndexType;

    type ReferenceMap<'a> = ReferenceMap<'a, GridImpl>
    where
        Self: 'a;

    type Point<'a> = Point<'a, T, GridImpl::Geometry> where Self: 'a;
    type Cell<'a> = Cell<'a, T, GridImpl> where Self: 'a;

    fn number_of_points(&self) -> usize {
        self.geometry().point_count()
    }
    fn number_of_vertices(&self) -> usize {
        self.topology().entity_count(ReferenceCellType::Point)
    }
    fn number_of_cells(&self) -> usize {
        self.geometry().cell_count()
    }

    fn point_index_from_id(&self, id: usize) -> usize {
        id
    }
    fn point_id_from_index(&self, index: usize) -> usize {
        index
    }

    fn cell_index_from_id(&self, id: usize) -> usize {
        id
    }
    fn cell_id_from_index(&self, index: usize) -> usize {
        index
    }

    fn point_from_index(&self, index: usize) -> Self::Point<'_> {
        Self::Point {
            geometry: self.geometry(),
            index,
            _t: std::marker::PhantomData,
        }
    }

    fn cell_from_index(&self, index: usize) -> Self::Cell<'_> {
        Self::Cell {
            grid: self,
            index,
            _t: std::marker::PhantomData,
        }
    }

    fn reference_to_physical_map<'a>(
        &'a self,
        reference_points: &'a [Self::T],
    ) -> Self::ReferenceMap<'a> {
        Self::ReferenceMap {
            grid: self,
            evaluator: self.geometry().get_evaluator(reference_points),
        }
    }

    fn vertex_to_cells(
        &self,
        vertex_index: Self::IndexType,
    ) -> &[CellLocalIndexPair<Self::IndexType>] {
        self.topology().entity_to_cells(0, vertex_index).unwrap()
    }

    fn edge_to_cells(&self, edge_index: Self::IndexType) -> &[CellLocalIndexPair<Self::IndexType>] {
        self.topology().entity_to_cells(1, edge_index).unwrap()
    }

    fn face_to_cells(&self, face_index: Self::IndexType) -> &[CellLocalIndexPair<Self::IndexType>] {
        self.topology().entity_to_cells(2, face_index).unwrap()
    }
}

#[cfg(test)]
mod test {
    use crate::grid_impl::grid::*;
    use crate::grid_impl::mixed_grid::SerialMixedGrid;
    use crate::grid_impl::single_element_grid::SerialSingleElementGrid;
    use rlst_dense::{rlst_dynamic_array2, traits::RandomAccessMut};

    #[test]
    fn test_grid_mixed_cell_type() {
        let mut points = rlst_dynamic_array2!(f64, [9, 3]);
        *points.get_mut([0, 0]).unwrap() = -1.0;
        *points.get_mut([0, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = -0.5;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([1, 2]).unwrap() = 0.2;
        *points.get_mut([2, 0]).unwrap() = 0.0;
        *points.get_mut([2, 1]).unwrap() = 0.0;
        *points.get_mut([2, 2]).unwrap() = 0.0;
        *points.get_mut([3, 0]).unwrap() = 1.0;
        *points.get_mut([3, 1]).unwrap() = 0.0;
        *points.get_mut([3, 2]).unwrap() = 0.0;
        *points.get_mut([4, 0]).unwrap() = 2.0;
        *points.get_mut([4, 1]).unwrap() = 0.0;
        *points.get_mut([4, 2]).unwrap() = 0.0;
        *points.get_mut([5, 0]).unwrap() = -std::f64::consts::FRAC_1_SQRT_2;
        *points.get_mut([5, 1]).unwrap() = std::f64::consts::FRAC_1_SQRT_2;
        *points.get_mut([5, 2]).unwrap() = 0.0;
        *points.get_mut([6, 0]).unwrap() = 0.0;
        *points.get_mut([6, 1]).unwrap() = 0.5;
        *points.get_mut([6, 2]).unwrap() = 0.0;
        *points.get_mut([7, 0]).unwrap() = 0.0;
        *points.get_mut([7, 1]).unwrap() = 1.0;
        *points.get_mut([7, 2]).unwrap() = 0.0;
        *points.get_mut([8, 0]).unwrap() = 1.0;
        *points.get_mut([8, 1]).unwrap() = 1.0;
        *points.get_mut([8, 2]).unwrap() = 0.0;
        let grid = SerialMixedGrid::<f64>::new(
            points,
            &[0, 2, 7, 6, 5, 1, 2, 3, 7, 8, 3, 4, 8],
            &[
                ReferenceCellType::Triangle,
                ReferenceCellType::Quadrilateral,
                ReferenceCellType::Triangle,
            ],
            &[2, 1, 1],
        );

        assert_eq!(grid.number_of_vertices(), 6);
        assert_eq!(grid.number_of_points(), 9);
        assert_eq!(grid.number_of_cells(), 3);

        let mut coords = vec![0.0; grid.geometry().dim()];
        for point in grid.iter_all_points() {
            point.coords(coords.as_mut_slice());
            println!("{:#?}", coords);
        }

        for cell in grid.iter_all_cells() {
            println!("{:?}", cell.index());
        }
        for cell in grid.iter_all_cells() {
            for (local_index, (vertex_index, edge_index)) in cell
                .topology()
                .vertex_indices()
                .zip(cell.topology().edge_indices())
                .enumerate()
            {
                println!(
                    "Cell: {}, Vertex: {}, {:?}, Edge: {}, {:?}, Volume: {}",
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

    #[test]
    fn test_grid_single_element() {
        let mut points = rlst_dynamic_array2!(f64, [9, 3]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.5;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([1, 2]).unwrap() = 0.2;
        *points.get_mut([2, 0]).unwrap() = 1.0;
        *points.get_mut([2, 1]).unwrap() = 0.0;
        *points.get_mut([2, 2]).unwrap() = 0.0;
        *points.get_mut([3, 0]).unwrap() = 0.0;
        *points.get_mut([3, 1]).unwrap() = 0.5;
        *points.get_mut([3, 2]).unwrap() = 0.0;
        *points.get_mut([4, 0]).unwrap() = 0.5;
        *points.get_mut([4, 1]).unwrap() = 0.5;
        *points.get_mut([4, 2]).unwrap() = 0.0;
        *points.get_mut([5, 0]).unwrap() = 1.0;
        *points.get_mut([5, 1]).unwrap() = 0.5;
        *points.get_mut([5, 2]).unwrap() = 0.0;
        *points.get_mut([6, 0]).unwrap() = 0.0;
        *points.get_mut([6, 1]).unwrap() = 1.0;
        *points.get_mut([6, 2]).unwrap() = 0.0;
        *points.get_mut([7, 0]).unwrap() = 0.5;
        *points.get_mut([7, 1]).unwrap() = 1.0;
        *points.get_mut([7, 2]).unwrap() = 0.0;
        *points.get_mut([8, 0]).unwrap() = 1.0;
        *points.get_mut([8, 1]).unwrap() = 1.0;
        *points.get_mut([8, 2]).unwrap() = 0.0;
        let grid = SerialSingleElementGrid::<f64>::new(
            points,
            &[0, 2, 6, 4, 3, 1, 2, 8, 6, 7, 4, 5],
            ReferenceCellType::Triangle,
            2,
        );

        assert_eq!(grid.number_of_vertices(), 4);
        assert_eq!(grid.number_of_points(), 9);
        assert_eq!(grid.number_of_cells(), 2);

        let mut coords = vec![0.0; grid.geometry().dim()];
        for point in grid.iter_all_points() {
            point.coords(coords.as_mut_slice());
            println!("{:#?}", coords);
        }

        for cell in grid.iter_all_cells() {
            println!("{:?}", cell.index());
        }
        for cell in grid.iter_all_cells() {
            for (local_index, (vertex_index, edge_index)) in cell
                .topology()
                .vertex_indices()
                .zip(cell.topology().edge_indices())
                .enumerate()
            {
                println!(
                    "Cell: {}, Vertex: {}, {:?}, Edge: {}, {:?}, Volume: {}",
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
