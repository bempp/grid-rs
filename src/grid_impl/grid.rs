use crate::grid_impl::traits::{Geometry, Grid, Topology};
use crate::reference_cell::ReferenceCellType;
use crate::traits::{
    cell::CellType, geometry::GeometryType, grid::GridType, point::PointType,
    topology::TopologyType,
};
use crate::types::vertex_iterator::PointIterator;
use num::Float;
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
    grid: &'a GridImpl,
    index: <<GridImpl as Grid>::Topology as Topology>::IndexType,
}
pub struct CellGeometry<'a, T: Float, GridImpl: Grid> {
    grid: &'a GridImpl,
    index: <<GridImpl as Grid>::Geometry as Geometry>::IndexType,
    _t: std::marker::PhantomData<T>,
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

impl<'grid, T: Float, GridImpl: Grid<T = T>> CellType for Cell<'grid, T, GridImpl>
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
            grid: self.grid, // TODO: replace with just topology
            index: self.grid.topology().index_map()[self.index],
        }
    }

    fn grid(&self) -> &Self::Grid {
        self.grid
    }

    fn geometry(&self) -> Self::Geometry<'_> {
        CellGeometry::<'_, T, GridImpl> {
            grid: self.grid,
            index: self.grid.geometry().index_map()[self.index],
            _t: std::marker::PhantomData,
        }
    }
}

impl<'grid, GridImpl: Grid> TopologyType for CellTopology<'grid, GridImpl>
where
    GridImpl: 'grid,
{
    type Grid = GridImpl;
    type IndexType = <<GridImpl as Grid>::Topology as Topology>::IndexType;
    type VertexIndexIter<'a> = Copied<std::slice::Iter<'a, Self::IndexType>>
    where
        Self: 'a;

    type EdgeIndexIter<'a> = Copied<std::slice::Iter<'a, Self::IndexType>>
    where
        Self: 'a;

    fn vertex_indices(&self) -> Self::VertexIndexIter<'_> {
        self.grid
            .topology()
            .connectivity(self.index, 0)
            .unwrap()
            .iter()
            .copied()
    }

    fn edge_indices(&self) -> Self::EdgeIndexIter<'_> {
        self.grid
            .topology()
            .connectivity(self.index, 1)
            .unwrap()
            .iter()
            .copied()
    }
}

impl<'grid, T: Float, GridImpl: Grid<T = T>> GeometryType for CellGeometry<'grid, T, GridImpl>
where
    GridImpl: 'grid,
{
    type Grid = GridImpl;

    type VertexIterator<'iter> =
        PointIterator<'iter, Self::Grid, Copied<std::slice::Iter<'iter, usize>>> where Self: 'iter;

    type PointIterator<'iter> = Self::VertexIterator<'iter> where Self: 'iter;

    fn physical_dimension(&self) -> usize {
        self.grid.geometry().dim()
    }

    fn midpoint(&self, point: &mut [T]) {
        self.grid.geometry().midpoint(self.index, point)
    }

    fn diameter(&self) -> T {
        self.grid.geometry().diameter(self.index)
    }

    fn volume(&self) -> T {
        self.grid.geometry().volume(self.index)
    }

    fn points(&self) -> Self::PointIterator<'_> {
        panic!();
    }
    fn vertices(&self) -> Self::VertexIterator<'_> {
        panic!();
    }
}

impl<T: Float, GridImpl: Grid<T = T>> GridType for GridImpl {
    type T = T;
    //  PROPOSAL:
    //  Vertex = one of the corners of a cell
    //  Point = a point in the geometry

    type Point<'a> = Point<'a, T, GridImpl::Geometry> where Self: 'a;
    type Cell<'a> = Cell<'a, T, GridImpl> where Self: 'a;
    type Edge = ();
    type Face = ();

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
}

#[cfg(test)]
mod test {
    use crate::grid_impl::grid::*;
    use crate::grid_impl::mixed_grid::SerialMixedGrid;

    #[test]
    fn test_grid_mixed_cell_type() {
        let grid = SerialMixedGrid::<f64>::new(
            vec![
                -1.0,
                0.0,
                0.0,
                -0.5,
                0.0,
                0.2,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                2.0,
                0.0,
                0.0,
                -std::f64::consts::FRAC_1_SQRT_2,
                std::f64::consts::FRAC_1_SQRT_2,
                0.0,
                0.0,
                0.5,
                0.0,
                0.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
            ],
            3,
            &[0, 2, 7, 6, 5, 1, 2, 3, 7, 8, 3, 4, 8],
            &[
                ReferenceCellType::Triangle,
                ReferenceCellType::Quadrilateral,
                ReferenceCellType::Triangle,
            ],
            &[2, 1, 1],
        );

        let mut coords = vec![0.0; grid.geometry().dim()];
        for vertex in grid.iter_all_points() {
            vertex.coords(coords.as_mut_slice());
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
        // assert_eq!(0, 1);
    }
}
