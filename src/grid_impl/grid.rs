use crate::grid_impl::traits::{Geometry, Grid, Topology};
use crate::traits::{
    cell::CellType, geometry::GeometryType, grid::GridType, topology::TopologyType,
    vertex::VertexType,
};
use num::Float;
use std::iter::Copied;

pub struct Vertex<'a, T: Float, G: Geometry> {
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
    index: usize,
}
pub struct CellGeometry<'a, T: Float, GridImpl: Grid> {
    grid: &'a GridImpl,
    index: usize,
    _t: std::marker::PhantomData<T>,
}

impl<'a, T: Float, G: Geometry<T = T>> VertexType for Vertex<'a, T, G> {
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
    type VertexIndexIter<'a> = Copied<std::slice::Iter<'a, usize>>
    where
        Self: 'a;

    type EdgeIndexIter<'a> = Copied<std::slice::Iter<'a, usize>>
    where
        Self: 'a;

    fn vertex_indices(&self) -> Self::VertexIndexIter<'_> {
        self.grid
            .topology()
            .cell(self.index)
            .unwrap()
            .iter()
            .copied()
    }

    fn edge_indices(&self) -> Self::EdgeIndexIter<'_> {
        self.grid
            .topology()
            .cell(self.index)
            .unwrap()
            .iter()
            .copied()
    }
}

impl<'grid, T: Float, GridImpl: Grid<T = T>> GeometryType for CellGeometry<'grid, T, GridImpl>
where
    GridImpl: 'grid,
{
    type T = T;

    // This assumed that points are stored xyzxyzxyz. Would it be better to have
    // an iterator through the indices, then GridType::vertex_From_index
    type PointIterator<'a> = Copied<std::slice::Iter<'a, &'a [Self::T]>>
    where
        Self: 'a;

    fn physical_dimension(&self) -> usize {
        self.grid.geometry().dim()
    }

    fn midpoint(&self, point: &mut [Self::T]) {
        self.grid.geometry().midpoint(self.index, point)
    }

    fn diameter(&self) -> Self::T {
        self.grid.geometry().diameter(self.index)
    }

    fn volume(&self) -> Self::T {
        self.grid.geometry().volume(self.index)
    }

    fn corners(&self) -> Self::PointIterator<'_> {
        panic!();
    }
}

impl<T: Float, GridImpl: Grid<T = T>> GridType for GridImpl {
    //  PROPOSAL:
    //  Vertex = one of the corners of a cell
    //  Point = a point in the geometry

    type Vertex<'a> = Vertex<'a, T, GridImpl::Geometry> where Self: 'a;
    type Cell<'a> = Cell<'a, T, GridImpl> where Self: 'a;
    type Edge = ();
    type Face = ();

    //fn number_of_points(&self) -> usize {
    //    self.geometry().point_count()
    //}
    fn number_of_vertices(&self) -> usize {
        // self.topology().entity_types(self.topology().dim()).point_count()
        self.geometry().point_count()
    }
    fn number_of_cells(&self) -> usize {
        self.geometry().cell_count()
    }

    fn vertex_index_from_id(&self, id: usize) -> usize {
        id
    }
    fn vertex_id_from_index(&self, index: usize) -> usize {
        index
    }

    fn cell_index_from_id(&self, id: usize) -> usize {
        id
    }
    fn cell_id_from_index(&self, index: usize) -> usize {
        index
    }

    fn vertex_from_index(&self, index: usize) -> Self::Vertex<'_> {
        Self::Vertex {
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
    use crate::reference_cell::ReferenceCellType;

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
        for vertex in grid.iter_all_vertices() {
            vertex.coords(coords.as_mut_slice());
            println!("{:#?}", coords);
        }

        for cell in grid.iter_all_cells() {
            for (local_index, (vertex_index, edge_index)) in cell
                .topology()
                .vertex_indices()
                .zip(cell.topology().edge_indices())
                .enumerate()
            {
                println!(
                    "Cell: {}, Vertex: {}, {}, Edge: {}, {}, Volume: {}",
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
