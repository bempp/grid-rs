//! Map from reference to physical space.

use super::GridType;

pub trait ReferenceMapType {
    type Grid: GridType;

    /// Return the domain dimension.
    fn domain_dimension(&self) -> usize;

    /// Return the physical dimension.
    fn physical_dimension(&self) -> usize;

    fn number_of_reference_points(&self) -> usize;

    /// Defines an iterator that returns a slice with the value of the
    /// physical point for each reference point.
    fn reference_to_physical(&self, cell_index: usize, point_index: usize, value: &mut [<Self::Grid as GridType>::T]);

    /// Defines an iterator that returns a slice with the value of the
    /// Jacobian at the physical point for each reference point.
    fn jacobian(&self, cell_index: usize, point_index: usize, value: &mut [<Self::Grid as GridType>::T]);

    /// Defines an iterator that returns a slice with the normal direction
    /// at each point.
    fn normal(&self, cell_index: usize, point_index: usize, value: &mut [<Self::Grid as GridType>::T]);
}
