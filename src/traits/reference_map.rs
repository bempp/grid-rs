//! Map from reference to physical space.

use num::Float;

pub trait ReferenceMap {
    type T: Float;

    type SliceIterator<'a>: std::iter::Iterator<Item = &'a [Self::T]>
    where
        Self: 'a;

    /// Return the domain dimension.
    fn domain_dimension(&self) -> usize;

    /// Return the physical dimension.
    fn physical_dimension(&self) -> usize;

    /// Defines an iterator that returns a slice with the value of the
    /// physical point for each reference point.
    fn reference_to_physical(&self) -> Self::SliceIterator<'_>;

    /// Defines an iterator that returns a slice with the value of the
    /// Jacobian at the physical point for each reference point.
    fn jacobian(&self) -> Self::SliceIterator<'_>;

    /// Defines an iterator that returns a slice with the normal direction
    /// at each point.
    fn normal(&self) -> Self::SliceIterator<'_>;
}
