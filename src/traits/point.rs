//! Definition of a vertex

use crate::types::Float;

pub trait PointType {
    //! A point

    /// The floating point type used for coordinates
    type T: Float;

    /// Get the coordinates of the point
    fn coords(&self, data: &mut [Self::T]);

    // Get the point's index
    fn index(&self) -> usize;

    // Get the point's id
    fn id(&self) -> usize;
}
