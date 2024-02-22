//! Map from reference to physical space.

use num::Float;

pub trait Map {
    type T: Float;

    fn domain_dimension(&self);

    fn range_dimension(&self);

    fn map_to_physical(&self, coord: &[Self::T], out: &mut [Self::T]);
    fn gradient(&self, coord: &[Self::T], out: &mut [Self::T]);
}
