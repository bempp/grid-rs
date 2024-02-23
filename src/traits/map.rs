//! Map from reference to physical space.

use num::Float;

pub trait Map {
    type T: Float;

    fn domain_dimension(&self) -> usize;

    fn range_dimension(&self) -> usize;

    fn map(&self, coord: &[Self::T], out: &mut [Self::T]);

    fn map_from_index(&self, index: usize, out: &mut [Self::T]);

    fn gradient(&self, coord: &[Self::T], out: &mut [Self::T]);

    fn gradient_from_index(&self, index: usize, out: &mut [Self::T]);
}
