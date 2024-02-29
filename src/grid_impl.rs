//! Grid implementation

pub mod geometry;
pub mod grid;
pub mod serial_grid;
pub mod topology;
pub mod traits;

pub use crate::grid_impl::traits::Geometry;
pub use crate::grid_impl::traits::Grid;
pub use crate::grid_impl::traits::Topology;
