//! Grid implementation

pub mod common;
pub mod flat_triangle_grid;
pub mod grid;
pub mod mixed_grid;
pub mod single_element_grid;
pub mod traits;

pub use crate::grid::traits::Geometry;
pub use crate::grid::traits::Grid;
pub use crate::grid::traits::Topology;
