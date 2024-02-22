//! Definition of an element.

use num::Float;

use super::{cell::Cell, map::Map};

pub trait Element {
    type T: Float;

    type CellType: Cell<T = Self::T>;
    type CellMap: Map<T = Self::T>;
    type ElementMap: Map<T = Self::T>;

    fn cell_map(&self) -> &Self::CellMap;
    fn cell(&self) -> &Self::CellType;
    fn element_map(&self) -> &Self::ElementMap;

    fn local_index(&self) -> usize;

    fn global_index(&self) -> usize;
}
