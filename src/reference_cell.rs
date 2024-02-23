//! Information about reference cells

pub use crate::traits::cell::ReferenceCellType;
use num::Float;

/// The topological dimension of the cell
pub fn dim(cell: ReferenceCellType) -> usize {
    match cell {
        ReferenceCellType::Point => 0,
        ReferenceCellType::Interval => 1,
        ReferenceCellType::Triangle => 2,
        ReferenceCellType::Quadrilateral => 2,
        ReferenceCellType::Tetrahedron => 3,
        ReferenceCellType::Hexahedron => 3,
        ReferenceCellType::Prism => 3,
        ReferenceCellType::Pyramid => 3,
    }
}
/// Is the cell a simplex?
pub fn is_simplex(cell: ReferenceCellType) -> bool {
    match cell {
        ReferenceCellType::Point => true,
        ReferenceCellType::Interval => true,
        ReferenceCellType::Triangle => true,
        ReferenceCellType::Quadrilateral => false,
        ReferenceCellType::Tetrahedron => true,
        ReferenceCellType::Hexahedron => false,
        ReferenceCellType::Prism => false,
        ReferenceCellType::Pyramid => false,
    }
}

/// The vertices of the reference cell
pub fn vertices<T: Float>(cell: ReferenceCellType) -> Vec<T> {
    let zero = T::from(0.0).unwrap();
    let one = T::from(1.0).unwrap();
    match cell {
        ReferenceCellType::Point => vec![],
        ReferenceCellType::Interval => vec![zero, one],
        ReferenceCellType::Triangle => vec![zero, zero, one, zero, zero, one],
        ReferenceCellType::Quadrilateral => vec![zero, zero, one, zero, zero, one, one, one],
        ReferenceCellType::Tetrahedron => vec![
            zero, zero, zero, one, zero, zero, zero, one, zero, zero, zero, one,
        ],
        ReferenceCellType::Hexahedron => vec![
            zero, zero, zero, one, zero, zero, zero, one, zero, one, one, zero, zero, zero, one,
            one, zero, one, zero, one, one, one, one, one,
        ],
        ReferenceCellType::Prism => vec![
            zero, zero, zero, one, zero, zero, zero, one, zero, zero, zero, one, one, zero, one,
            zero, one, one,
        ],
        ReferenceCellType::Pyramid => vec![
            zero, zero, zero, one, zero, zero, zero, one, zero, one, one, zero, zero, zero, one,
        ],
    }
}
