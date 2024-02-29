//! Information about reference cells

pub use crate::types::ReferenceCellType;
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

/// The edges of the reference cell
pub fn edges(cell: ReferenceCellType) -> Vec<Vec<usize>> {
    match cell {
        ReferenceCellType::Point => vec![],
        ReferenceCellType::Interval => vec![vec![0, 1]],
        ReferenceCellType::Triangle => vec![vec![1, 2], vec![0, 2], vec![0, 1]],
        ReferenceCellType::Quadrilateral => vec![vec![0, 1], vec![0, 2], vec![1, 3], vec![2, 3]],
        _ => {
            panic!("Not implemented yet");
        }
    }
}

/// The types of the subentities of the reference cell
pub fn entity_types(cell: ReferenceCellType) -> Vec<Vec<ReferenceCellType>> {
    match cell {
        ReferenceCellType::Point => vec![vec![ReferenceCellType::Point], vec![], vec![], vec![]],
        ReferenceCellType::Interval => vec![
            vec![ReferenceCellType::Point; 2],
            vec![ReferenceCellType::Interval],
            vec![],
            vec![],
        ],
        ReferenceCellType::Triangle => vec![
            vec![ReferenceCellType::Point; 3],
            vec![ReferenceCellType::Interval; 3],
            vec![ReferenceCellType::Triangle],
            vec![],
        ],
        ReferenceCellType::Quadrilateral => vec![
            vec![ReferenceCellType::Point; 4],
            vec![ReferenceCellType::Interval; 4],
            vec![ReferenceCellType::Quadrilateral],
            vec![],
        ],
        _ => {
            panic!("Not implemented yet");
        }
    }
}

/// The number of subentities of each dimension
pub fn entity_counts(cell: ReferenceCellType) -> Vec<usize> {
    match cell {
        ReferenceCellType::Point => vec![1, 0, 0, 0],
        ReferenceCellType::Interval => vec![2, 1, 0, 0],
        ReferenceCellType::Triangle => vec![3, 3, 1, 0],
        ReferenceCellType::Quadrilateral => vec![4, 4, 1, 0],
        _ => {
            panic!("Not implemented yet");
        }
    }
}

/// The connectivity of the reference cell
///
/// The indices of the result are [i][j][k][l]
pub fn connectivity(cell: ReferenceCellType) -> Vec<Vec<Vec<Vec<usize>>>> {
    match cell {
        ReferenceCellType::Point => vec![vec![vec![vec![0]]]],
        ReferenceCellType::Interval => vec![
            vec![vec![vec![0], vec![0]], vec![vec![1], vec![0]]],
            vec![vec![vec![0, 1], vec![0]]],
        ],
        ReferenceCellType::Triangle => vec![
            vec![
                vec![vec![0], vec![1, 2], vec![0]],
                vec![vec![1], vec![0, 2], vec![0]],
                vec![vec![2], vec![0, 1], vec![0]],
            ],
            vec![
                vec![vec![1, 2], vec![0], vec![0]],
                vec![vec![0, 2], vec![1], vec![0]],
                vec![vec![0, 1], vec![2], vec![0]],
            ],
            vec![vec![vec![0, 1, 2], vec![0, 1, 2], vec![0]]],
        ],
        ReferenceCellType::Quadrilateral => vec![
            vec![
                vec![vec![0], vec![0, 1], vec![0]],
                vec![vec![1], vec![0, 2], vec![0]],
                vec![vec![2], vec![1, 3], vec![0]],
                vec![vec![3], vec![2, 3], vec![0]],
            ],
            vec![
                vec![vec![0, 1], vec![0], vec![0]],
                vec![vec![0, 2], vec![1], vec![0]],
                vec![vec![1, 3], vec![2], vec![0]],
                vec![vec![2, 3], vec![3], vec![0]],
            ],
            vec![vec![vec![0, 1, 2, 3], vec![0, 1, 2, 3], vec![0]]],
        ],
        _ => {
            panic!("Not implemented yet");
        }
    }
}
