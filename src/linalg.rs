#![allow(dead_code)]

use ndarray::Array2;
use num::Float;

/// Gauss-Jordan Elimination
pub fn gauss_jordan<F: Float>(arr: &mut Array2<F>) {
    for k in 0..arr.nrows() {
        let pivot = arr[(k, k)];
        arr.row_mut(k).mapv_inplace(|e| e / pivot);
        for j in 0..arr.nrows() {
            if j == k {
                continue;
            }
            let scale = arr[(j, k)];
            for i in 0..arr.ncols() {
                arr[(j, i)] = arr[(j, i)] - scale * arr[(k, i)];
            }
        }
    }
}

pub fn transpose<F: Float>(arr: &Array2<F>) -> Array2<F> {
    let ni = arr.shape()[0];
    let nj = arr.shape()[1];
    let mut res = Array2::zeros((nj, ni));

    for j in 0..ni {
        for i in 0..nj {
            res[(i, j)] = arr[(j, i)];
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array;

    #[test]
    fn test_gauss_jordan() {
        let mut input: Array2<f32> = Array::from_shape_vec(
            (2, 5),
            vec![-1.0, 2.0, -6.0, 1.0, 0.0, 3.0, -4.0, 14.0, 0.0, 1.0],
        )
        .unwrap();
        let output: Array2<f32> = Array::from_shape_vec(
            (2, 5),
            vec![1.0, 0.0, 2.0, 2.0, 1.0, 0.0, 1.0, -2.0, 1.5, 0.5],
        )
        .unwrap();

        gauss_jordan(&mut input);

        assert_eq!(input, output);
    }

    #[test]
    fn test_transpose() {
        let input: Array2<f32> =
            Array::from_shape_vec((2, 3), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).unwrap();
        let output: Array2<f32> =
            Array::from_shape_vec((3, 2), vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]).unwrap();

        assert_eq!(transpose(&input), output);
    }
}
