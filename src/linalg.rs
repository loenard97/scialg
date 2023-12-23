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
}
