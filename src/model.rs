//! Fitting of datasets with known errors to a known model

use ndarray::prelude::*;

use crate::linalg::{invert, transpose};

/// Calculate a linear regression model `f(x) = b * x + a` for the data set (xs, ys), where xs are exact and ys have
/// standard deviations of yerrs
/// Returns (a, b, siga, sigb, chi2) where siga and sigb are standard deviations of a and b and
/// chi2 is Chi squared of the fitted model
///
/// # Example
/// ```
/// use scialg::model::linear_regression;
///
/// let xs = vec![0.0, 1.0, 2.0];
/// let ys = vec![1.0, 3.0, 5.0];
/// let sigs = vec![1e-3, 1e-3, 1e-3];
///
/// let (a, b, siga, sigb, chi2) = linear_regression(&xs, &ys, &sigs);
///
/// assert_eq!(a, 1.0);
/// assert_eq!(b, 2.0);
/// ```
pub fn linear_regression(xs: &[f64], ys: &[f64], yerrs: &[f64]) -> (f64, f64, f64, f64, f64) {
    let mut ss = 0.0;
    let mut sx = 0.0;
    let mut sy = 0.0;

    for i in 0..xs.len() {
        let wt = 1.0 / (yerrs[i] * yerrs[i]);
        ss += wt;
        sx += wt * xs[i];
        sy += wt * ys[i];
    }

    let sxoss = sx / ss;

    let mut st2 = 0.0;
    let mut b = 0.0;
    for i in 0..xs.len() {
        let t = (xs[i] - sxoss) / yerrs[i];
        st2 += t * t;
        b += t * ys[i] / yerrs[i];
    }
    b /= st2;
    let a = (sy - sx * b) / ss;
    let siga = ((1.0 + sx * sx / (ss * st2)) / ss).sqrt();
    let sigb = (1.0 / st2).sqrt();

    let mut chi2 = 0.0;
    for i in 0..xs.len() {
        chi2 += ((ys[i] - a - b * xs[i]) / yerrs[i]).powi(2);
    }

    (a, b, siga, sigb, chi2)
}

/// Compute solution to the general linear least squares problem
///
/// # Reference
///  - [Wikipedia: Generalized least squares](https://en.wikipedia.org/wiki/Generalized_least_squares)
pub fn general_linear_least_squares(
    xs: &[f64],
    ys: &[f64],
    fns: &Vec<fn(f64) -> f64>,
) -> (Vec<f64>, Vec<f64>) {
    let n = xs.len();
    let m = fns.len();
    let ys: Array1<f64> = Array1::from(ys.to_vec());

    let design: Array2<f64> = Array2::from_shape_fn((n, m), |(j, i)| fns[i](xs[j]));
    let design_transposed = transpose(&design);
    let covariance = invert(&design_transposed.dot(&design));

    let res = covariance.dot(&design_transposed);
    let vals = res.dot(&ys);

    let vals: Vec<_> = vals.iter().cloned().collect();
    let vars: Vec<_> = covariance.diag().iter().cloned().collect();

    (vals, vars)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_general_linear_least_squares() {
        let xs = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        let ys = vec![3.0, 1.0, 1.0, 3.0, 7.0];

        let f0: fn(f64) -> f64 = |_| 1.0;
        let f1: fn(f64) -> f64 = |x| x;
        let f2: fn(f64) -> f64 = |x| x * x;
        let fs = vec![f0, f1, f2];

        let (coeffs, _vars) = general_linear_least_squares(&xs, &ys, &fs);

        for c in coeffs {
            assert!((c - 1.0) < 1e-5);
        }
    }
}
