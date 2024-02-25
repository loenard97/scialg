//! Interpolation (and extrapolation) of errorless datasets

use num::Float;

/// Interpolate a dataset *f(xs)->ys* at point *x* using Neville interpolation
///
/// # Example
/// ```
/// use scialg::interpolation::neville;
///
/// let xs = vec![-1.0, 0.0, 1.0];
/// let ys = vec![2.0, 0.0, 2.0];
///
/// assert_eq!(neville(&xs, &ys, 1.0), 2.0);
/// ```
///
/// # References
///  - [Wikipedia: Neville's Algorithm](https://en.wikipedia.org/wiki/Neville%27s_algorithm)
///  - [Wolfram MathWorld: Neville's Algorithm](https://mathworld.wolfram.com/NevillesAlgorithm.html)
pub fn neville<F: Float>(xs: &[F], ys: &[F], x: F) -> F {
    let n = xs.len();
    let mut q = ys.to_vec();

    for k in 1..n {
        for i in 0..n - k {
            q[i] = ((x - xs[i + k]) * q[i] + (xs[i] - x) * q[i + 1]) / (xs[i] - xs[i + k]);
        }
    }

    q[0]
}

/// Linear interpolation of dataset *f(xs) -> ys* at point *x*
///
/// # Example
/// ```
/// use scialg::interpolation::linear;
///
/// let xs = vec![0.0, 2.0];
/// let ys = vec![1.0, 0.5];
///
/// assert_eq!(linear(&xs, &ys, 1.0), 0.75);
/// ```
pub fn linear<F: Float>(xs: &[F], ys: &[F], x: F) -> F {
    let idx = xs.binary_search_by(|v| v.partial_cmp(&x).unwrap());

    match idx {
        Ok(i) => ys[i],
        Err(i) => ys[i - 1] + x * (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interpolation_linear() {
        let xs = vec![0.0, 2.0];
        let ys = vec![1.0, 2.0];

        assert_eq!(linear(&xs, &ys, 1.0), 1.5);
    }
}
