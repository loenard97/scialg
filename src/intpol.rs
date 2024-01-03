/// Interpolation (and extrapolation) algorithms
use num::Float;

/// Interpolate a dataset *f(xs)->ys* at point *x* using Neville interpolation
///
/// # Example
/// ```
/// use scialg::intpol::neville;
///
/// let xs = vec![-1.0, 0.0, 1.0];
/// let ys = vec![2.0, 0.0, 2.0];
///
/// assert_eq!(neville(&xs, &ys, 1.0), 2.0);
/// ```
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
