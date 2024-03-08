//! Minimization of functions

use num::Zero;

use crate::sort::co_sort;
use crate::vector::Vector;

/// Find a local minimum of *f* in *[a, c]* using the Golden-section search
///
/// # Example
/// ```
/// use scialg::optimize::golden_section;
///
/// let min = golden_section(f32::sin, 3.0, 6.0, 1e-4);
///
/// assert!((min - 1.5 * std::f32::consts::PI).abs() < 5e-4);
/// ```
///
/// # References
///  - [Wikipedia: Golden-section search](https://en.wikipedia.org/wiki/Golden-section_search)
pub fn golden_section(f: fn(f32) -> f32, a: f32, c: f32, tol: f32) -> f32 {
    let phi = 1.618033988749894848204586834365638118_f32;

    let mut a = a;
    let mut b = c;

    while (b - a).abs() > tol {
        let c = b - (b - a) / phi;
        let d = a + (b - a) / phi;
        if f(c) < f(d) {
            b = d;
        } else {
            a = c;
        }
    }

    (b + a) / 2.0
}

/// Find a local minima of *f* using the Nelder-Mead algorithm starting at *p0*
///
/// # Example
/// ```
/// use num::Zero;
///
/// use scialg::vector::Vector;
/// use scialg::optimize::nelder_mead;
///
/// let func = |v: Vector<2>| v.coeff[0] * v.coeff[0] + v.coeff[1] * v.coeff[1];
/// let minima = nelder_mead(func, &Vector::new(&[3.0, 5.0]), 100);
///
/// assert!((minima - Vector::zero()).length() < 1e-5);
/// ```
///
/// # References
///  - [Wikipedia: Nelder-Mead method](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)
///  - [Wikipedia: Centroid](https://en.wikipedia.org/wiki/Centroid)
pub fn nelder_mead<const N: usize>(
    f: fn(Vector<N>) -> f64,
    p0: &Vector<N>,
    max_iter: usize,
) -> Vector<N> {
    let alpha = 1.0;
    let gamma = 2.0;
    let rho = 0.5;
    let sigma = 0.5;

    // starting simplex
    let mut xs: Vec<Vector<N>> = vec![*p0];
    for i in 0..N {
        let mut coeffs = p0.coeff;
        coeffs[i] += 1.0;
        xs.push(Vector { coeff: coeffs });
    }

    for _ in 0..max_iter {
        let mut fs: Vec<f64> = xs.iter().map(|x| f(*x)).collect();
        co_sort(&mut fs, &mut xs);

        // compute centoid
        let mut x0 = Vector::zero();
        for i in 0..N {
            x0 = x0 + xs[i];
        }
        x0 = x0 / N as f64;

        // reflect
        let xr = x0 + (x0 - xs[N]) * alpha;
        let fr = f(xr);
        if fs[0] <= fr && fr <= fs[N - 1] {
            xs[N] = xr;
            continue;
        }

        // expand
        if fr < fs[0] {
            let xe = x0 + (xr - x0) * gamma;
            let fe = f(xe);
            if fe < fr {
                xs[N] = xe;
            } else {
                xs[N] = xr;
            }
            continue;
        }

        // contract
        if fr < fs[N] {
            // contract outside
            let xc = x0 + (xr - x0) * rho;
            let fc = f(xc);
            if fc < fr {
                xs[N] = xc;
                continue;
            }
        } else {
            // contract inside
            let xc = x0 + (xs[N] - x0) * rho;
            let fc = f(xc);
            if fc < fs[N] {
                xs[N] = xc;
                continue;
            }
        }

        // shrink
        for i in 0..xs.len() {
            xs[i] = xs[0] + (xs[i] - xs[0]) * sigma;
        }
    }

    // compute centoid
    let mut x0 = Vector::zero();
    for i in 0..N {
        x0 = x0 + xs[i];
    }
    x0 / N as f64
}
