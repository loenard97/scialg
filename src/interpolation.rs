//! Interpolation (and extrapolation) of errorless datasets

use num::Float;

use crate::{function::Polynomial, sort::co_sort};

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

/// Linear spline interpolation
pub struct LinearSplines {
    xs: Vec<f64>,
    ys: Vec<f64>,
    splines: Vec<Polynomial<2>>,
}

impl LinearSplines {
    pub fn new(xs: &[f64], ys: &[f64]) -> Self {
        assert_eq!(xs.len(), ys.len());
        let n = xs.len();

        let mut ls = LinearSplines {
            xs: xs.to_vec(),
            ys: ys.to_vec(),
            splines: Vec::new(),
        };

        co_sort(&mut ls.xs, &mut ls.ys);

        for i in 1..n {
            ls.splines.push(Polynomial::new([
                ys[i - 1],
                (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]),
            ]));
        }

        ls
    }

    pub fn eval(&self, x: f64) -> f64 {
        let idx = self.xs.binary_search_by(|v| v.partial_cmp(&x).unwrap());

        match idx {
            Ok(i) => self.ys[i],
            Err(i) => self.splines[i - 1].eval(x - self.xs[i - 1]),
        }
    }
}

/// Natural cubic spline interpolation
pub struct CubicSplines {
    xs: Vec<f64>,
    ys: Vec<f64>,
    splines: Vec<Polynomial<4>>,
}

impl CubicSplines {
    pub fn new(xs: &[f64], ys: &[f64]) -> Self {
        let n = xs.len();
        let a = ys.to_vec();
        let mut b = vec![0.0; n - 1];
        let mut d = vec![0.0; n - 1];

        let h: Vec<f64> = xs.windows(2).map(|x| x[1] - x[0]).collect();
        let mut alpha: Vec<f64> = vec![0.0; n - 1];
        for i in 1..n - 1 {
            alpha[i] = 3.0 / h[i] * (a[i + 1] - a[i]) - 3.0 / h[i - 1] * (a[i] - a[i - 1]);
        }

        let mut c = vec![0.0; n];
        let mut l = vec![0.0; n];
        let mut mu = vec![0.0; n];
        let mut z = vec![0.0; n];
        l[0] = 1.0;
        mu[0] = 0.0;
        z[0] = 0.0;

        for i in 1..n - 1 {
            l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        l[n - 1] = 1.0;
        z[n - 1] = 0.0;
        c[n - 1] = 0.0;

        for i in (0..n - 1).rev() {
            c[i] = z[i] - mu[i] * c[i + 1];
            b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
            d[i] = (c[i + 1] - c[i]) / 3.0 / h[i];
        }

        let mut vec_splices = Vec::new();
        for i in 0..n - 1 {
            let mut coeffs: Vec<f64> = vec![0.0; 4];
            coeffs[0] = a[i];
            coeffs[1] = b[i];
            coeffs[2] = c[i];
            coeffs[3] = d[i];

            let splice: Polynomial<4> = Polynomial::from_slice(&coeffs[..]);
            vec_splices.push(splice);
        }

        CubicSplines {
            xs: xs.to_vec(),
            ys: ys.to_vec(),
            splines: vec_splices,
        }
    }

    pub fn eval(&self, x: f64) -> f64 {
        let idx = self.xs.binary_search_by(|v| v.partial_cmp(&x).unwrap());

        match idx {
            Ok(i) => self.ys[i],
            Err(i) => self.splines[i - 1].eval(x - self.xs[i - 1]),
        }
    }
}

#[cfg(test)]
mod tests {
    use ndarray::Array1;

    use super::*;

    #[test]
    fn test_linear_splines() {
        let xs = Array1::linspace(0.0, 4.0 * std::f64::consts::PI, 10);
        let ys = xs.map(|x| f64::sin(*x));

        let linear_spline = LinearSplines::new(&xs.as_slice().unwrap(), &ys.as_slice().unwrap());

        for i in 0..6000 {
            let x: f64 = i as f64 * 0.001;
            println!("{} {}", x, linear_spline.eval(x));
        }
    }

    #[test]
    fn test_cubic_splines() {
        let xs = Array1::linspace(0.0, 4.0 * std::f64::consts::PI, 8);
        let ys = xs.map(|x| f64::sin(*x));

        let cubic_spline = CubicSplines::new(&xs.as_slice().unwrap(), &ys.as_slice().unwrap());

        for i in 0..6000 {
            let x: f64 = i as f64 * 0.001;
            println!("{} {}", x, cubic_spline.eval(x));
        }
    }
}
