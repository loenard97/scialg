//! Evaluation algorithms for mathematical constants

use super::function;

/// Compute *pi* using the Gauss-Legendre algorithm
///
/// # Example
/// ```
/// use scialg::consts::pi_gauss_legendre;
///
/// let pi = pi_gauss_legendre(3);
/// assert_eq!(pi, 3.141592653589794);
/// ```
///
/// # References
///  - [Wikipedia: Gauss-Legendre algorithm](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_algorithm)
pub fn pi_gauss_legendre(iter: usize) -> f64 {
    let mut a = 1.0;
    let mut b = 1.0 / 2.0_f64.sqrt();
    let mut t = 0.25;
    let mut p = 1.0;

    for _ in 0..iter {
        let a_n = 0.5 * (a + b);
        b = (a * b).sqrt();
        t = t - p * (a - a_n) * (a - a_n);
        p = 2.0 * p;
        a = a_n;
    }
    return (a + b) * (a + b) / (4.0 * t);
}

/// Calculate *e* using a series
///
/// # References
///  - [Math is Fun: Eulers Number](https://www.mathsisfun.com/numbers/e-eulers-number.html)
pub fn e_series(iter: i32) -> f64 {
    (1.0 + 1.0 / (iter as f64)).powi(iter)
}

/// Calculate *e* using a series of factorials
///
/// # References
///  - [Math is Fun: Eulers Number](https://www.mathsisfun.com/numbers/e-eulers-number.html)
pub fn e_factorial(iter: i64) -> f64 {
    let mut sum = 1.0;

    for i in 0..iter {
        sum += 1.0 / function::factorial(i) as f64;
    }

    return sum;
}
