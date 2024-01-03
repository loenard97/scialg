#![allow(dead_code)]

use num::Complex;

/// Calcuate the gamma function of z0 using the Lanczos approximation
///
/// # Example
/// ```
/// use num::Complex;
///
/// use scialg::func::gamma;
///
/// let tol = 1e-3;
/// let z = Complex::new(1.0, 0.0);
/// let res = gamma(z);
///
/// assert!((gamma(Complex::new(1.0, 0.0)) - Complex::new(1.0, 0.0)).norm() < tol);
/// assert!((gamma(Complex::new(0.5, 0.0)) - Complex::new(1.7724, 0.0)).norm() < tol);
/// assert!((gamma(Complex::new(1.0, 0.0)) - Complex::new(1.0, 0.0)).norm() < tol);
/// assert!((gamma(Complex::new(5.0, 0.0)) - Complex::new(24.0, 0.0)).norm() < tol);
/// ```
pub fn gamma(z0: Complex<f64>) -> Complex<f64> {
    let pi = std::f64::consts::PI;
    let coeffs = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    if z0.re < 0.5 {
        return pi / ((pi * z0).sin() * gamma(1.0 - z0));
    }

    let z = z0 - 1.0;
    let mut x = Complex::new(coeffs[0], 0.0);
    for i in 1..coeffs.len() {
        x += coeffs[i] / (z + i as f64);
    }
    let t = z + 7.5;
    return f64::sqrt(2.0 * pi) * t.powc(z + 0.5) * (-t).exp() * x;
}

/// Calcuate the beta function of (z, w)
pub fn beta(z: Complex<f64>, w: Complex<f64>) -> Complex<f64> {
    gamma(z) * gamma(w) / gamma(z + w)
}

/// Calcuate the factorial n!
///
/// # Examples
/// ```
/// use scialg::func::factorial;
///
/// assert_eq!(factorial(0), 1);
/// assert_eq!(factorial(1), 1);
/// assert_eq!(factorial(2), 2);
/// assert_eq!(factorial(3), 6);
/// assert_eq!(factorial(4), 24);
/// assert_eq!(factorial(5), 120);
/// assert_eq!(factorial(6), 720);
/// ```
pub fn factorial(n: i64) -> i64 {
    let z: Complex<f64> = Complex::new((n + 1) as f64, 0.0);
    let res = gamma(z);

    return res.re.round() as i64;
}

/// Calculate the binomial distribution of (n k)
///
/// # Examples
/// ```
/// use scialg::func::binomial;
///
/// assert_eq!(binomial(0, 0), 1);
/// assert_eq!(binomial(0, 1), 0);
/// assert_eq!(binomial(1, 0), 1);
/// assert_eq!(binomial(12, 6), 924);
/// assert_eq!(binomial(3, 5), 0);
/// assert_eq!(binomial(5, 3), 10);
/// ```
pub fn binomial(n: i64, k: i64) -> i64 {
    if k > n {
        return 0;
    }
    factorial(n) / (factorial(k) * factorial(n - k))
}
