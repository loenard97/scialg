//! Evaluation algorithms for functions

use num::Complex;

/// Polynomial of degree N-1
#[derive(Debug, Copy, Clone)]
pub struct Polynomial<const N: usize> {
    coeff: [f64; N],
}

impl<const N: usize> Polynomial<N> {
    pub fn new(coeff: [f64; N]) -> Self {
        Polynomial { coeff }
    }

    pub fn from_slice(coeff: &[f64]) -> Self {
        assert_eq!(coeff.len(), N);

        let mut arr = [0.0; N];

        for i in 0..coeff.len() {
            arr[i] = coeff[i];
        }

        Polynomial { coeff: arr }
    }

    /// Evaluate polynomial at *x* using Horner's method
    ///
    /// # References
    ///  - [Wikipedia: Horner's method](https://en.wikipedia.org/wiki/Horner%27s_method)
    pub fn eval(&self, x: f64) -> f64 {
        let mut val = 0.0;

        for i in (0..N).rev() {
            val = self.coeff[i] + x * val;
        }

        val
    }
}

/// Calcuate the gamma function of z0 using the Lanczos approximation
///
/// # Example
/// ```
/// use num::Complex;
///
/// use scialg::function::gamma;
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
///
/// # References
///  - [Wikipedia: Lanczos approximation](https://en.wikipedia.org/wiki/Lanczos_approximation)
///  - [Wolfram MathWorld: Lanczos
///  approximation](https://mathworld.wolfram.com/LanczosApproximation.html)
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

    f64::sqrt(2.0 * pi) * t.powc(z + 0.5) * (-t).exp() * x
}

/// Calcuate the beta function of (z, w)
///
/// The beta function can be computed as: beta(z, w) = gamma(z) * gamma(w) / gamma(z + w)
///
/// # References
///  - [Wikipedia: Beta function](https://en.wikipedia.org/wiki/Beta_function)
pub fn beta(z: Complex<f64>, w: Complex<f64>) -> Complex<f64> {
    gamma(z) * gamma(w) / gamma(z + w)
}

/// Calcuate the factorial n!
///
/// # Examples
/// ```
/// use scialg::function::factorial;
///
/// assert_eq!(factorial(0), 1);
/// assert_eq!(factorial(1), 1);
/// assert_eq!(factorial(2), 2);
/// assert_eq!(factorial(3), 6);
/// assert_eq!(factorial(4), 24);
/// assert_eq!(factorial(5), 120);
/// assert_eq!(factorial(6), 720);
/// ```
///
/// # References
///  - [Wikipedia: Factorial](https://en.wikipedia.org/wiki/Factorial)
pub fn factorial(n: i64) -> i64 {
    let z: Complex<f64> = Complex::new((n + 1) as f64, 0.0);
    let res = gamma(z);

    res.re.round() as i64
}

/// Calculate the binomial distribution of (n k)
///
/// # Examples
/// ```
/// use scialg::function::binomial;
///
/// assert_eq!(binomial(0, 0), 1);
/// assert_eq!(binomial(0, 1), 0);
/// assert_eq!(binomial(1, 0), 1);
/// assert_eq!(binomial(12, 6), 924);
/// assert_eq!(binomial(3, 5), 0);
/// assert_eq!(binomial(5, 3), 10);
/// ```
///
/// # References
///  - [Wikipedia: Binomial distribution](https://en.wikipedia.org/wiki/Binomial_distribution)
pub fn binomial(n: i64, k: i64) -> i64 {
    if k > n {
        return 0;
    }

    factorial(n) / (factorial(k) * factorial(n - k))
}

/// Calculate the n-th element of the Fibonacci sequence
///
/// # Examples
/// ```
/// use scialg::function::fibonacci;
///
/// assert_eq!(fibonacci(0), 0.0);
/// assert_eq!(fibonacci(1), 1.0);
/// assert_eq!(fibonacci(2), 1.0);
/// assert_eq!(fibonacci(3), 2.0);
/// assert_eq!(fibonacci(4), 3.0);
/// assert_eq!(fibonacci(5), 5.0);
/// assert_eq!(fibonacci(20), 6765.0);
/// ```
///
/// # References
///  - [Wikipedia: Fibonacci sequence](https://en.wikipedia.org/wiki/Fibonacci_sequence#Computation_by_rounding)
pub fn fibonacci(n: i32) -> f64 {
    let phi = 1.618033988749894848204586834365638118_f64;

    (phi.powi(n) / 5.0_f64.sqrt()).round().try_into().unwrap()
}

/// Calculate the natural logarithm *ln(x)*
pub fn ln(x: f64, iter: i32) -> f64 {
    let mut sum = 0.0;

    for i in 0..iter {
        sum += (-1.0_f64).powi(i + 1) * (x - 1.0).powi(i) / i as f64;
    }

    sum
}

/// Calculate the square root of *x*
pub fn sqrt(x: f64, iter: usize) -> f64 {
    // Newton Method
    // f(x) = y^2 - x
    // f'(x) = 2*y
    let mut y = 1.0;

    for _ in 0..iter {
        y = y - (y * y - x) / (2.0 * y);
    }

    y
}

pub fn pow(x: f64, n: i32) -> f64 {
    if n == 0 {
        return 1.0;
    }

    let t = pow(x, n / 2);

    if n % 2 == 0 {
        t * t
    } else {
        x * t * t
    }
}

pub fn sin(x: f64, iter: i32) -> f64 {
    let pi_2 = 2.0 * std::f64::consts::PI;
    let mut y = 0.0;
    let mut x = x;

    while x > pi_2 {
        x -= pi_2;
    }

    for i in (1..iter).step_by(2) {
        y += (-1.0_f64).powf(i as f64 / 2.0) / factorial(i as i64) as f64 * x.powi(i);
    }

    y
}

pub fn cos(x: f64, iter: i32) -> f64 {
    let pi_2 = 2.0 * std::f64::consts::PI;
    let mut y = 1.0;
    let mut x = x;

    while x > pi_2 {
        x -= pi_2;
    }

    for i in (2..iter).step_by(2) {
        y += (-1.0_f64).powf(i as f64 / 2.0) / factorial(i as i64) as f64 * x.powi(i);
    }

    y
}
