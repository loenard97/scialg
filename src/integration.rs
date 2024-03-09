//! Evaluation algorithms for integrals of functions

use std::mem::swap;

/// Evaluate the integral of f in [a, b] with step size h using the Trapezoid rule
///
/// # Panics
/// Panics when b < a or h <= 0.
///
/// # Examples
/// ```
/// use scialg::integration::trapezoid;
///
/// let area = trapezoid(f32::sin, 0.0, std::f32::consts::PI, 0.001);
///
/// assert!((2.0 - area).abs() < 0.001);
/// ```
pub fn trapezoid(f: fn(f32) -> f32, a: f32, b: f32, h: f32) -> f32 {
    assert!(a <= b);
    assert!(h > 0.0);

    let mut sum = 0.0;
    let mut xn = a;
    while xn < b {
        sum += 0.5 * h * (f(xn) + f(xn + h));
        xn += h;
    }

    sum
}

/// Evaluate the integral of f in [a, b] using Romberg integration
///
/// # Examples
/// ```
/// use scialg::integration::romberg;
///
/// let area = romberg(f64::sin, 0.0, std::f64::consts::PI, 1e-5, None);
///
/// assert!((2.0 - area).abs() < 0.001);
/// ```
///
/// # References
///  - [Wikipedia: Romberg's method](https://en.wikipedia.org/wiki/Romberg%27s_method)
pub fn romberg(f: fn(f64) -> f64, a: f64, b: f64, acc: f64, max_steps: Option<usize>) -> f64 {
    let max_steps = max_steps.unwrap_or(100);

    let mut rp = vec![0.0; max_steps];
    let mut rc = vec![0.0; max_steps];

    let mut h = b - a;
    rp[0] = 0.5 * h * (f(a) + f(b));

    for i in 1..max_steps {
        h /= 2.0;
        let mut c = 0.0;
        let ep = 2_i32.pow((i - 1) as u32);
        for j in 1..ep + 1 {
            c += f(a + (2 * j - 1) as f64 * h);
        }
        rc[0] = h * c + 0.5 * rp[0];

        for j in 1..i + 1 {
            let nk = 4_i32.pow(j as u32);
            rc[j] = (nk as f64 * rc[j - 1] - rp[j - 1]) / (nk - 1) as f64;
        }

        if i > 2 && (rp[i - 1] - rc[i]).abs() < acc {
            return rc[i];
        }

        swap(&mut rp, &mut rc);
    }

    rp[max_steps - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn test_trapezoid_panic() {
        trapezoid(f32::sin, 1.0, 0.0, -0.01);
    }

    #[test]
    fn test_romberg() {
        let area = romberg(
            |x: f64| x.powi(4) * (x + (x * x + 1.0).sqrt()).log(10.0),
            0.0,
            2.0,
            1e-10,
            None,
        );

        assert!((area - 3.5410).abs() < 1e-3);
    }
}
