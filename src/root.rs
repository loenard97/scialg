use std::mem::swap;

/// Find the root of f in [a, b] using the Bisection method
///
/// # References
///  - [Wikipedia](https://en.wikipedia.org/wiki/Bisection_method)
pub fn bisection(f: fn(f32) -> f32, a: f32, b: f32, epsilon: Option<f32>) -> f32 {
    let (mut rbr, mut lbr) = if f(a) < 0.0 && f(b) > 0.0 {
        (b, a)
    } else {
        (a, b)
    };

    loop {
        let mid = (rbr + lbr) / 2.0;
        if f(mid) > 0.0 {
            rbr = mid;
        } else {
            lbr = mid;
        }

        if (rbr - lbr).abs() < epsilon.unwrap() {
            return (rbr + lbr) / 2.0;
        }
    }
}

/// Find the root of f in [a, b] using the Secant method.
///
/// # References
///  - [Wikipedia](https://en.wikipedia.org/wiki/Secant_method)
pub fn secant(f: fn(f32) -> f32, a: f32, b: f32, epsilon: Option<f32>) -> f32 {
    let mut fl = f(a);
    let mut fr = f(b);

    let (mut xl, mut rts) = if fl.abs() < fr.abs() {
        swap(&mut fl, &mut fr);
        (b, a)
    } else {
        (a, b)
    };

    loop {
        let dx = (xl - rts) * fr / (fr - fl);
        xl = rts;
        fl = fr;
        rts += dx;
        fr = f(rts);

        if dx.abs() < epsilon.unwrap() {
            return rts;
        }
    }
}

/// Find the root of f in [a, b] using the Regula falsi method.
///
/// # References
///  - [Wikipedia](https://en.wikipedia.org/wiki/Regula_falsi)
pub fn regula_falsi(
    f: fn(f64) -> f64,
    a0: f64,
    b0: f64,
    epsilon: Option<f64>,
    max_iter: Option<usize>,
) -> Option<f64> {
    let mut a = a0;
    let mut b = b0;

    for _ in 0..max_iter.unwrap_or(20) {
        let fa = f(a);
        let fb = f(b);
        let c = (a * fb - b * fa) / (fb - fa);
        let fc = f(c);

        if fa * fc > 0.0 {
            a = c;
        } else {
            b = c;
        }

        if (a - b).abs() < epsilon.unwrap_or(1e-10) {
            return Some(c);
        }
    }

    None
}

/// Find the root of f in [a, b] using Ridder's method.
///
/// # References
///  - [Wikipedia](https://en.wikipedia.org/wiki/Ridders'_method)
pub fn ridder(
    f: fn(f64) -> f64,
    a0: f64,
    b0: f64,
    epsilon: Option<f64>,
    max_iter: Option<usize>,
) -> Option<f64> {
    let mut a = a0;
    let mut b = b0;

    for _ in 0..max_iter.unwrap_or(20) {
        let fa = f(a);
        let fb = f(b);
        let c = (a + b) / 2.0;
        let fc = f(c);
        let d = c + (c - a) * (fa - fb).signum() * fc / (fc.powi(2) - fa * fb).sqrt();
        let fd = f(d);

        if fa * fd > 0.0 {
            a = d;
        } else {
            b = d;
        }

        if fd.abs() < epsilon.unwrap_or(1e-10) {
            return Some(d);
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bisection() {
        let epsilon = 1e-6;
        let root = bisection(f32::sin, 1.5, 4.5, Some(epsilon));

        assert!((root - std::f32::consts::PI).abs() < epsilon);
    }

    #[test]
    fn test_secant() {
        let epsilon = 1e-6;
        let root = secant(f32::sin, 1.5, 4.5, Some(epsilon));

        assert!((root - std::f32::consts::PI).abs() < epsilon);
    }

    #[test]
    fn test_regula_falsi() {
        let epsilon = 1e-15;
        let pi = std::f64::consts::PI;
        let root = regula_falsi(f64::sin, 1.5, 4.5, Some(epsilon), None).expect("did not converge");

        assert!((root - pi).abs() < epsilon);
    }

    #[test]
    fn test_ridder() {
        let epsilon = 1e-10;
        let pi = std::f64::consts::PI;
        let root = ridder(f64::sin, 1.5, 4.5, Some(epsilon), None).expect("did not converge");

        assert!((root - pi).abs() < epsilon);
    }
}
