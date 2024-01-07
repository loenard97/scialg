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
        sum += h * (0.5 * f(xn) + 0.5 * f(xn + h));
        xn += h;
    }

    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn test_trapezoid_panic() {
        trapezoid(f32::sin, 1.0, 0.0, -0.01);
    }
}
