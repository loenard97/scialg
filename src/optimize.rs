/// Find a local minimum of *f* in *[a, c]* using the Golden-section search
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

    return (b + a) / 2.0;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_golden_section() {
        let pi = std::f32::consts::PI;
        let tolerance = 1e-4_f32;
        let min = golden_section(f32::sin, 3.0, 6.0, tolerance);

        assert!((min - 1.5 * pi).abs() < 5.0 * tolerance);
    }
}
