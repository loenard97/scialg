//! Filter onedimensional data

/// Compute the moving average of *xs* with a window size of *k*
///
/// # Reference
///  - [Wikipedia: Moving average](https://en.wikipedia.org/wiki/Moving_average)
pub fn moving_average(xs: &[f64], k: usize) -> Vec<f64> {
    xs.windows(k)
        .map(|w| w.iter().sum::<f64>() / k as f64)
        .collect()
}
