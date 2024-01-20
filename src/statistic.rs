//! Statistical modeling of datasets

/// Calculate the mean of a dataset *xs*
pub fn mean(xs: &[f64]) -> f64 {
    let n = xs.len() as f64;

    xs.iter().sum::<f64>() / n
}

/// Calculate the *n*-th moment of a dataset *xs*
///
/// # References
///  - [Wikiepdia: Moment](https://en.wikipedia.org/wiki/Moment_(mathematics))
pub fn moment(xs: &[f64], n: i32) -> f64 {
    let len = xs.len() as f64;
    let mean = mean(xs);

    xs.iter().map(|x| (x - mean).powi(n)).sum::<f64>() / len
}
