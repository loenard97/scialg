//! Fourier analysis

use ndarray::Array1;
use num::complex::Complex;

/// Replace *arr* with its Fourier transform
///
/// # Note
/// This algorithm takes O(N^2) computations.
/// Use the function scialg::fft::fft instead, which takes O(n log2(n)) computations.
///
/// # References
///  - [Fourier transforms (Harvard
///  lecture)](https://scholar.harvard.edu/files/schwartz/files/lecture8-fouriertransforms.pdf)
///  - [Wikipedia: Fourier transform](https://en.wikipedia.org/wiki/Fourier_transform)
///  - [Wikipedia: Discrete Fourier
///  transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform)
pub fn ft(arr: &Array1<f32>) -> Vec<Complex<f32>> {
    let pi = std::f32::consts::PI;
    let cn = arr.len() as f32;
    let mut res = Vec::new();

    for n in 0..arr.len() {
        let mut val = Complex::new(0.0, 0.0);
        for k in 0..arr.len() {
            let x = Complex::new(0.0, 2.0 * pi * k as f32 * n as f32 / cn);
            val += Complex::new(arr[k], 0.0) * x.exp();
        }

        res.push(val);
    }

    res
}

/// Bit reverse number *a* with a given *bit_size*
/// This is used in the FFT algorithm.
///
/// # Example
/// ```
/// use scialg::fourier::bit_reversed;
///
/// assert_eq!(bit_reversed(0b100, 3), 0b001);
/// assert_eq!(bit_reversed(0b110, 3), 0b011);
/// assert_eq!(bit_reversed(0b11000101, 8), 0b10100011);
/// assert_eq!(bit_reversed(0b1, 8), 0b10000000);
/// ```
pub fn bit_reversed(a: u8, bit_size: u8) -> u8 {
    let mut result = 0;
    let mut n = a;

    for _ in 0..bit_size {
        result = (result << 1) | (n & 1);
        n >>= 1;
    }

    result
}

/// Replace *arr* with its Fourier Transform using the Cooley-Tukey algorithm
///
/// # Panics
/// Panics if the length of the input array is not a multiple of 2.
/// If your input data does not have the right length, consider padding it with zeros.
///
/// # References
///  - [Fourier Transforms and the Fast Fourier Transform (FFT) Algorithm](https://www.cs.cmu.edu/afs/andrew/scs/cs/15-463/2001/pub/www/notes/fourier/fourier.pdf)
///  - [Wikipedia: Cooley-Tukey FFT algorithm](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm)
///  - [Wikipedia: Butterfly Diagram](https://en.wikipedia.org/wiki/Butterfly_diagram)
pub fn fft(arr: &mut [Complex<f64>]) {
    let n = arr.len();
    let order = (n as f32).log2().round() as u8;

    // bit-reversal
    for j in 0..n {
        let nj = bit_reversed(j as u8, order) as usize;
        if j < nj {
            arr.swap(j, nj);
        }
    }

    // butterfly computations
    let w = Complex::new(0.0, -2.0 * std::f64::consts::PI / n as f64).exp();
    // runs log2(n) times
    for k in 0..order {
        let offset = 2_usize.pow(k.into());
        let step = 2 * offset;
        let multiplier: usize = (n / step).try_into().unwrap();
        // both inner loops together run n/2 times
        for i in (0..n).step_by(step) {
            for j in 0..offset {
                let p = arr[i + j];
                let aq = w.powi((multiplier * j).try_into().unwrap()) * arr[i + j + offset];
                arr[i + j] = p + aq;
                arr[i + j + offset] = p - aq;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use ndarray::Array;
    use num::complex::ComplexFloat;

    use super::*;

    #[test]
    fn test_tf() {
        let pi = std::f32::consts::PI;
        let x = Array::linspace(0.0, 4.0 * pi, 10);
        let y = x.mapv(f32::cos);

        let res = ft(&y);

        for val in res {
            println!("{:?}", val.abs());
        }
    }

    #[test]
    fn test_fft_delta() {
        let mut input = vec![
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ];
        let output = vec![
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
        ];

        fft(&mut input);

        assert_eq!(input, output);
    }

    #[test]
    fn test_fft_constant() {
        let mut input = vec![
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
        ];
        let output = vec![
            Complex::new(8.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ];

        fft(&mut input);

        assert_eq!(input, output);
    }

    #[test]
    #[should_panic]
    fn test_fft_arr_length() {
        let mut input = vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ];

        fft(&mut input);
    }
}
