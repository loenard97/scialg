#![allow(dead_code)]

use ndarray::{Array, Array1};
use num::Float;

#[derive(Debug, PartialEq)]
pub struct NevilleCoeff<F: Float> {
    coeff: Vec<Array1<F>>
}

pub fn neville<F: Float + std::fmt::Display>(xi: Array1<F>, yi: Array1<F>, x: F) -> NevilleCoeff<F> {
    let mut coeff: NevilleCoeff<F> = NevilleCoeff { coeff: Vec::new() };
    coeff.coeff.push(yi.clone());

    for i in 1..xi.len() {
        let pre_p = coeff.coeff[i-1].clone();
        let mut cur_p: Array1<F> = Array::zeros(coeff.coeff[i-1].len() - 1);
        for j in 0..cur_p.len() {
            cur_p[j] = ((x - xi[j+i]) * pre_p[j] + (xi[j] - x) * pre_p[j+1]) / (xi[j] - xi[j+i])
        }
        coeff.coeff.push(cur_p.clone());
    }

    coeff
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_neville() {
        let x: Array1<f32> = Array::from_vec(vec![-1.0, 0.0, 2.0]);
        let y: Array1<f32> = Array::from_vec(vec![ 2.0, 0.0, 4.0]);

        let output = neville(x, y, 1.0);
        println!("{:?}", output.coeff.last().unwrap()[0]);

        assert_eq!(output.coeff.last().unwrap()[0], 2.0);
    }
}
