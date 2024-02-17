pub mod dormand_prince;
pub mod euler;
pub mod midpoint;
pub mod runge_kutta;

use crate::vector::Vector;

pub enum StepperMethod {
    Euler,
    Midpoint,
    RungeKutta,
    DormandPrince,
}

struct StepperData<const N: usize> {
    h_cur: f64,
    h_new: f64,
    x_cur: f64,
    x_old: f64,
    y_cur: Vector<N>,
    y_out: Vector<N>,
    derive: fn(f64, Vector<N>) -> Vector<N>,
}

impl<const N: usize> StepperData<N> {
    fn new(h: f64, p0: Vector<N>, df: fn(f64, Vector<N>) -> Vector<N>) -> Self {
        StepperData {
            h_cur: h,
            h_new: h,
            x_cur: 0.0,
            x_old: 0.0,
            y_cur: p0,
            y_out: p0,
            derive: df,
        }
    }
}

pub trait Stepper<const N: usize> {
    fn step(&mut self) -> Vector<N>;
}

