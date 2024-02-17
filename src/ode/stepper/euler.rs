use crate::ode::stepper::{Stepper, StepperData};
use crate::vector::Vector;

/// Euler method
pub struct Euler<const N: usize> {
    data: StepperData<N>,
}

impl<const N: usize> Euler<N> {
    pub fn new(h: f64, p0: Vector<N>, df: fn(f64, Vector<N>) -> Vector<N>) -> Self {
        Euler {
            data: StepperData::new(h, p0, df),
        }
    }

    fn dy(&self, h: f64) -> Vector<N> {
        self.data.y_cur + (self.data.derive)(self.data.x_cur, self.data.y_cur) * h
    }
}

impl<const N: usize> Stepper<N> for Euler<N> {
    fn step(&mut self) -> Vector<N> {
        self.data.y_cur = self.dy(self.data.h_cur);
        self.data.y_cur
    }
}


