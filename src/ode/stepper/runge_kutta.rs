use crate::ode::stepper::{Stepper, StepperData};
use crate::vector::Vector;

/// Runge-Kutta method
pub struct RungeKutta<const N: usize> {
    data: StepperData<N>,
}

impl<const N: usize> RungeKutta<N> {
    pub fn new(h: f64, p0: Vector<N>, df: fn(f64, Vector<N>) -> Vector<N>) -> Self {
        RungeKutta {
            data: StepperData::new(h, p0, df),
        }
    }

    fn dy(&self, h: f64) -> Vector<N> {
        let k1 = (self.data.derive)(self.data.x_cur, self.data.y_cur) * h;
        let k2 = (self.data.derive)(self.data.x_cur, self.data.y_cur + k1 / 2.0) * h;
        let k3 = (self.data.derive)(self.data.x_cur, self.data.y_cur + k2 / 2.0) * h;
        let k4 = (self.data.derive)(self.data.x_cur, self.data.y_cur + k3) * h;

        self.data.y_cur + k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0
    }
}

impl<const N: usize> Stepper<N> for RungeKutta<N> {
    fn step(&mut self) -> Vector<N> {
        self.data.y_cur = self.dy(self.data.h_cur);
        self.data.y_cur
    }
}
