//! Ordinary differential equation solver

mod controller;
pub mod stepper;

use crate::ode::stepper::dormand_prince::DormandPrince;
use crate::ode::stepper::euler::Euler;
use crate::ode::stepper::midpoint::Midpoint;
use crate::ode::stepper::runge_kutta::RungeKutta;
use crate::ode::stepper::{Stepper, StepperMethod};
use crate::vector::Vector;

/// Interface for solving ordinary differential equations
pub struct ODESolver<const N: usize> {
    pub steps: usize,
    pub step_size: f64,
    pub cur_step: Vector<N>,
    pub ode: fn(f64, Vector<N>) -> Vector<N>,
    pub stepper: Box<dyn Stepper<N>>,
    pub data: Vec<Vector<N>>,
}

impl<const N: usize> ODESolver<N> {
    pub fn new(
        steps: usize,
        step_size: f64,
        p0: Vector<N>,
        ode: fn(f64, Vector<N>) -> Vector<N>,
        stepper: StepperMethod,
    ) -> Self {
        let stepper_struct: Box<dyn Stepper<N>> = match stepper {
            StepperMethod::Euler => Box::new(Euler::new(step_size, p0, ode)),
            StepperMethod::Midpoint => Box::new(Midpoint::new(step_size, p0, ode)),
            StepperMethod::RungeKutta => Box::new(RungeKutta::new(step_size, p0, ode)),
            StepperMethod::DormandPrince => Box::new(DormandPrince::new(step_size, p0, ode)),
        };

        ODESolver {
            steps,
            step_size,
            cur_step: p0,
            ode,
            stepper: stepper_struct,
            data: Vec::new(),
        }
    }

    pub fn run(&mut self) {
        while self.steps > 0 {
            self.data.push(self.stepper.step());
            self.steps -= 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solver() {
        let steps = 300;
        let step_size = 0.001;
        let p0 = Vector::new(&[0.0, 1.5, 1.0, 1.0]);
        let gravity = |_: f64, x: Vector<4>| Vector::new(&[x[2], x[3], 0.0, -9.81]);
        let stepper_method = StepperMethod::Euler;

        let mut solver = ODESolver::new(steps, step_size, p0, gravity, stepper_method);
        solver.run();

        println!("{:?}", solver.data);
    }
}
