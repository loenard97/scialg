mod stepper;

use crate::vector::Vector;
use stepper::*;


pub struct ODESolver<const N: usize> {
    pub steps: usize,
    pub step_size: f64,
    pub cur_step: Vector<N>,
    pub ode: fn(Vector<N>) -> Vector<N>,
    pub stepper: StepperMethod,
}

impl<const N: usize> ODESolver<N> {
    pub fn new(
        steps: usize,
        step_size: f64,
        p0: Vector<N>,
        ode: fn(Vector<N>) -> Vector<N>,
        stepper: StepperMethod,
    ) -> Self {
        ODESolver {
            steps,
            step_size,
            cur_step: p0,
            ode,
            stepper,
        }
    }
}

impl<const N: usize> Iterator for ODESolver<N> {
    type Item = Vector<N>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.steps == 0 {
            return None;
        }

        let stepper_fn = match self.stepper {
            StepperMethod::Euler => euler::<N>,
            StepperMethod::Midpoint => midpoint::<N>,
            StepperMethod::RungeKutta => runge_kutta::<N>,
        };
        self.cur_step = stepper_fn(self.ode, self.cur_step, self.step_size);
        self.steps -= 1;

        Some(self.cur_step)
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
        let gravity = |x: Vector<4>| Vector::new(&[x.coeff[2], x.coeff[3], 0.0, -9.81]);
        let stepper_method = StepperMethod::RungeKutta;

        let solver = ODESolver::new(steps, step_size, p0, gravity, stepper_method);

        for x in solver {
            println!("{:?}", x);
        }
    }
}
