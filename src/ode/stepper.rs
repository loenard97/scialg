use crate::vector::Vector;

pub enum StepperMethod {
    Euler,
    Midpoint,
    RungeKutta,
}

/// Euler method
pub fn euler<const N: usize>(f: fn(Vector<N>) -> Vector<N>, p0: Vector<N>, h: f64) -> Vector<N> {
    p0 + f(p0) * h
}

/// Midpoint method
pub fn midpoint<const N: usize>(f: fn(Vector<N>) -> Vector<N>, p0: Vector<N>, h: f64) -> Vector<N> {
    let k1 = f(p0) * h;
    let k2 = f(p0 + k1 / 2.0) * h;

    p0 + k2
}

/// Runge-Kutta method
pub fn runge_kutta<const N: usize>(f: fn(Vector<N>) -> Vector<N>, p0: Vector<N>, h: f64) -> Vector<N> {
    let k1 = f(p0) * h;
    let k2 = f(p0 + k1 / 2.0) * h;
    let k3 = f(p0 + k2 / 2.0) * h;
    let k4 = f(p0 + k3) * h;

    p0 + k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0
}
