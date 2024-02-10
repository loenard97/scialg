use num::Zero;

use crate::ode::controller::Controller;
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

pub trait Stepper {
    fn step(&mut self);
}

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

impl<const N: usize> Stepper for Euler<N> {
    fn step(&mut self) {
        self.data.y_cur = self.dy(self.data.h_cur);
    }
}

/// Midpoint method
pub struct Midpoint<const N: usize> {
    data: StepperData<N>,
}

impl<const N: usize> Midpoint<N> {
    pub fn new(h: f64, p0: Vector<N>, df: fn(f64, Vector<N>) -> Vector<N>) -> Self {
        Midpoint {
            data: StepperData::new(h, p0, df),
        }
    }

    fn dy(&self, h: f64) -> Vector<N> {
        let k1 = (self.data.derive)(self.data.x_cur, self.data.y_cur) * h;
        let k2 = (self.data.derive)(self.data.x_cur, self.data.y_cur + k1 / 2.0) * h;

        self.data.y_cur + k2
    }
}

impl<const N: usize> Stepper for Midpoint<N> {
    fn step(&mut self) {
        self.data.y_cur = self.dy(self.data.h_cur);
    }
}

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

impl<const N: usize> Stepper for RungeKutta<N> {
    fn step(&mut self) {
        self.data.y_cur = self.dy(self.data.h_cur);
    }
}

/// Dormand Prince method
///
/// # References
///  - [Wikipedia: Dormand-Prince
///  method](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method)
pub struct DormandPrince<const N: usize> {
    x: f64,
    x_old: f64,
    y: Vector<N>,
    dydx: Vector<N>,
    h_did: f64,
    h_next: f64,
    y_out: Vector<N>,
    y_err: Vector<N>,
    k2: Vector<N>,
    k3: Vector<N>,
    k4: Vector<N>,
    k5: Vector<N>,
    k6: Vector<N>,
    rcont1: Vector<N>,
    rcont2: Vector<N>,
    rcont3: Vector<N>,
    rcont4: Vector<N>,
    rcont5: Vector<N>,
    dydxnew: Vector<N>,
    atol: f64,
    rtol: f64,
    controller: Controller,
    data: StepperData<N>,
}

impl<const N: usize> DormandPrince<N> {
    pub fn new(h: f64, p0: Vector<N>, df: fn(f64, Vector<N>) -> Vector<N>) -> Self {
        let x = 0.0;
        let x_old = 0.0;
        let y = p0;
        let dydx = Vector::zero();
        let y_out = Vector::zero();
        let y_err = Vector::zero();
        let h_did = 0.0;
        let h_next = 0.0;
        let k2 = Vector::zero();
        let k3 = Vector::zero();
        let k4 = Vector::zero();
        let k5 = Vector::zero();
        let k6 = Vector::zero();
        let rcont1 = Vector::zero();
        let rcont2 = Vector::zero();
        let rcont3 = Vector::zero();
        let rcont4 = Vector::zero();
        let rcont5 = Vector::zero();
        let dydxnew = Vector::zero();
        let atol = 0.01;
        let rtol = 0.01;
        let controller = Controller::new();

        let data = StepperData::new(h, p0, df);

        DormandPrince {
            x,
            x_old,
            y,
            dydx,
            y_out,
            y_err,
            h_did,
            h_next,
            k2,
            k3,
            k4,
            k5,
            k6,
            rcont1,
            rcont2,
            rcont3,
            rcont4,
            rcont5,
            dydxnew,
            atol,
            rtol,
            controller,
            data,
        }
    }

    fn error(&mut self) -> f64 {
        let mut err = 0.0;

        for i in 0..self.y.dim() {
            let sk = self.atol + self.rtol * f64::max(self.y[i].abs(), self.y_out[i].abs());
            err += (self.y_err[i] / sk).powi(2);
        }
        (err / self.y.dim() as f64).sqrt()
    }

    fn dy(&mut self) {
        // Dormand-Prince 5(4) parameters
        static C2: f64 = 0.2;
        static C3: f64 = 0.3;
        static C4: f64 = 0.8;
        static C5: f64 = 8.0 / 9.0;

        static A21: f64 = 0.2;
        static A31: f64 = 3.0 / 40.0;
        static A32: f64 = 9.0 / 40.0;
        static A41: f64 = 44.0 / 45.0;
        static A42: f64 = -56.0 / 15.0;
        static A43: f64 = 32.0 / 9.0;
        static A51: f64 = 19372.0 / 6561.0;
        static A52: f64 = -25360.0 / 2187.0;
        static A53: f64 = 64448.0 / 6561.0;
        static A54: f64 = -212.0 / 729.0;
        static A61: f64 = 9017.0 / 3168.0;
        static A62: f64 = -355.0 / 33.0;
        static A63: f64 = 46732.0 / 5247.0;
        static A64: f64 = 49.0 / 176.0;
        static A65: f64 = -5103.0 / 18656.0;
        static A71: f64 = 35.0 / 384.0;
        static A72: f64 = 0.0;
        static A73: f64 = 500.0 / 1113.0;
        static A74: f64 = 125.0 / 192.0;
        static A75: f64 = -2187.0 / 6784.0;
        static A76: f64 = 11.0 / 84.0;

        static E1: f64 = 71.0 / 57600.0;
        static E3: f64 = -71.0 / 16695.0;
        static E4: f64 = 71.0 / 1920.0;
        static E5: f64 = -17253.0 / 339200.0;
        static E6: f64 = 22.0 / 525.0;
        static E7: f64 = -1.0 / 40.0;

        let mut y_temp = self.y + self.dydx * A21 * self.data.h_cur;
        self.k2 = (self.data.derive)(self.x + C2 * self.data.h_cur, y_temp);

        y_temp = self.y + (self.dydx * A31 + self.k2 * A32) * self.data.h_cur;
        self.k3 = (self.data.derive)(self.x + C3 * self.data.h_cur, y_temp);

        y_temp = self.y + (self.dydx * A41 + self.k2 * A42 + self.k3 * A43) * self.data.h_cur;
        self.k4 = (self.data.derive)(self.x + C4 * self.data.h_cur, y_temp);

        y_temp = self.y
            + (self.dydx * A51 + self.k2 * A52 + self.k3 * A53 + self.k4 * A54) * self.data.h_cur;
        self.k5 = (self.data.derive)(self.x + C5 * self.data.h_cur, y_temp);

        y_temp = self.y
            + (self.dydx * A61 + self.k2 * A62 + self.k3 * A63 + self.k4 * A64 + self.k5 * A65)
                * self.data.h_cur;
        let xph = self.x + self.data.h_cur;
        self.k6 = (self.data.derive)(xph, y_temp);

        self.y_out = self.y
            + (self.dydx * A71
                + self.k2 * A72
                + self.k3 * A73
                + self.k4 * A74
                + self.k5 * A75
                + self.k6 * A76)
                * self.data.h_cur;
        self.dydxnew = (self.data.derive)(xph, self.y_out);

        self.y_err = (self.dydx * E1
            + self.k3 * E3
            + self.k4 * E4
            + self.k5 * E5
            + self.k6 * E6
            + self.dydxnew * E7)
            * self.data.h_cur;
    }
}

impl<const N: usize> Stepper for DormandPrince<N> {
    fn step(&mut self) {
        loop {
            self.dy();
            let err = self.error();
            let (success, h_new) = self.controller.success(err, self.data.h_cur);
            if success {
                break;
            }
            self.data.h_cur = h_new;
        }

        self.dydx = self.dydxnew;
        self.y = self.y_out;
        self.x_old = self.x;
        self.h_did = self.data.h_cur;
        self.x += self.h_did;
        self.h_next = self.controller.h_next;
    }
}
