use std::ops::{Add, Div, Mul, Sub};

use num::Zero;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vector<const N: usize> {
    pub coeff: [f64; N],
}

impl<const N: usize> Vector<N> {
    pub fn dim() -> usize {
        N
    }

    pub fn new(coeff: &[f64]) -> Self {
        if coeff.len() != N {
            panic!("length of input does not match dimensionality of vector");
        }

        let mut arr = [0.0; N];
        for i in 0..N {
            arr[i] = coeff[i];
        }

        Vector { coeff: arr }
    }

    /// Return the vector length
    ///
    /// # Example
    /// ```
    /// use scialg::vector::Vector;
    ///
    /// let v: Vector<3> = Vector::new(&[1.0, -1.0, 1.0]);
    ///
    /// assert_eq!(v.length(), 3.0_f64.sqrt());
    /// ```
    pub fn length(self) -> f64 {
        self.coeff.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    /// Return a normalized vector from self
    ///
    /// # Example
    /// ```
    /// use scialg::vector::Vector;
    ///
    /// let mut v: Vector<3> = Vector::new(&[3.0, -1.0, 2.0]);
    /// v = v.normalize();
    ///
    /// assert_eq!(v.length(), 1.0);
    /// ```
    pub fn normalize(self) -> Self {
        self / self.length()
    }

    /// Return scalar product between *self* and *rhs*
    pub fn scalar_product(self, rhs: &Self) -> f64 {
        let mut sum = 0.0;

        for i in 0..N {
            sum += self.coeff[i] * rhs.coeff[i];
        }

        sum
    }

    /// Return angle between *self* and *rhs* in radians
    ///
    /// # Example
    /// ```
    /// use scialg::vector::Vector;
    ///
    /// let pi_half = std::f64::consts::FRAC_PI_2;
    /// let x: Vector<3> = Vector::new(&[1.0, 0.0, 0.0]);
    /// let y: Vector<3> = Vector::new(&[0.0, 1.0, 0.0]);
    /// let angle = x.angle(&y);
    ///
    /// assert_eq!(angle, pi_half);
    /// ```
    pub fn angle(self, other: &Self) -> f64 {
        let v1 = self.normalize();
        let v2 = other.normalize();

        v1.scalar_product(&v2).acos()
    }

    /// Return a new vector rotated around *axis* by angle *theta* (in radians)
    ///
    /// # Example
    /// ```
    /// use scialg::vector::Vector;
    ///
    /// let x: Vector<3> = Vector::new(&[1.0, 0.0, 0.0]);
    /// let z: Vector<3> = Vector::new(&[0.0, 0.0, 1.0]);
    /// let theta = 0.5 * std::f64::consts::PI;
    ///
    /// // let r = x.rotate(z, theta);
    ///
    /// // assert!((r - Vector::new(&[0.0, 1.0, 0.0])).length() < 1e-5);
    /// ```
    pub fn rotate(self, _axis: Vector<N>, _theta: f64) -> Self {
        assert_eq!(N, 3);
        todo!()

        // let u = axis.normalize();
        // let t_cos = theta.cos();
        // let t_sin = theta.sin();

        // Vector<N> {
        //     x: (t_cos + u.x.powi(2) * (1.0 - t_cos)) * self.x
        //         + (u.x * u.y * (1.0 - t_cos) - u.z * t_sin) * self.y
        //         + (u.x * u.z * (1.0 - t_cos) + u.y * t_sin) * self.z,
        //     y: (u.y * u.x * (1.0 - t_cos) + u.z * t_sin) * self.x
        //         + (t_cos + u.y.powi(2) * (1.0 - t_cos)) * self.y
        //         + (u.y * u.z * (1.0 - t_cos) - u.x * t_sin) * self.z,
        //     z: (u.z * u.x * (1.0 - t_cos) - u.y * t_sin) * self.x
        //         + (u.z * u.y * (1.0 - t_cos) + u.x * t_sin) * self.y
        //         + (t_cos + u.z.powi(2) * (1.0 - t_cos)) * self.z,
        // }
    }

    /// Return cross product between *self* and *other*
    pub fn cross_product(self, _other: &Self) -> Self {
        assert_eq!(N, 3);
        todo!()

        // Vector<N> {
        //     x: self.y * other.z - self.z * other.y,
        //     y: self.z * other.x - self.x * other.z,
        //     z: self.x * other.y - self.y * other.x,
        // }
    }
}

impl<const N: usize> Zero for Vector<N> {
    fn zero() -> Self {
        Vector { coeff: [0.0; N] }
    }

    fn is_zero(&self) -> bool {
        self.coeff.iter().all(|x| *x == 0.0)
    }
}

impl<const N: usize> Add for Vector<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut cs = [0.0; N];
        for i in 0..N {
            cs[i] = self.coeff[i] + rhs.coeff[i];
        }

        Vector { coeff: cs }
    }
}

impl<const N: usize> Sub for Vector<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut cs = [0.0; N];
        for i in 0..N {
            cs[i] = self.coeff[i] - rhs.coeff[i];
        }

        Vector { coeff: cs }
    }
}

impl<const N: usize> Mul<f64> for Vector<N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        let mut cs = [0.0; N];
        for i in 0..N {
            cs[i] = self.coeff[i] * rhs;
        }

        Vector { coeff: cs }
    }
}

impl<const N: usize> Div<f64> for Vector<N> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        let mut cs = [0.0; N];
        for i in 0..N {
            cs[i] = self.coeff[i] / rhs;
        }

        Vector { coeff: cs }
    }
}

