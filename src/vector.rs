#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3D {
    pub const ZERO: Vector3D = Vector3D {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    pub const X: Vector3D = Vector3D {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    pub const Y: Vector3D = Vector3D {
        x: 0.0,
        y: 1.0,
        z: 0.0,
    };
    pub const Z: Vector3D = Vector3D {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };

    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vector3D { x, y, z }
    }

    pub fn near_equal(self, other: &Vector3D, tolerance: f64) -> bool {
        (self.x - other.x).abs() < tolerance
            && (self.y - other.y).abs() < tolerance
            && (self.z - other.z).abs() < tolerance
    }

    /// Return the vector length
    ///
    /// # Example
    /// ```
    /// use scialg::vector::Vector3D;
    ///
    /// let v = Vector3D::new(1.0, -1.0, 1.0);
    ///
    /// assert_eq!(v.length(), 3.0_f64.sqrt());
    /// ```
    pub fn length(self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    /// Return a normalized vector from self
    ///
    /// # Example
    /// ```
    /// use scialg::vector::Vector3D;
    ///
    /// let mut v = Vector3D::new(3.0, -1.0, 2.0);
    /// v = v.normalize();
    ///
    /// assert_eq!(v.length(), 1.0);
    /// ```
    pub fn normalize(self) -> Self {
        let len = self.length();
        Vector3D {
            x: self.x / len,
            y: self.y / len,
            z: self.z / len,
        }
    }

    /// Return angle between *self* and *other* in radians
    ///
    /// # Example
    /// ```
    /// use scialg::vector::Vector3D;
    ///
    /// let pi_half = std::f64::consts::FRAC_PI_2;
    /// let angle = Vector3D::X.angle(&Vector3D::Y);
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
    /// use scialg::vector::Vector3D;
    ///
    /// let v = Vector3D::X;
    /// let u = Vector3D::Z;
    /// let theta = 0.5 * std::f64::consts::PI;
    ///
    /// let r = v.rotate(u, theta);
    ///
    /// assert!(r.near_equal(&Vector3D::Y, 1e-7));
    /// ```
    pub fn rotate(self, axis: Vector3D, theta: f64) -> Self {
        let u = axis.normalize();
        let t_cos = theta.cos();
        let t_sin = theta.sin();

        Vector3D {
            x: (t_cos + u.x.powi(2) * (1.0 - t_cos)) * self.x
                + (u.x * u.y * (1.0 - t_cos) - u.z * t_sin) * self.y
                + (u.x * u.z * (1.0 - t_cos) + u.y * t_sin) * self.z,
            y: (u.y * u.x * (1.0 - t_cos) + u.z * t_sin) * self.x
                + (t_cos + u.y.powi(2) * (1.0 - t_cos)) * self.y
                + (u.y * u.z * (1.0 - t_cos) - u.x * t_sin) * self.z,
            z: (u.z * u.x * (1.0 - t_cos) - u.y * t_sin) * self.x
                + (u.z * u.y * (1.0 - t_cos) + u.x * t_sin) * self.y
                + (t_cos + u.z.powi(2) * (1.0 - t_cos)) * self.z,
        }
    }

    /// Return scalar product between *self* and *other*
    pub fn scalar_product(self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Return cross product between *self* and *other*
    pub fn cross_product(self, other: &Self) -> Self {
        Vector3D {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
}
