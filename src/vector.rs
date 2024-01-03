#[derive(Debug)]
pub enum Axis {
    X,
    Y,
    Z,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3D {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vector3D { x, y, z }
    }

    pub fn from_axis(axis: Axis) -> Self {
        match axis {
            Axis::X => Vector3D {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            Axis::Y => Vector3D {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
            Axis::Z => Vector3D {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
        }
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

    /// Return a new vector rotated around *axis* by angle *theta* (in radians)
    ///
    /// # Example
    /// ```
    /// use scialg::vector::{Axis, Vector3D};
    ///
    /// let v = Vector3D::from_axis(Axis::X);
    /// let u = Vector3D::from_axis(Axis::Z);
    /// let theta = 0.5 * std::f64::consts::PI;
    ///
    /// let r = v.rotate(u, theta);
    ///
    /// assert!(r.near_equal(&Vector3D::from_axis(Axis::Y), 1e-7));
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
}

/*
    def scalar_product(self, other) -> float:
        """
        Scalar Product between self and other
        """
        return Vector3D(np.dot(self.coords, other.coords))

    def cross_product(self, other) -> "Vector3D":
        """
        Cross Product between self and other
        """
        return Vector3D(np.cross(self.coords, other.coords))

    def angle(self, other) -> float:
        """
        Angle between self and other in radians
        """
        v1 = self.normalize()
        v2 = other.normalize()
        return np.arccos(np.clip(np.dot(v1.coords, v2.coords), -1.0, 1.0))

    def project(self, plane) -> "Vector3D":
        """
        Project self onto plane
        """
        if plane == 'xy':
            return Vector3D(self.x, self.y, 0)

        if plane == 'xz':
            return Vector3D(self.x, 0, self.z)

        if plane == 'yz':
            return Vector3D(0, self.y, self.z)

        raise ValueError('Plane has to be xy, xz or yz')

*/
