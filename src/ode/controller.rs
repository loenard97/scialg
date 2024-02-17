pub struct Controller {
    pub h_next: f64,
    pub err_old: f64,
    pub rejected: bool,
}

impl Controller {
    pub fn new() -> Self {
        Controller {
            h_next: 0.0,
            err_old: 1.0e-4,
            rejected: false,
        }
    }

    pub fn success(&mut self, err: f64, h: f64) -> (bool, f64) {
        let beta = 0.0;
        let alpha = 0.2 - beta * 0.75;
        let safe = 0.9;
        let min_scale = 0.2;
        let max_scale = 10.0;

        if err <= 1.0 {
            let scale = if err == 0.0 {
                max_scale
            } else {
                let temp = safe * err.powf(-alpha) * self.err_old.powf(beta);
                f64::max(min_scale, f64::min(max_scale, temp))
            };
            if self.rejected {
                self.h_next = h * f64::min(scale, 1.0);
            } else {
                self.h_next = h * scale;
            }
            self.err_old = f64::max(err, 1.0e-4);
            self.rejected = false;
            return (true, h);
        }

        // truncation error too large
        let scale = f64::max(safe * err.powf(-alpha), min_scale);
        self.rejected = true;
        
        (false, h * scale)
    }
}
