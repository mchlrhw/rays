use std::fmt;

use rand::{thread_rng, Rng};

use crate::{Vector, SAMPLES_PER_PIXEL};

#[derive(Clone)]
pub struct Rgb(Vector);

impl Rgb {
    pub fn new(r: f64, g: f64, b: f64) -> Self {
        Self(Vector::new(r, g, b))
    }

    pub fn random_in_range(min: f64, max: f64) -> Self {
        let mut rng = thread_rng();

        Self::new(
            rng.gen_range(min..max),
            rng.gen_range(min..max),
            rng.gen_range(min..max),
        )
    }

    pub fn random() -> Self {
        Self::random_in_range(0.0, 1.0)
    }
}

impl fmt::Display for Rgb {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut r = self.0.x;
        let mut g = self.0.y;
        let mut b = self.0.z;

        let scale = 1.0 / SAMPLES_PER_PIXEL as f64;
        r = (r * scale).sqrt();
        g = (g * scale).sqrt();
        b = (b * scale).sqrt();

        let ir = (256.0 * r.clamp(0.0, 0.999)) as u64;
        let ig = (256.0 * g.clamp(0.0, 0.999)) as u64;
        let ib = (256.0 * b.clamp(0.0, 0.999)) as u64;

        write!(f, "{} {} {}", ir, ig, ib)
    }
}

impl std::ops::Mul<Rgb> for f64 {
    type Output = Rgb;

    fn mul(self, other: Rgb) -> Self::Output {
        Rgb(self * other.0)
    }
}

impl std::ops::Mul<Rgb> for Rgb {
    type Output = Rgb;

    fn mul(self, other: Rgb) -> Self::Output {
        Rgb(self.0.component_mul(&other.0))
    }
}

impl std::ops::Add<Rgb> for Rgb {
    type Output = Rgb;

    fn add(self, other: Rgb) -> Self::Output {
        Rgb(self.0 + other.0)
    }
}

impl std::ops::AddAssign for Rgb {
    fn add_assign(&mut self, other: Self) {
        *self = Self(self.0 + other.0);
    }
}
