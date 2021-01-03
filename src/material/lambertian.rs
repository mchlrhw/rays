use crate::{random_unit_vector, ray::Ray, rgb::Rgb, Vector};

use super::{HitRecord, Material};

fn near_zero(vec: Vector) -> bool {
    const S: f64 = 1e-8;

    vec.x.abs() < S && vec.y.abs() < S && vec.z.abs() < S
}

pub struct Lambertian {
    pub albedo: Rgb,
}

impl Material for Lambertian {
    fn scatter(&self, _ray: &Ray, hit_rec: &HitRecord) -> Option<(Rgb, Ray)> {
        let mut direction = hit_rec.normal + random_unit_vector();
        if near_zero(direction) {
            direction = hit_rec.normal;
        }

        let scattered = Ray {
            origin: hit_rec.p,
            direction,
        };
        let attenuation = self.albedo.clone();

        Some((attenuation, scattered))
    }
}
