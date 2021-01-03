use rand::{thread_rng, Rng};

use crate::{ray::Ray, rgb::Rgb};

use super::{reflect, reflectance, refract, HitRecord, Material};

pub struct Dielectric {
    pub ir: f64,
}

impl Material for Dielectric {
    fn scatter(&self, ray: &Ray, hit_rec: &HitRecord) -> Option<(Rgb, Ray)> {
        let attenuation = Rgb::new(1.0, 1.0, 1.0);
        let refraction_ratio = if hit_rec.front_face {
            1.0 / self.ir
        } else {
            self.ir
        };

        let unit_direction = ray.direction.normalize();
        let cos_theta = (-unit_direction).dot(&hit_rec.normal).min(1.0);
        let sin_theta = (1.0 - cos_theta.powi(2)).sqrt();

        let cannot_refract = refraction_ratio * sin_theta > 1.0
            || reflectance(cos_theta, refraction_ratio) > thread_rng().gen_range(0.0..1.0);

        let direction = if cannot_refract {
            reflect(unit_direction, hit_rec.normal)
        } else {
            refract(unit_direction, hit_rec.normal, refraction_ratio)
        };

        let scattered = Ray {
            origin: hit_rec.p,
            direction,
        };

        Some((attenuation, scattered))
    }
}
