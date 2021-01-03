use crate::{random_unit_vector, ray::Ray, rgb::Rgb};

use super::{reflect, HitRecord, Material};

pub struct Metal {
    pub albedo: Rgb,
    pub fuzz: f64,
}

impl Material for Metal {
    fn scatter(&self, ray: &Ray, hit_rec: &HitRecord) -> Option<(Rgb, Ray)> {
        let reflected = reflect(ray.direction.normalize(), hit_rec.normal);

        let scattered = Ray {
            origin: hit_rec.p,
            direction: reflected + (self.fuzz * random_unit_vector()),
        };
        let attenuation = self.albedo.clone();

        if scattered.direction.dot(&hit_rec.normal) > 0.0 {
            Some((attenuation, scattered))
        } else {
            None
        }
    }
}
