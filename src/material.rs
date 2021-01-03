mod dielectric;
mod lambertian;
mod metal;



use crate::{hittable::HitRecord, ray::Ray, rgb::Rgb, Vector};

pub use dielectric::Dielectric;
pub use lambertian::Lambertian;
pub use metal::Metal;

pub trait Material: Sync + Send {
    fn scatter(&self, ray: &Ray, hit_rec: &HitRecord) -> Option<(Rgb, Ray)>;
}

fn reflect(v: Vector, n: Vector) -> Vector {
    v - (2.0 * v.dot(&n) * n)
}

fn refract(uv: Vector, n: Vector, etai_over_etat: f64) -> Vector {
    let cos_theta = (-uv).dot(&n).min(1.0);
    let r_out_perp = etai_over_etat * (uv + (cos_theta * n));
    let r_out_parallel = -(1.0 - r_out_perp.magnitude_squared()).abs().sqrt() * n;

    r_out_perp + r_out_parallel
}

fn reflectance(cosine: f64, ref_idx: f64) -> f64 {
    let r0 = ((1.0 - ref_idx) / (1.0 + ref_idx)).powi(2);

    r0 + ((1.0 - r0) * (1.0 - cosine).powi(5))
}
