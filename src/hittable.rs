mod sphere;

use std::sync::Arc;

use crate::{material::Material, ray::Ray, Point, Vector};

pub use sphere::Sphere;

pub struct HitRecord {
    pub p: Point,
    pub normal: Vector,
    pub material: Arc<dyn Material>,
    pub t: f64,
    pub front_face: bool,
}

pub trait Hittable: Sync {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

pub struct HitList(pub Vec<Box<dyn Hittable>>);

impl Hittable for HitList {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut maybe_hit = None;
        let mut closest_so_far = t_max;

        for hittable in &self.0 {
            if let Some(hit_rec) = hittable.hit(ray, t_min, closest_so_far) {
                closest_so_far = hit_rec.t;
                maybe_hit = Some(hit_rec);
            }
        }

        maybe_hit
    }
}
