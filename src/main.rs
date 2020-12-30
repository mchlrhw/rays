use std::{fmt, sync::Arc};

use indicatif::ProgressBar;
use na::{Point3, Vector3};
use nalgebra as na;
use rand::{thread_rng, Rng};
use rayon::prelude::*;

const ASPECT_RATIO: f64 = 16.0 / 9.0;
const IMAGE_WIDTH: u64 = 400;
const IMAGE_HEIGHT: u64 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u64;
const SAMPLES_PER_PIXEL: u64 = 100;
const MAX_DEPTH: i64 = 50;

type Point = Point3<f64>;
type Vector = Vector3<f64>;

#[derive(Clone)]
struct Rgb(Vector3<f64>);

impl Rgb {
    fn new(r: f64, g: f64, b: f64) -> Self {
        Self(Vector3::new(r, g, b))
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

struct Ray {
    origin: Point,
    direction: Vector,
}

impl Ray {
    fn at(&self, t: f64) -> Point {
        self.origin + (t * self.direction)
    }
}

struct HitRecord {
    p: Point,
    normal: Vector,
    material: Arc<dyn Material>,
    t: f64,
    front_face: bool,
}

trait Material: Sync + Send {
    fn scatter(&self, ray: &Ray, hit_rec: &HitRecord) -> Option<(Rgb, Ray)>;
}

fn near_zero(vec: Vector) -> bool {
    const S: f64 = 1e-8;

    vec.x.abs() < S && vec.y.abs() < S && vec.z.abs() < S
}

struct Lambertian {
    albedo: Rgb,
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

fn reflect(v: Vector, n: Vector) -> Vector {
    v - (2.0 * v.dot(&n) * n)
}

struct Metal {
    albedo: Rgb,
    fuzz: f64,
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

fn refract(uv: Vector, n: Vector, etai_over_etat: f64) -> Vector {
    let cos_theta = (-uv).dot(&n).min(1.0);
    let r_out_perp = etai_over_etat * (uv + (cos_theta * n));
    let r_out_parallel = -(1.0 - r_out_perp.magnitude_squared()).abs().sqrt() * n;

    r_out_perp + r_out_parallel
}

struct Dielectric {
    ir: f64,
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
        let refracted = refract(unit_direction, hit_rec.normal, refraction_ratio);
        let scattered = Ray { origin: hit_rec.p, direction: refracted };

        Some((attenuation, scattered))
    }
}

trait Hittable: Sync {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

struct HitList(Vec<Box<dyn Hittable>>);

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

struct Sphere {
    center: Point,
    radius: f64,
    material: Arc<dyn Material>,
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = ray.origin - self.center;
        let a = ray.direction.magnitude_squared();
        let half_b = oc.dot(&ray.direction);
        let c = oc.magnitude_squared() - self.radius.powi(2);

        let discriminant = half_b.powi(2) - (a * c);
        if discriminant < 0.0 {
            return None;
        }

        let sqrtd = discriminant.sqrt();
        let mut root = (-half_b - sqrtd) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrtd) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }

        let t = root;
        let p = ray.at(t);
        let outward_normal = (p - self.center) / self.radius;
        let front_face = ray.direction.dot(&outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };

        Some(HitRecord {
            t,
            normal,
            material: self.material.clone(),
            p,
            front_face,
        })
    }
}

fn ray_colour(ray: &Ray, world: &dyn Hittable, depth: i64) -> Rgb {
    if depth < 0 {
        return Rgb::new(0.0, 0.0, 0.0);
    }

    if let Some(hit_rec) = world.hit(ray, 0.001, f64::INFINITY) {
        if let Some((attenuation, scattered_ray)) = hit_rec.material.scatter(ray, &hit_rec) {
            return attenuation * ray_colour(&scattered_ray, world, depth - 1);
        }
        return Rgb::new(0.0, 0.0, 0.0);
    }

    let t = 0.5 * (ray.direction.normalize().y + 1.0);

    ((1.0 - t) * Rgb::new(1.0, 1.0, 1.0)) + (t * Rgb::new(0.5, 0.7, 1.0))
}

struct Camera {
    origin: Point,
    lower_left_corner: Point,
    horizontal: Vector,
    vertical: Vector,
}

impl Camera {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const VIEWPORT_HEIGHT: f64 = 2.0;
    const VIEWPORT_WIDTH: f64 = Self::ASPECT_RATIO * Self::VIEWPORT_HEIGHT;
    const FOCAL_LENGTH: f64 = 1.0;

    fn new() -> Self {
        let origin = Point::new(0.0, 0.0, 0.0);
        let horizontal = Vector::new(Self::VIEWPORT_WIDTH, 0.0, 0.0);
        let vertical = Vector::new(0.0, Self::VIEWPORT_HEIGHT, 0.0);
        let lower_left_corner = origin
            - (horizontal / 2.0)
            - (vertical / 2.0)
            - Vector::new(0.0, 0.0, Self::FOCAL_LENGTH);

        Self {
            origin,
            horizontal,
            vertical,
            lower_left_corner,
        }
    }

    fn get_ray(&self, u: f64, v: f64) -> Ray {
        Ray {
            origin: self.origin,
            direction: self.lower_left_corner + (u * self.horizontal) + (v * self.vertical)
                - self.origin,
        }
    }
}

fn random_vector(min: f64, max: f64) -> Vector {
    let mut rng = thread_rng();

    Vector::new(
        rng.gen_range(min..max),
        rng.gen_range(min..max),
        rng.gen_range(min..max),
    )
}

fn random_vector_in_unit_sphere() -> Vector {
    loop {
        let p = random_vector(-1.0, 1.0);
        if p.magnitude_squared() > 1.0 {
            continue;
        }

        return p;
    }
}

fn random_unit_vector() -> Vector {
    random_vector_in_unit_sphere().normalize()
}

fn main() {
    // World

    let material_ground = Arc::new(Lambertian {
        albedo: Rgb::new(0.8, 0.8, 0.0),
    });
    let material_center = Arc::new(Dielectric {
        ir: 1.5,
    });
    let material_left = Arc::new(Dielectric {
        ir: 1.5,
    });
    let material_right = Arc::new(Metal {
        albedo: Rgb::new(0.8, 0.6, 0.2),
        fuzz: 1.0,
    });

    let world = HitList(vec![
        Box::new(Sphere {
            center: Point::new(0.0, -100.5, -1.0),
            radius: 100.0,
            material: material_ground,
        }),
        Box::new(Sphere {
            center: Point::new(0.0, 0.0, -1.0),
            radius: 0.5,
            material: material_center,
        }),
        Box::new(Sphere {
            center: Point::new(-1.0, 0.0, -1.0),
            radius: 0.5,
            material: material_left,
        }),
        Box::new(Sphere {
            center: Point::new(1.0, 0.0, -1.0),
            radius: 0.5,
            material: material_right,
        }),
    ]);

    // Camera

    let camera = Camera::new();

    // Render

    println!("P3\n{} {}\n255", IMAGE_WIDTH, IMAGE_HEIGHT);

    let progress = ProgressBar::new(IMAGE_HEIGHT);
    for j in (0..IMAGE_HEIGHT).rev() {
        progress.inc(1);
        for i in 0..IMAGE_WIDTH {
            let pixel_colour = (0..SAMPLES_PER_PIXEL)
                .into_par_iter()
                .map(|_s| {
                    let mut rng = thread_rng();

                    let u = (i as f64 + rng.gen_range(0.0..1.0)) / (IMAGE_WIDTH - 1) as f64;
                    let v = (j as f64 + rng.gen_range(0.0..1.0)) / (IMAGE_HEIGHT - 1) as f64;

                    let ray = camera.get_ray(u, v);

                    ray_colour(&ray, &world, MAX_DEPTH)
                })
                .reduce(|| Rgb::new(0.0, 0.0, 0.0), |pc, rc| pc + rc);

            println!("{}", pixel_colour);
        }
    }
    progress.finish();
}
