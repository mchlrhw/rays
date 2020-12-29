use std::fmt;

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
    t: f64,
    front_face: bool,
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
        let target = hit_rec.p + hit_rec.normal + random_vector_in_unit_sphere();
        let recursion = ray_colour(
            &Ray {
                origin: hit_rec.p,
                direction: target - hit_rec.p,
            },
            world,
            depth - 1,
        );

        return 0.5 * recursion;
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

fn main() {
    // World

    let world = HitList(vec![
        Box::new(Sphere {
            center: Point::new(0.0, 0.0, -1.0),
            radius: 0.5,
        }),
        Box::new(Sphere {
            center: Point::new(0.0, -100.5, -1.0),
            radius: 100.0,
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
