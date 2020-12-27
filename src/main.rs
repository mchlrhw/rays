use std::fmt;

use indicatif::ProgressBar;
use na::{Point3, Vector3};
use nalgebra as na;

const ASPECT_RATIO: f64 = 16.0 / 9.0;
const IMAGE_WIDTH: u64 = 400;
const IMAGE_HEIGHT: u64 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u64;

const VIEWPORT_HEIGHT: f64 = 2.0;
const VIEWPORT_WIDTH: f64 = ASPECT_RATIO * VIEWPORT_HEIGHT;
const FOCAL_LENGTH: f64 = 1.0;

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
        let ir = (255.999 * self.0.x) as u64;
        let ig = (255.999 * self.0.y) as u64;
        let ib = (255.999 * self.0.z) as u64;

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

struct Ray {
    origin: Point,
    direction: Vector,
}

impl Ray {
    fn at(&self, t: f64) -> Point {
        self.origin + (t * self.direction)
    }
}

fn hit_sphere(center: Point, radius: f64, ray: &Ray) -> f64 {
    let oc = ray.origin - center;

    let a = ray.direction.dot(&ray.direction);
    let b = 2.0 * oc.dot(&ray.direction);
    let c = oc.dot(&oc) - (radius * radius);

    let discriminant = (b * b) - (4.0 * a * c);

    if discriminant < 0.0 {
        -1.0
    } else {
        (-b - discriminant.sqrt()) / (2.0 * a)
    }
}

fn ray_colour(ray: &Ray) -> Rgb {
    let t = hit_sphere(Point::new(0.0, 0.0, -1.0), 0.5, ray);
    if t > 0.0 {
        let n = (ray.at(t) - Vector::new(0.0, 0.0, -1.0)).coords.normalize();

        return 0.5 * Rgb::new(n.x + 1.0, n.y + 1.0, n.z + 1.0);
    }

    let t = 0.5 * (ray.direction.normalize().y + 1.0);

    ((1.0 - t) * Rgb::new(1.0, 1.0, 1.0)) + (t * Rgb::new(0.5, 0.7, 1.0))
}

fn main() {
    // Camera

    let origin = Point::new(0.0, 0.0, 0.0);
    let horizontal = Vector::new(VIEWPORT_WIDTH, 0.0, 0.0);
    let vertical = Vector::new(0.0, VIEWPORT_HEIGHT, 0.0);
    let lower_left_corner =
        origin - (horizontal / 2.0) - (vertical / 2.0) - Vector::new(0.0, 0.0, FOCAL_LENGTH);

    // Render

    println!("P3\n{} {}\n255", IMAGE_WIDTH, IMAGE_HEIGHT);

    let progress = ProgressBar::new(IMAGE_HEIGHT);
    for j in (0..IMAGE_HEIGHT).rev() {
        progress.inc(1);
        for i in 0..IMAGE_WIDTH {
            let u = i as f64 / (IMAGE_WIDTH - 1) as f64;
            let v = j as f64 / (IMAGE_HEIGHT - 1) as f64;

            let ray = Ray {
                origin,
                direction: lower_left_corner + (u * horizontal) + (v * vertical) - origin,
            };
            let pixel_colour = ray_colour(&ray);

            println!("{}", pixel_colour);
        }
    }
    progress.finish();
}
