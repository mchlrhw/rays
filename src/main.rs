mod camera;
mod hittable;
mod material;
mod ray;
mod rgb;

use std::sync::Arc;

use indicatif::ProgressBar;
use na::{Point3, Vector3};
use nalgebra as na;
use rand::{thread_rng, Rng};
use rayon::prelude::*;

use crate::{
    camera::Camera,
    hittable::{HitList, Hittable, Sphere},
    material::{Dielectric, Lambertian, Material, Metal},
    ray::Ray,
    rgb::Rgb,
};

const ASPECT_RATIO: f64 = 3.0 / 2.0;
const IMAGE_WIDTH: u64 = 1200;
const IMAGE_HEIGHT: u64 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u64;
const SAMPLES_PER_PIXEL: u64 = 500;
const MAX_DEPTH: i64 = 50;

type Point = Point3<f64>;
type Vector = Vector3<f64>;

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

    let mut hit_list: Vec<Box<dyn Hittable>> = vec![];

    let material_ground = Arc::new(Lambertian {
        albedo: Rgb::new(0.5, 0.5, 0.5),
    });
    hit_list.push(Box::new(Sphere {
        center: Point::new(0.0, -1000.0, 0.0),
        radius: 1000.0,
        material: material_ground,
    }));

    for a in -11..11 {
        for b in -11..11 {
            let a = a as f64;
            let b = b as f64;

            let mut rng = thread_rng();

            let choose_mat = rng.gen_range(0.0..1.0);
            let center = Point::new(
                a + (0.9 * rng.gen_range(0.0..1.0)),
                0.2,
                b + (0.9 * rng.gen_range(0.0..1.0)),
            );

            if (center - Point::new(4.0, 0.2, 0.0)).magnitude() > 0.9 {
                let material: Arc<dyn Material> = if choose_mat < 0.8 {
                    // diffuse
                    let albedo = Rgb::random() * Rgb::random();
                    Arc::new(Lambertian { albedo })
                } else if choose_mat < 0.95 {
                    // metal
                    let albedo = Rgb::random_in_range(0.5, 1.0);
                    let fuzz = rng.gen_range(0.0..0.5);
                    Arc::new(Metal { albedo, fuzz })
                } else {
                    // glass
                    Arc::new(Dielectric { ir: 1.5 })
                };

                hit_list.push(Box::new(Sphere {
                    center,
                    radius: 0.2,
                    material,
                }));
            }
        }
    }

    let material1 = Arc::new(Dielectric { ir: 1.5 });
    hit_list.push(Box::new(Sphere {
        center: Point::new(0.0, 1.0, 0.0),
        radius: 1.0,
        material: material1,
    }));

    let material2 = Arc::new(Lambertian {
        albedo: Rgb::new(0.4, 0.2, 0.1),
    });
    hit_list.push(Box::new(Sphere {
        center: Point::new(-4.0, 1.0, 0.0),
        radius: 1.0,
        material: material2,
    }));

    let material3 = Arc::new(Metal {
        albedo: Rgb::new(0.7, 0.6, 0.5),
        fuzz: 0.0,
    });
    hit_list.push(Box::new(Sphere {
        center: Point::new(4.0, 1.0, 0.0),
        radius: 1.0,
        material: material3,
    }));

    let world = HitList(hit_list);

    // Camera

    let lookfrom = Point::new(13.0, 2.0, 3.0);
    let lookat = Point::new(0.0, 0.0, 0.0);
    let camera = Camera::new(
        lookfrom,
        lookat,
        Vector::new(0.0, 1.0, 0.0),
        20.0,
        ASPECT_RATIO,
        0.1,
        10.0,
    );

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
