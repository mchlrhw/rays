use std::fmt;

use indicatif::ProgressBar;

const IMAGE_WIDTH: u64 = 256;
const IMAGE_HEIGHT: u64 = 256;

struct Rgb {
    r: f64,
    g: f64,
    b: f64,
}

impl fmt::Display for Rgb {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self { r, g, b } = self;

        let ir = (255.999 * r) as u64;
        let ig = (255.999 * g) as u64;
        let ib = (255.999 * b) as u64;

        write!(f, "{} {} {}", ir, ig, ib)
    }
}

fn main() {
    let progress = ProgressBar::new(IMAGE_HEIGHT);

    println!("P3\n{} {}\n255", IMAGE_WIDTH, IMAGE_HEIGHT);

    for j in (0..IMAGE_HEIGHT).rev() {
        progress.inc(1);
        for i in 0..IMAGE_WIDTH {
            let pixel_colour = Rgb {
                r: i as f64 / (IMAGE_WIDTH - 1) as f64,
                g: j as f64 / (IMAGE_HEIGHT - 1) as f64,
                b: 0.25,
            };

            println!("{}", pixel_colour);
        }
    }

    progress.finish();
}
