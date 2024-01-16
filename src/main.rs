/// Based on Bartosz Zacynski's "Draw the Mandelbrot Set in Python"
/// https://realpython.com/mandelbrot-set-python/#color-palette
use clap::Parser;
use image::{Rgb, RgbImage};
use num::complex::Complex;
use rayon::prelude::*;

#[derive(Debug)]
struct ColorInterpolator {
    colors: Vec<Rgb<u8>>,
}

impl ColorInterpolator {
    fn new(colors: Vec<Rgb<u8>>) -> Self {
        Self { colors }
    }

    /// Takes a value `t` in the domain \[0.0, 1.0\] and returns a color along the linearly
    /// interpolated curve of `self.colors`.
    fn color_at(&self, t: f64) -> Rgb<u8> {
        fn lerp(a: f64, b: f64, t: f64) -> f64 {
            a + t * (b - a)
        }

        let t = t.clamp(0.0, 1.0);
        let n = self.colors.len() as f64 - 1.0;
        let x = t * n;
        let i = x.floor();
        let j = x.ceil();
        if i == j {
            return self.colors[i as usize];
        }
        // scaling t to be in between color[i] (t_min) and color[j] (t_max)
        let t_min = i / n;
        let t_max = j / n;
        let scaled_t = (t - t_min) / (t_max - t_min);
        let [r1, g1, b1] = self.colors[i as usize].0.map(|c| c as f64);
        let [r2, g2, b2] = self.colors[j as usize].0.map(|c| c as f64);
        let r = lerp(r1, r2, scaled_t) as u8;
        let g = lerp(g1, g2, scaled_t) as u8;
        let b = lerp(b1, b2, scaled_t) as u8;

        Rgb([r, g, b])
    }
}

#[derive(Debug)]
struct MandelbrotSet {
    max_iterations: u32,
    escape_radius: f64,
}

impl MandelbrotSet {
    fn new(max_iterations: u32, escape_radius: f64) -> Self {
        Self {
            max_iterations,
            escape_radius,
        }
    }

    /// Count how many iterations it takes for a value to "escape" the Mandelbrot set
    fn escape_count(&self, c: Complex<f64>, smooth: bool) -> f64 {
        let mut z = Complex::new(0.0, 0.0);
        for i in 0..self.max_iterations {
            z = z * z + c;
            let znorm2 = z.re * z.re + z.im * z.im;
            if znorm2 > self.escape_radius * self.escape_radius {
                if smooth {
                    return i as f64 + 1.0 - znorm2.sqrt().ln().log(2.0);
                }
            }
        }
        self.max_iterations as f64
    }

    /// The ratio of escape count / max iterations
    fn stability(&self, c: Complex<f64>, smooth: bool) -> f64 {
        (self.escape_count(c, smooth) / self.max_iterations as f64).clamp(0.0, 1.0)
    }
}

#[derive(Debug)]
struct Viewport {
    image_width: usize,
    image_height: usize,
    center: Complex<f64>,
    view_width: f64,
}

impl Viewport {
    fn new(image_width: usize, image_height: usize, center: Complex<f64>, view_width: f64) -> Self {
        Self {
            image_width,
            image_height,
            center,
            view_width,
        }
    }

    fn scale(&self) -> f64 {
        self.view_width / self.image_width as f64
    }

    fn view_height(&self) -> f64 {
        self.scale() * self.image_height as f64
    }

    fn offset(&self) -> Complex<f64> {
        self.center + Complex::new(-self.view_width, self.view_height()) / 2.0
    }

    /// Given a pixel in the image plane, return the coordinates in the complex plane
    fn pixel_loc(&self, row: usize, col: usize) -> Complex<f64> {
        Complex::new(col as f64, -(row as f64)) * self.scale() + self.offset()
    }
}

#[derive(Parser)]
struct Args {
    /// Output image width
    #[arg(long)]
    image_width: usize,
    /// Output image height
    #[arg(long)]
    image_height: usize,
    /// Real part of center of image
    #[arg(long)]
    center_real: f64,
    /// Imaginary part of center of image
    #[arg(long)]
    center_imag: f64,
    /// Width of the viewport
    #[arg(short, long)]
    viewport_width: f64,
    /// Maximum iterations of pixel calculation
    #[arg(short, long)]
    max_iterations: u32,
    /// Color map for the image in hex strings. Defaults to grayscale
    #[arg(long, num_args=2.., value_delimiter=' ')]
    colors: Option<Vec<String>>,
    /// Output path for the resulting image
    #[arg(short, long)]
    output: String,
}

fn hex_to_rgb(hex: &str) -> Rgb<u8> {
    let n = i64::from_str_radix(hex, 16).unwrap();
    let r = (n >> 16) & 0xFF;
    let g = (n >> 8) & 0xFF;
    let b = n & 0xFF;
    Rgb([r as u8, g as u8, b as u8])
}

fn main() {
    let args = Args::parse();
    let colors = match args.colors {
        Some(cs) => cs.iter().map(|c| hex_to_rgb(c)).collect(),
        None => vec![Rgb([0, 0, 0]), Rgb([0xFF, 0xFF, 0xFF])],
    };
    let color_interpolator = ColorInterpolator::new(colors);
    let mandelbrot_set = MandelbrotSet::new(args.max_iterations, 2.0);
    let viewport = Viewport::new(
        args.image_width,
        args.image_height,
        Complex::new(args.center_real, args.center_imag),
        args.viewport_width,
    );
    let mut pixels: Vec<Vec<Rgb<u8>>> = Vec::new();
    for row in 0..args.image_height {
        pixels.push(
            (0..args.image_width)
                .into_par_iter()
                .map(|col| {
                    let c = viewport.pixel_loc(row, col);
                    let instability = 1.0 - mandelbrot_set.stability(c, true);
                    let color = color_interpolator.color_at(instability);
                    color
                })
                .collect(),
        );
    }
    let mut image = RgbImage::new(args.image_width as u32, args.image_height as u32);
    for row in 0..args.image_height {
        for col in 0..args.image_width {
            image.put_pixel(col as u32, row as u32, pixels[row][col]);
        }
    }
    let _ = image.save(args.output);
}
