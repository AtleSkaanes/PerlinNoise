use std::ops::{Add, Div, Mul, Sub};

use cgmath::num_traits::Float;
use cgmath::prelude::*;
use cgmath::{vec2, Rad, Vector2};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

pub struct Perlin {
    size: Vector2<u32>,
    width: usize,
    height: usize,
    octaves: usize,
    lacunarity: f32,
    persistance: f32,
    seed: Option<u64>,
    rng: ChaCha20Rng,
    vectors: Vec<Vec<Vec<Vector2<f32>>>>,
}

impl Perlin {
    pub fn new(
        size: Vector2<u32>,
        width: usize,
        height: usize,
        octaves: usize,
        lacunarity: f32,
        persistance: f32,
        seed: Option<u64>,
    ) -> Self {
        let mut noise = Self {
            size,
            width,
            height,
            octaves,
            lacunarity,
            persistance,
            seed,
            rng: ChaCha20Rng::from_entropy(),
            vectors: vec![],
        };
        noise.generate_grid();
        noise
    }
    pub fn get(&self, x: f32, y: f32) -> f32 {
        let mut val = 0.0;
        let mut full_height = 0.0;

        for octave in 0..self.octaves {
            // Map x and y to noisemap coordinates
            let point = vec2(
                (x / self.size.x as f32) * (self.width * self.get_frequency(octave)) as f32,
                (y / self.size.y as f32) * (self.height * self.get_frequency(octave)) as f32,
            );
            // Get the corner coordinates
            let l_top = vec2(point.x.floor(), point.y.floor());
            let r_top = vec2(point.x.ceil(), point.y.floor());
            let l_bot = vec2(point.x.floor(), point.y.ceil());
            let r_bot = vec2(point.x.ceil(), point.y.ceil());

            let corners = self.get_corners(point.x, point.y, octave);
            // Get the dot products
            let dots = vec![
                corners[0].dot(l_top - point),
                corners[1].dot(r_top - point),
                corners[2].dot(l_bot - point),
                corners[3].dot(r_bot - point),
            ];
            // let dots = vec![
            //     corners[0].dot(point - l_top),
            //     corners[1].dot(point - r_top),
            //     corners[2].dot(point - l_bot),
            //     corners[3].dot(point - r_bot),
            // ];
            // Compute interpolation weights
            let weights = vec2(point.x - l_top.x as f32, point.y - l_top.y as f32);

            let x_top_val = interpolate(dots[0], dots[1], weights.x);
            let x_bot_val = interpolate(dots[2], dots[3], weights.x);

            full_height += self.get_amplitude(octave);
            val += interpolate(x_top_val, x_bot_val, weights.y) * self.get_amplitude(octave);
        }
        // Get max value for a 2D perlin noise, which is the same as sin(0.25pi)
        let val_max = f32::sqrt(2.0 / 4.0) * full_height;
        map_range((-val_max, val_max), (0.0, 1.0), val)
    }
    fn get_corners(&self, x: f32, y: f32, octave: usize) -> [Vector2<f32>; 4] {
        [
            self.vectors[octave][x.floor() as usize][y.floor() as usize],
            self.vectors[octave][x.ceil() as usize][y.floor() as usize],
            self.vectors[octave][x.floor() as usize][y.ceil() as usize],
            self.vectors[octave][x.ceil() as usize][y.ceil() as usize],
        ]
    }
    fn generate_grid(&mut self) {
        if let Some(seed) = self.seed {
            self.rng = ChaCha20Rng::seed_from_u64(seed);
        }
        for o in 0..self.octaves as usize {
            self.vectors.push(vec![]);
            for i in 0..((self.width + 1) * self.get_frequency(o)) as usize {
                self.vectors[o].push(vec![]);
                for _ in 0..((self.height + 1) * self.get_frequency(o)) as usize {
                    let angle = Rad(self.rng.gen_range(0.0..360.0));
                    self.vectors[o][i].push(vec2(Rad::cos(angle), Rad::sin(angle)))
                }
            }
        }
    }
    fn get_frequency(&self, octave: usize) -> usize {
        f32::powi(self.lacunarity, octave as i32) as usize
    }
    fn get_amplitude(&self, octave: usize) -> f32 {
        f32::powf(self.persistance as f32, octave as f32)
    }
}

pub struct Value {
    size: Vector2<u32>,
    width: usize,
    height: usize,
    octaves: usize,
    lacunarity: f32,
    persistance: f32,
    seed: Option<u64>,
    rng: ChaCha20Rng,
    vectors: Vec<Vec<Vec<f32>>>,
}

impl Value {
    pub fn new(
        size: Vector2<u32>,
        width: usize,
        height: usize,
        octaves: usize,
        lacunarity: f32,
        persistance: f32,
        seed: Option<u64>,
    ) -> Self {
        let mut noise = Self {
            size,
            width,
            height,
            octaves,
            lacunarity,
            persistance,
            seed,
            rng: ChaCha20Rng::from_entropy(),
            vectors: vec![],
        };
        noise.generate_grid();
        noise
    }
    pub fn get(&self, x: f32, y: f32) -> f32 {
        let mut val = 0.0;
        for octave in 0..self.octaves {
            // Map x and y to noisemap coordinates
            let point = vec2(
                (x / self.size.x as f32) * (self.width * self.get_frequency(octave)) as f32,
                (y / self.size.y as f32) * (self.height * self.get_frequency(octave)) as f32,
            );
            // Get the corner coordinates
            let l_top = vec2(point.x.floor(), point.y.floor());

            let corners = self.get_corners(point.x, point.y, octave);

            // Compute interpolation weights
            let weights = vec2(point.x - l_top.x as f32, point.y - l_top.y as f32);

            let x_top_val = lerp(corners[0], corners[1], weights.x);
            let x_bot_val = lerp(corners[2], corners[3], weights.x);

            val += lerp(x_top_val, x_bot_val, weights.y) * self.get_amplitude(octave);
        }
        val /= self.octaves as f32;
        let constraint = f32::sqrt(2.0 / 4.0);
        (val + constraint) / 2.0
    }
    fn get_corners(&self, x: f32, y: f32, octave: usize) -> [f32; 4] {
        [
            self.vectors[octave][x.floor() as usize][y.floor() as usize],
            self.vectors[octave][x.ceil() as usize][y.floor() as usize],
            self.vectors[octave][x.floor() as usize][y.ceil() as usize],
            self.vectors[octave][x.ceil() as usize][y.ceil() as usize],
        ]
    }
    fn generate_grid(&mut self) {
        if let Some(seed) = self.seed {
            self.rng = ChaCha20Rng::seed_from_u64(seed);
        }

        for o in 0..self.octaves as usize {
            self.vectors.push(vec![]);
            for i in 0..((self.width + 1) * self.get_frequency(o)) as usize {
                self.vectors[o].push(vec![]);
                for _ in 0..((self.height + 1) * self.get_frequency(o)) as usize {
                    self.vectors[o][i].push(self.rng.gen_range(-1.0..1.0))
                }
            }
        }
    }
    pub fn get_seed(&self) -> u64 {
        self.rng.get_stream()
    }
    fn get_frequency(&self, octave: usize) -> usize {
        f32::powi(self.lacunarity, octave as i32) as usize
    }
    fn get_amplitude(&self, octave: usize) -> f32 {
        f32::powf(self.persistance as f32, octave as f32)
    }
}

pub fn get_color(value: f32) -> [u8; 3] {
    return if value > 0.6 {
        [50, 100, 255]
    } else if value > 0.3 {
        [50, 255, 100]
    } else if value > 0.2 {
        [150, 150, 150]
    } else {
        [230, 230, 255]
    };
}

pub fn get_cave(value: f32) -> [u8; 3] {
    // return if (0.45..0.5).contains(&value) {
    //     [100, 100, 255]
    // } else if (0.5..0.55).contains(&value) {
    //     [255, 100, 100]
    // } else {
    //     [0, 0, 0]
    // };

    return if (0.45..0.5).contains(&value) {
        [255, 255, 255]
    } else if (0.5..0.55).contains(&value) {
        [255, 255, 255]
    } else {
        [0, 0, 0]
    };
}

pub fn get_mod(value: f32) -> [u8; 3] {
    return if value % 0.02 <= 0.007 {
        [255, 100, 100]
    } else {
        [0, 0, 0]
    };
}

fn interpolate(a0: f32, a1: f32, weight: f32) -> f32 {
    (a1 - a0) * (3.0 - weight * 2.0) * weight * weight + a0
}
fn lerp(a0: f32, a1: f32, weight: f32) -> f32 {
    a0 + weight * (a1 - a0)
}

fn map_range<T: Copy>(from_range: (T, T), to_range: (T, T), s: T) -> T
where
    T: Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T>,
{
    to_range.0 + (s - from_range.0) * (to_range.1 - to_range.0) / (from_range.1 - from_range.0)
}
