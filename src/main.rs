use std::time::{Duration, Instant};

use cgmath::{num_traits::Float, vec2, Vector2};

mod noise;

fn main() {
    const N: usize = 11;
    let before = Instant::now();
    time(N);
    let after = before.elapsed();
    println!("Took {} ns", after.as_nanos());
    //
    return;
    let name = "2x2_o4.png";
    let size = vec2(512, 512);
    let mut imgbuf = image::ImageBuffer::new(size.x, size.y);

    let timer = Instant::now();
    let n = noise::Perlin::new(size, 2, 2, 6, 2.0, 0.5, Some(444));
    let elapsed_time = timer.elapsed();

    println!("Took {} ns", elapsed_time.as_nanos());
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let mut val = (n.get(x as f32, y as f32) * 255.0) as u8;
        // let point = vec2(
        //     (x as f32 / size.x as f32) * 10.0,
        //     (y as f32 / size.y as f32) * 10.0,
        // );
        // if point.x % 1.0 <= 0.1 || point.y % 1.0 <= 0.1 {
        //     val = 255;
        // }
        let rgb: [u8; 3] = [val, val, val];
        *pixel = image::Rgb(rgb);

        // *pixel = image::Rgb(noise::get_color(n.get(x as f32, y as f32)))

        // *pixel = image::Rgb(noise::get_cave(n.get(x as f32, y as f32)))
    }

    println!("Saved to {}", name);
    imgbuf.save(name).unwrap();
}

fn time(o: usize) {
    let size = vec2(512, 512);
    let mut imgbuf = image::ImageBuffer::new(size.x, size.y);

    let n = noise::Perlin::new(size, 1, 1, o, 2.0, 0.5, Some(444));

    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let val = (n.get(x as f32, y as f32) * 255.0) as u8;
        let rgb: [u8; 3] = [val, val, val];
        *pixel = image::Rgb(rgb);
    }
    imgbuf.save(&format!("timer-{}.png", o)).unwrap();
}
