//! Dev-only: measure pure inflate cost of a .gz file.
use std::io::Read;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: inflate_bench <file.gz>");
    let f = std::fs::File::open(&path).unwrap();
    let mut dec = flate2::read::MultiGzDecoder::new(f);
    let mut out = Vec::new();
    let t = std::time::Instant::now();
    dec.read_to_end(&mut out).unwrap();
    let dt = t.elapsed();
    println!(
        "{} bytes in {:.3}s = {:.0} MB/s",
        out.len(),
        dt.as_secs_f64(),
        out.len() as f64 / dt.as_secs_f64() / 1e6
    );
}
