//! 100 M-name memory measurement for the §8 performance acceptance
//! ("100 M-name membership ≤ 2.5 GB RSS", measured and documented).
//!
//! Usage: `mem_100m <exact|bloom>` — inserts 100 M distinct synthetic
//! name0s (Illumina-shaped) into the structure and reports its own byte
//! estimate plus the process resident set size (VmRSS).

use fastq_purge::membership::{Bloom64, ExactSet};
use std::io::Read;

fn rss_kb() -> u64 {
    let mut s = String::new();
    std::fs::File::open("/proc/self/status")
        .unwrap()
        .read_to_string(&mut s)
        .unwrap();
    s.lines()
        .find(|l| l.starts_with("VmRSS:"))
        .unwrap()
        .split_whitespace()
        .nth(1)
        .unwrap()
        .parse()
        .unwrap()
}

fn main() {
    let mode = std::env::args()
        .nth(1)
        .expect("usage: mem_100m <exact|bloom>");
    let n = 100_000_000u64;
    let gb = |kb: u64| kb as f64 / 1024.0 / 1024.0; // KB -> GB
    let gbb = |b: u64| b as f64 / 1024.0 / 1024.0 / 1024.0; // bytes -> GB
    match mode.as_str() {
        "exact" => {
            let mut set = ExactSet::with_capacity(n as usize);
            for i in 0..n {
                set.insert_name(format!("M00150:2:AABCCC:1:1107:{i}:1082").as_bytes());
            }
            let kb = rss_kb();
            println!("exact: {} unique name0s", set.len());
            println!(
                "exact: estimated_bytes = {} ({:.2} GB)",
                set.estimated_bytes(),
                gbb(set.estimated_bytes() as u64)
            );
            println!("exact: VmRSS = {kb} KB ({:.2} GB)", gb(kb));
        }
        "bloom" => {
            let mut bloom = Bloom64::new(n, 1e-4);
            for i in 0..n {
                bloom.insert_name(format!("M00150:2:AABCCC:1:1107:{i}:1082").as_bytes());
            }
            let kb = rss_kb();
            println!("bloom: {} name0s inserted", bloom.observed());
            println!(
                "bloom: {} bits, k = {}, estimated_bytes = {} ({:.2} GB)",
                bloom.bits(),
                bloom.k(),
                bloom.estimated_bytes(),
                gbb(bloom.estimated_bytes() as u64)
            );
            println!("bloom: VmRSS = {kb} KB ({:.2} GB)", gb(kb));
        }
        _ => eprintln!("unknown mode {mode:?}"),
    }
}
