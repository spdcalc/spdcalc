#[macro_use]
extern crate criterion;

use criterion::{black_box, Criterion};

extern crate spdcalc;
use spdcalc::{crystal::CrystalType, dim::ucum::M};

fn create_interpolated_crystal() -> CrystalType {
  // Create an interpolated crystal with realistic data covering visible to NIR range
  // This mimics a typical BBO-like crystal
  let json = r#"{
    "name": "InterpolatedUniaxial",
    "wavelengths_nm": [400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0],
    "no": [1.6749, 1.6725, 1.6701, 1.6677, 1.6653, 1.6629, 1.6605, 1.6581, 1.6557, 1.6533, 1.6509],
    "ne": [1.5555, 1.5540, 1.5525, 1.5510, 1.5495, 1.5480, 1.5465, 1.5450, 1.5435, 1.5420, 1.5405]
  }"#;
  serde_json::from_str(json).expect("Failed to create interpolated crystal")
}

fn criterion_benchmark(c: &mut Criterion) {
  let temperature = 293.0 * spdcalc::dim::ucum::K;
  let single_wavelength = 720e-9 * M;
  let varied_wavelengths = [450e-9 * M, 550e-9 * M, 650e-9 * M, 750e-9 * M, 850e-9 * M];

  let interpolated_crystal = create_interpolated_crystal();

  // Single wavelength benchmarks - apples to apples comparison
  c.bench_function("BBO single wavelength", |b| {
    b.iter(|| CrystalType::BBO_1.get_indices(black_box(single_wavelength), temperature))
  });

  #[allow(non_snake_case)]
  c.bench_function("AgGaS2 single wavelength", |b| {
    b.iter(|| CrystalType::AgGaS2_1.get_indices(black_box(single_wavelength), temperature))
  });

  c.bench_function("Interpolated single wavelength", |b| {
    b.iter(|| interpolated_crystal.get_indices(black_box(single_wavelength), temperature))
  });

  // Varied wavelengths benchmarks - realistic usage pattern
  c.bench_function("BBO varied wavelengths", |b| {
    b.iter(|| {
      for &wl in &varied_wavelengths {
        CrystalType::BBO_1.get_indices(black_box(wl), temperature);
      }
    })
  });

  #[allow(non_snake_case)]
  c.bench_function("AgGaS2 varied wavelengths", |b| {
    b.iter(|| {
      for &wl in &varied_wavelengths {
        CrystalType::AgGaS2_1.get_indices(black_box(wl), temperature);
      }
    })
  });

  c.bench_function("Interpolated varied wavelengths", |b| {
    b.iter(|| {
      for &wl in &varied_wavelengths {
        interpolated_crystal.get_indices(black_box(wl), temperature);
      }
    })
  });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
