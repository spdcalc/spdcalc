#[macro_use]
extern crate criterion;

use criterion::{black_box, Criterion};
use serde_json::json;
use spdcalc::utils::from_celsius_to_kelvin;
use spdcalc::{crystal::CrystalType, jsa::JointSpectrum, math::Integrator, *};
use spdcalc::dim::{f64prefixes::*, ucum::*};

fn create_spdc_with_crystal(crystal: CrystalType) -> SPDC {
  let json_cfg = json!({
      "crystal": {
        "kind": crystal,
        "pm_type": "Type2_e_eo",
        "theta_deg": 90.0,
        "length_um": 2000.0,
        "temperature_c": 20.0,
        "counter_propagation": false
      },
      "pump": {
        "wavelength_nm": 775.0,
        "waist_um": 100.0,
        "bandwidth_nm": 5.35,
        "average_power_mw": 1.0,
        "spectrum_threshold": 0.01
      },
      "signal": {
        "wavelength_nm": 1550.0,
        "phi_deg": 0.0,
        "theta_deg": 0.0,
        "waist_um": 100.0
      },
      "idler": "auto",
      "periodic_poling": {
        "poling_period_um": "auto"
      },
      "deff_pm_per_volt": 1.0
    });
  let config: SPDCConfig = serde_json::from_value(json_cfg).unwrap();
  config.try_as_spdc().unwrap()
}

fn create_interpolated_crystal() -> CrystalType {
  use spdcalc::crystal::InterpolatedCrystal;

  // Build from BBO_1 data points
  let wavelengths = vec![1500.0, 1520.0, 1540.0, 1560.0, 1580.0, 1600.0].into_iter()
    .map(|nm| nm * NANO * M)
    .collect::<Vec<_>>();
  let indices = wavelengths.iter().map(|lambda| {
    CrystalType::BBO_1.get_indices(*lambda, from_celsius_to_kelvin(20.0))
  })
  .collect::<Vec<_>>();

  let no_values = indices.iter().map(|ind| ind.x).collect::<Vec<f64>>();
  let ne_values = indices.iter().map(|ind| ind.z).collect::<Vec<f64>>();

  let crystal = InterpolatedCrystal::new_uniaxial(wavelengths, no_values, ne_values)
    .expect("Failed to create interpolated crystal");

  CrystalType::Interpolated(crystal)
}

fn benchmark_jsi_calculation(c: &mut Criterion) {
  let mut group = c.benchmark_group("JSI-calculation");
  group.sample_size(10);

  // Create SPDC objects with different crystals
  let spdc_bbo = create_spdc_with_crystal(CrystalType::BBO_1);
  let spdc_aggas2 = create_spdc_with_crystal(CrystalType::AgGaS2_1);
  let spdc_interpolated = create_spdc_with_crystal(create_interpolated_crystal());

  // Define a realistic JSI calculation range (20x20 grid)
  let bbo_range = spdc_bbo.optimum_range(200);
  let aggas2_range = spdc_aggas2.optimum_range(200);
  let interpolated_range = spdc_interpolated.optimum_range(200);

  // Benchmark JSI calculation with BBO (hand-optimized Sellmeier)
  group.bench_function("JSI with BBO", |b| {
    b.iter(|| {
      let jsi = JointSpectrum::new(black_box(spdc_bbo.clone()), Integrator::default());
      let result = jsi.jsi_range(black_box(bbo_range));
      black_box(result)
    })
  });

  // Benchmark JSI calculation with AgGaS2 (generic Sellmeier framework)
  group.bench_function("JSI with AgGaS2", |b| {
    b.iter(|| {
      let jsi = JointSpectrum::new(black_box(spdc_aggas2.clone()), Integrator::default());
      let result = jsi.jsi_range(black_box(aggas2_range));
      black_box(result)
    })
  });

  // Benchmark JSI calculation with interpolated crystal
  group.bench_function("JSI with Interpolated", |b| {
    b.iter(|| {
      let jsi = JointSpectrum::new(black_box(spdc_interpolated.clone()), Integrator::default());
      let result = jsi.jsi_range(black_box(interpolated_range));
      black_box(result)
    })
  });

  group.finish();
}

criterion_group!(benches, benchmark_jsi_calculation);
criterion_main!(benches);
