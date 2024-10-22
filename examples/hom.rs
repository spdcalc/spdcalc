use spdcalc::prelude::*;
use textplots::{Chart, Plot, Shape};

fn main() {
  let spdc = SPDC::from_json(
    r#"
{
  "crystal": {
    "kind": "KTP",
    "pm_type": "e->eo",
    "phi_deg": 0,
    "theta_deg": 90,
    "length_um": 14000,
    "temperature_c": 20
  },
  "pump": {
    "wavelength_nm": 775,
    "waist_um": 200,
    "bandwidth_nm": 0.5,
    "average_power_mw": 300
  },
  "signal": {
    "wavelength_nm": 1550,
    "phi_deg": 0,
    "theta_external_deg": 0,
    "waist_um": 100,
    "waist_position_um": "auto"
  },
  "idler": "auto",
  "periodic_poling": {
    "poling_period_um": "auto"
  },
  "deff_pm_per_volt": 7.6
}
  "#,
  )
  .unwrap();

  let x_values = Steps(380. * PICO * S, 390. * PICO * S, 80);
  let data = spdc.hom_rate_series(x_values, spdc.optimum_range(100), Integrator::default());
  let values: Vec<_> = x_values
    .into_iter()
    .map(|s| *(s / PICO / S))
    .map(|y| y as f32)
    .zip(data.iter().map(|x| *x as f32))
    .collect();

  println!("HOM Series");
  Chart::new(
    160,
    60,
    *(x_values.0 / PICO / S) as f32,
    *(x_values.1 / PICO / S) as f32,
  )
  .lineplot(&Shape::Points(values.as_slice()))
  .display();
}
