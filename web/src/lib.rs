mod utils;
extern crate num;
extern crate spdcalc;
extern crate wasm_bindgen;

use num::traits::Pow;
use spdcalc::{crystal::Crystal, dim::f64prefixes::NANO};
use wasm_bindgen::prelude::*;
// use std::time::{Duration, SystemTime, UNIX_EPOCH};
// use std::f64::consts::PI;

// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
// allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC : wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

// lifted from the `console_log` example
#[wasm_bindgen]
extern "C" {
  #[wasm_bindgen(js_namespace = console)]
  fn log(a : &str);
}

#[allow(unused_macros)]
macro_rules! console_log {
    ($($t:tt)*) => (log(&format_args!($($t)*).to_string()))
}

#[wasm_bindgen]
pub fn browser_debug() {
  utils::set_panic_hook();
}

#[wasm_bindgen]
#[no_mangle]
pub fn speed_test(amount : usize) -> String {
  let window = web_sys::window().expect("should have a window in this context");
  let performance = window
    .performance()
    .expect("performance should be available");

  let test_values = init(amount);

  // console_log!("WASM: Testing with {} calculations", test_values.len());

  let start = performance.now();

  let mut total = 0.0;
  for _i in 0..test_values.len() {
    total += spdcalc::junk::calc(1.0 / test_values[_i]);
  }

  let end = performance.now();

  let elapsed = end - start;
  // console_log!("WASM: Calculation took {} ms", elapsed);
  return format_args!("[{}, {}]", total, elapsed).to_string();
}

#[wasm_bindgen]
pub fn reserve_array_test(size : usize) -> String {
  let window = web_sys::window().expect("should have a window in this context");
  let performance = window
    .performance()
    .expect("performance should be available");

  let start = performance.now();

  let _arr = init(size);

  let end = performance.now();
  let elapsed = end - start;

  return format_args!("[{}, {}]", size, elapsed).to_string();
}

fn gaussian(x : f64, y : f64, sigma : f64, x_o : f64, y_o : f64, a : f64) -> f64 {
  let s2 = 2.0 * sigma;
  let x2 : f64 = (x - x_o).pow(2);
  let y2 : f64 = (y - y_o).pow(2);

  a * (-x2 / s2 - y2 / s2).exp()
}

#[wasm_bindgen]
pub fn get_gaussian(width : usize, height : usize) -> Vec<f64> {
  let len = width * height;
  let mut arr = vec![];
  let sigma = width as f64 / 3.0;
  let x_o = width as f64 / 2.0;
  let y_o = height as f64 / 2.0;
  let mut x = 0_f64;
  let mut y = 0_f64;

  for _ in 0..len {
    let z = gaussian(x, y, sigma, x_o, y_o, 2.0);
    arr.push(z);

    x = x + 1.0;

    if (x as usize) >= width {
      x = 0.0;
      y += 1.0;
    }
  }

  arr
}

#[wasm_bindgen]
pub fn get_gaussian_ptr(width : usize, height : usize) -> *const f64 {
  let len = width * height;
  let mut arr = vec![];
  let sigma = width as f64 / 3.0;
  let x_o = width as f64 / 2.0;
  let y_o = height as f64 / 2.0;
  let mut x = 0_f64;
  let mut y = 0_f64;

  for _ in 0..len {
    let z = gaussian(x, y, sigma, x_o, y_o, 2.0);
    arr.push(z);

    x = x + 1.0;

    if (x as usize) >= width {
      x = 0.0;
      y += 1.0;
    }
  }

  arr.as_ptr()
}

#[wasm_bindgen]
pub fn get_indices(name : String, wavelength : f64, temperature : f64) -> Vec<f64> {
  let lamda = wavelength * spdcalc::dim::ucum::M;
  let kelvin = temperature * spdcalc::dim::ucum::K;
  let crystal = match name.as_ref() {
    "bbo" => Crystal::BBO_1,
    "ktp" => Crystal::KTP,
    "bibo" => Crystal::BiBO_1,
    "aggas2" => Crystal::AgGaS2_1,
    "liio3" => Crystal::LiIO3_1,
    _ => panic!(),
  };

  let indices = crystal.get_indices(lamda, kelvin);
  indices.iter().map(|i| *i).collect()
}

#[wasm_bindgen]
pub fn get_jsi_data(width : usize, height : usize) -> Vec<f64> {
  let defaults = spdcalc::spd::SPD::default();
  let mut params = spdcalc::spd::SPD { ..defaults };
  // params.crystal_setup.crystal = spdcalc::crystal::Crystal::BiBO_1;
  params.pp = Some(params.calc_periodic_poling());
  params.crystal_setup.theta = 0.5515891191131287 * spdcalc::dim::ucum::RAD;
  // params.assign_optimum_theta();

  let cfg = spdcalc::plotting::HistogramConfig {
    x_range : (1490.86 * NANO, 1609.14 * NANO),
    y_range : (1495.05 * NANO, 1614.03 * NANO),

    x_count : width,
    y_count : height,
  };

  spdcalc::plotting::plot_JSI(&params, &cfg)
}

// fn perf_to_system(amt: f64) -> SystemTime {
//     let secs = (amt as u64) / 1_000;
//     let nanos = ((amt as u32) % 1_000) * 1_000_000;
//     UNIX_EPOCH + Duration::new(secs, nanos)
// }

fn init(amount : usize) -> Vec<f64> {
  let mut test_values = vec![];
  test_values.reserve(amount);

  for _i in 0..amount {
    test_values.push((_i + 1) as f64);
  }

  return test_values;
}
