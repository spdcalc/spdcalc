//! # KTP CrystalType
//!
//! [Reference](http://dx.doi.org/10.1063/1.1668320)
//! [More Information](https://web.archive.org/web/20240225020609/http://www.redoptronics.com/KTP-crystal.html)
//!
//! # Example
//! ```
//! use spdcalc::{crystal::*, dim::ucum, na::Vector3, utils::*};
//! let nm = spdcalc::dim::f64prefixes::NANO * ucum::M;
//! let indices = CrystalType::KTP.get_indices(720.0 * nm, from_celsius_to_kelvin(30.));
//! let expected = ucum::Unitless::new(Vector3::new(
//!   1.7540699746332105,
//!   1.7625839942396933,
//!   1.8533562248650441,
//! ));
//! assert_eq!(indices, expected)
//! ```

use super::*;
use crate::utils::*;
use dim::{
  f64prefixes::MICRO,
  ucum::{Kelvin, K, M},
};

pub const META: CrystalMeta = CrystalMeta {
  id: "KTP",
  name: "KTP ref 1",
  reference_url: "http://dx.doi.org/10.1063/1.1668320",
  axis_type: OpticAxisType::PositiveBiaxial,
  point_group: PointGroup::HM_mm2,
  transmission_range: Some(ValidWavelengthRange(350e-9, 3500e-9)),
  temperature_dependence_known: true,
};

// from Newlight Photonics
const DNX: f64 = 1.1e-5;
const DNY: f64 = 1.3e-5;
const DNZ: f64 = 1.6e-5;

/// Get refractive Indices
#[allow(clippy::unreadable_literal)]
pub fn get_indices(wavelength: Wavelength, temperature: Kelvin<f64>) -> Indices {
  let lambda_sq = (wavelength / (MICRO * M)).powi(2);

  // https://web.archive.org/web/20240225020609/http://www.redoptronics.com/KTP-crystal.html
  let mut nx = (2.10468 + 0.89342 * lambda_sq / (lambda_sq - 0.04438) - 0.01036 * lambda_sq).sqrt();

  let mut ny = if *(wavelength / (MICRO * M)) < 1.2 {
    (2.14559 + 0.87629 * lambda_sq / (lambda_sq - 0.0485) - 0.01173 * lambda_sq).sqrt()
  } else {
    (2.0993 + 0.922683 * lambda_sq / (lambda_sq - 0.0467695) - 0.0138408 * lambda_sq).sqrt()
  };

  let mut nz = (1.9446 + 1.3617 * lambda_sq / (lambda_sq - 0.047) - 0.01491 * lambda_sq).sqrt();

  let f = *((temperature - from_celsius_to_kelvin(20.0)) / K);

  nx += f * DNX;
  ny += f * DNY;
  nz += f * DNZ;

  Indices::new(na::Vector3::new(nx, ny, nz))
}
