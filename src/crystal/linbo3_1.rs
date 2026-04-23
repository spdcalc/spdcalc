//! # LiNbO3 ref 1
//!
//! [More Information](http://www.newlightphotonics.com/v1/LN-crystal.html)
//!

use super::*;
use crate::utils::from_celsius_to_kelvin;
use dim::{
  f64prefixes::MICRO,
  ucum::{self, K, M},
};

pub const META: CrystalMeta = CrystalMeta {
  id: "LiNbO3_1",
  name: "LiNbO3 1",
  reference_url: "http://www.newlightphotonics.com/v1/LN-crystal.html",
  axis_type: OpticAxisType::NegativeUniaxial,
  point_group: PointGroup::HM_3m,
  // Corrected from spdcalc v1 based on newlight photonics website
  transmission_range: Some(ValidWavelengthRange(420e-9, 5_200e-9)),
  temperature_dependence_known: true,
};

// from Newlight Photonics
const DNX: f64 = -0.874e-6;
const DNY: f64 = DNX;
const DNZ: f64 = 39.073e-6;

/// Get refractive Indices
///
/// n_o^2 = 4.9048+0.11768/(\lambda^2-0.04750)-0.027169\lambda^2
/// n_e^2 = 4.5820+0.099169/(\lambda^2-0.04443)-0.021950\lambda^2
///
/// # Example
/// ```
/// use spdcalc::{crystal::*, dim::ucum, na::Vector3, utils::from_celsius_to_kelvin};
/// let nm = spdcalc::dim::f64prefixes::NANO * ucum::M;
/// let indices = CrystalType::LiNbO3_1.get_indices(720.0 * nm, from_celsius_to_kelvin(30.));
/// let expected = ucum::Unitless::new(Vector3::new(
///   2.267284807097424,
///   2.267284807097424,
///   2.18667811845247,
/// ));
/// assert_eq!(indices, expected)
/// ```
#[allow(clippy::unreadable_literal)]
pub fn get_indices(wavelength: Wavelength, temperature: ucum::Kelvin<f64>) -> Indices {
  let lambda_sq = (wavelength / (MICRO * M)).powi(2);

  let mut nx = (4.9048 + 0.11768 / (lambda_sq - 0.04750) - 0.027169 * lambda_sq).sqrt();
  let mut ny = nx;
  let mut nz = (4.5820 + 0.099169 / (lambda_sq - 0.044432) - 0.021950 * lambda_sq).sqrt();

  let f = *((temperature - from_celsius_to_kelvin(20.0)) / K);

  nx += f * DNX;
  ny += f * DNY;
  nz += f * DNZ;

  Indices::new(na::Vector3::new(nx, ny, nz))
}
