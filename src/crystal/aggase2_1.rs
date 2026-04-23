//! # AgGaSe2-1
//!
//! H. Kildal, J. Mikkelsen, Opt. Commun. 9, 315 (1973)

use super::*;
use crate::utils::from_celsius_to_kelvin;
use dim::{
  f64prefixes::MICRO,
  ucum::{self, K, M},
};

pub const META: CrystalMeta = CrystalMeta {
  id: "AgGaSe2_1",
  name: "AgGaSe2 Ref 1",
  reference_url: "https://www.sciencedirect.com/science/article/pii/0030401873903167",
  axis_type: OpticAxisType::NegativeUniaxial,
  point_group: PointGroup::HM_3m,
  transmission_range: Some(ValidWavelengthRange(1_000e-9, 13_500e-9)),
  temperature_dependence_known: true,
};

// from Newlight Photonics
const DNX: f64 = 15e-5;
const DNY: f64 = DNX;
const DNZ: f64 = 15e-5;

/// Get refractive Indices
///
/// # Example
/// ```
/// use spdcalc::{crystal::*, dim::ucum, na::Vector3, utils::from_celsius_to_kelvin};
/// let nm = spdcalc::dim::f64prefixes::NANO * ucum::M;
/// let indices = CrystalType::AgGaSe2_1.get_indices(720.0 * nm, from_celsius_to_kelvin(30.));
/// let expected = ucum::Unitless::new(Vector3::new(
///   2.837020760678037,
///   2.837020760678037,
///   2.8283867598339847,
/// ));
/// assert_eq!(indices, expected)
/// ```
#[allow(clippy::unreadable_literal)]
pub fn get_indices(wavelength: Wavelength, temperature: ucum::Kelvin<f64>) -> Indices {
  let lambda = wavelength / (MICRO * M);

  let mut nx = (3.9362
    + 2.9113 / (1.0 - (0.38821 / lambda).powi(2))
    + 1.7954 / (1.0 - (40.0 / lambda).powi(2)))
  .sqrt();
  let mut ny = nx;
  let mut nz = (3.3132
    + 3.3616 / (1.0 - (0.38201 / lambda).powi(2))
    + 1.7677 / (1.0 - (40.0 / lambda).powi(2)))
  .sqrt();

  let f = *((temperature - from_celsius_to_kelvin(20.0)) / K);

  nx += f * DNX;
  ny += f * DNY;
  nz += f * DNZ;

  Indices::new(na::Vector3::new(nx, ny, nz))
}
