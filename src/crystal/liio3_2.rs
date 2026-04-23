//! # LiIO3-2 CrystalType
//!
//! K. Takizawa, M. Okada, S. Ieiri, Opt. Commun., 23, 279 (1977)
//!
//! ## Example
//! ```
//! use spdcalc::{crystal::*, dim::ucum, na::Vector3, utils::*};
//! let nm = spdcalc::dim::f64prefixes::NANO * ucum::M;
//! let indices = CrystalType::LiIO3_2.get_indices(720.0 * nm, from_celsius_to_kelvin(30.));
//! let expected = ucum::Unitless::new(Vector3::new(
//!   1.872937850009036,
//!   1.872937850009036,
//!   1.7283279244482976,
//! ));
//! assert_eq!(indices, expected)
//! ```

use super::*;
use dim::{
  f64prefixes::MICRO,
  ucum::{Kelvin, M},
};

pub const META: CrystalMeta = CrystalMeta {
  id: "LiIO3_2",
  name: "LiIO3 ref 2",
  reference_url: "https://www.sciencedirect.com/science/article/abs/pii/0030401877903261",
  axis_type: OpticAxisType::NegativeUniaxial,
  point_group: PointGroup::HM_622,
  transmission_range: Some(ValidWavelengthRange(300e-9, 5_000e-9)),
  temperature_dependence_known: false,
};

/// Get refractive Indices
#[allow(clippy::unreadable_literal)]
pub fn get_indices(wavelength: Wavelength, _temperature: Kelvin<f64>) -> Indices {
  let l_sq = (wavelength / (MICRO * M)).powi(2); // Convert for Sellmeier Coefficients

  let no = (3.4095 + 0.047664 / (l_sq - 0.033991)).sqrt();
  let ne = (2.9163 + 0.034514 / (l_sq - 0.031034)).sqrt();

  Indices::new(na::Vector3::new(no, no, ne))
}
