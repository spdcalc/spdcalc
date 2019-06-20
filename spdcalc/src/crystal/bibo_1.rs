//! # BiBO ref 1
//!
//! [More Information](http://www.newlightphotonics.com/v1/bibo-properties.html)

use super::*;
use dim::ucum;
use dim::ucum::{M};
use dim::f64prefixes::{MICRO};

pub const META :CrystalMeta = CrystalMeta {
  name: "BiBO",
  reference_url: "http://www.newlightphotonics.com/v1/bibo-properties.html",
  axis_type: OpticAxisType::PositiveBiaxial,
  point_group: PointGroup::HM_2,
  temperature_dependence_known: false,
};

/// Get refractive Indices
///
/// # Example
/// ```
/// use spdcalc::crystal::*;
/// use spdcalc::utils::dim_vector3;
/// use spdcalc::dim::ucum;
/// let nm = spdcalc::dim::f64prefixes::NANO * ucum::M;
/// let indices = bibo_1::get_indices( 720.0 * nm, 293.0 * ucum::K );
/// let expected = dim_vector3(ucum::ONE, &[1.770147077637903, 1.7990347340642352, 1.93622182289392]);
/// assert_eq!(indices, expected)
/// ```
#[allow(clippy::unreadable_literal)]
pub fn get_indices( wavelength :Wavelength, _temperature :ucum::Kelvin<f64> ) -> Indices {
  let wl_sq = (wavelength / (MICRO * M)).powi(2);

  let nx = ucum::ONE * (
    3.0740
    + 0.0323 / (wl_sq - 0.0316)
    - 0.01337 * wl_sq
  ).sqrt();

  let ny = ucum::ONE * (
    3.1685
    + 0.0373 / (wl_sq - 0.0346)
    - 0.01750 * wl_sq
  ).sqrt();

  let nz = ucum::ONE * (
    3.6545
    + 0.0511 / (wl_sq - 0.0371)
    - 0.0226 * wl_sq
  ).sqrt();

  Indices::new(nx, ny, nz)
}
