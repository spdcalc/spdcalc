//! # AgGaS2_1 CrystalType
//!
//! [Reference](https://web.archive.org/web/20240225020634/http://www.redoptronics.com/AgGaS2-AgGaSe2.html)
//!
//! ## Example
//! ```
//! use spdcalc::{crystal::*, dim::ucum, na::Vector3, utils::*};
//! let nm = spdcalc::dim::f64prefixes::NANO * ucum::M;
//! let indices = CrystalType::AgGaS2_1.get_indices(720.0 * nm, from_celsius_to_kelvin(30.));
//! let expected = ucum::Unitless::new(Vector3::new(
//!   2.5146353236904937,
//!   2.5146353236904937,
//!   2.46323810043117,
//! ));
//! assert_eq!(indices, expected)
//! ```
use super::*;

use crate::crystal::CrystalMeta;
use sellmeier::{equations::SellmeierStandard, temperature_dependence::Standard, SellmeierCrystal};

#[allow(non_upper_case_globals)]
pub const AgGaS2_1: SellmeierCrystal<SellmeierStandard, Standard> = SellmeierCrystal {
  meta: CrystalMeta {
    id: "AgGaS2_1",
    name: "AgGaS2 ref 1",
    reference_url: "https://www.semanticscholar.org/paper/Refractive-index-interpolation-in-phase-matching.-Bhar/63a04ace5f192a7224faa7c36c57f6ea03b99215/figure/0",
    axis_type: OpticAxisType::NegativeUniaxial,
    point_group: PointGroup::HM_4,
    transmission_range: Some(ValidWavelengthRange(500e-9, 13_000e-9)),
    temperature_dependence_known: true,
  },

  eqn: SellmeierStandard {
    a: [3.628, 3.628, 4.0172],
    b1: [2.1686, 2.1686, 1.5274],
    b2: [2.1753, 2.1753, 2.1699],
    b3: [0., 0., 0.],

    c1: [0.1003, 0.1003, 0.131],
    c2: [950.0, 950.0, 950.0],

    c3: [0., 0., 0.],
  },

  temperature_dependence: Standard {
    dn: [15.4e-5, 15.4e-5, 15.5e-5],
  },
};
