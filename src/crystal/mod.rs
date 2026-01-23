//! Defines types of crystals used in spdc.
use crate::*;

mod meta;
pub use self::meta::*;

// crystals
mod aggas2_1;
mod aggase2_1;
mod aggase2_2;
mod bbo_1;
mod bibo_1;
mod kdp_1;
mod ktp;
mod liio3_1;
mod liio3_2;
mod linb_mgo;
mod linbo3_1;

/// useful for custom crystals
pub mod sellmeier;
// use sellmeier::equations::SellmeierEquation;
// use sellmeier::temperature_dependence::TemperatureDependence;

mod crystal_type;
pub use crystal_type::*;

mod crystal_setup;
pub use crystal_setup::*;

mod polarization_type;
pub use polarization_type::*;

mod interpolated;
pub use interpolated::*;
