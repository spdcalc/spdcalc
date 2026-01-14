//! # Interpolated Crystal Types
//!
//! Support for user-defined crystals specified via data points
//! that are linearly interpolated.
//!
//! ## Example
//! ```
//! use spdcalc::crystal::InterpolatedCrystal;
//! use serde_json;
//!
//! let json = r#"{
//!   "name": "InterpolatedUniaxial",
//!   "wavelengths_nm": [400, 500, 600],
//!   "no": [1.66, 1.65, 1.64],
//!   "ne": [1.55, 1.54, 1.53]
//! }"#;
//!
//! let crystal: InterpolatedCrystal = serde_json::from_str(json).unwrap();
//! ```

use super::*;
use crate::math::Interpolator;
use crate::SPDCError;
use dim::f64prefixes::NANO;
use dim::ucum::{Kelvin, M};
use serde_derive::{Deserialize, Serialize};

/// Data struct for serialization (plain f64 values)
#[derive(Debug, Clone, Serialize, Deserialize)]
struct InterpolatedUniaxialData {
  wavelengths_nm: Vec<f64>,
  no: Vec<f64>,
  ne: Vec<f64>,
}

/// Runtime struct with dimensional types and cached interpolator
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(
  try_from = "InterpolatedUniaxialData",
  into = "InterpolatedUniaxialData"
)]
pub struct InterpolatedUniaxialImpl {
  wavelengths: Vec<Wavelength>,
  no: Vec<f64>, // Unitless refractive indices
  ne: Vec<f64>,
  #[serde(skip)]
  interpolator: Interpolator<2>, // Cached
}

impl PartialEq for InterpolatedUniaxialImpl {
  fn eq(&self, other: &Self) -> bool {
    // Compare only the data, not the cached interpolator
    self.no == other.no && self.ne == other.ne && self.wavelengths == other.wavelengths
  }
}

impl TryFrom<InterpolatedUniaxialData> for InterpolatedUniaxialImpl {
  type Error = SPDCError;

  fn try_from(data: InterpolatedUniaxialData) -> Result<Self, Self::Error> {
    // Validate array lengths
    if data.wavelengths_nm.len() != data.no.len() || data.wavelengths_nm.len() != data.ne.len() {
      return Err(SPDCError::new(format!(
        "Array length mismatch: {} wavelengths, {} no values, {} ne values",
        data.wavelengths_nm.len(),
        data.no.len(),
        data.ne.len()
      )));
    }

    if data.wavelengths_nm.is_empty() {
      return Err(SPDCError::new(
        "At least one data point is required".to_string(),
      ));
    }

    // Validate refractive index ranges
    for (i, &n) in data.no.iter().enumerate() {
      if n <= 1.0 {
        return Err(SPDCError::new(format!(
          "no at wavelength {} nm must be greater than 1.0: {}",
          data.wavelengths_nm[i], n
        )));
      }
    }

    for (i, &n) in data.ne.iter().enumerate() {
      if n <= 1.0 {
        return Err(SPDCError::new(format!(
          "ne at wavelength {} nm must be greater than 1.0: {}",
          data.wavelengths_nm[i], n
        )));
      }
    }

    // Convert wavelengths to dimensional types
    let wavelengths = data
      .wavelengths_nm
      .iter()
      .map(|&wl_nm| wl_nm * NANO * M)
      .collect();

    // Construct interpolator with paired (no, ne) outputs
    let outputs: Vec<[f64; 2]> = data
      .no
      .iter()
      .zip(&data.ne)
      .map(|(&n_o, &n_e)| [n_o, n_e])
      .collect();

    let interpolator = Interpolator::<2>::new(data.wavelengths_nm.clone(), outputs)
      .map_err(|e| SPDCError::new(format!("Wavelength validation failed: {}", e)))?;

    Ok(Self {
      wavelengths,
      no: data.no,
      ne: data.ne,
      interpolator,
    })
  }
}

impl From<&InterpolatedUniaxialImpl> for InterpolatedUniaxialData {
  fn from(value: &InterpolatedUniaxialImpl) -> Self {
    Self {
      wavelengths_nm: value
        .wavelengths
        .iter()
        .map(|&wl| *(wl / (NANO * M)))
        .collect(),
      no: value.no.clone(),
      ne: value.ne.clone(),
    }
  }
}

impl From<InterpolatedUniaxialImpl> for InterpolatedUniaxialData {
  fn from(value: InterpolatedUniaxialImpl) -> Self {
    Self::from(&value)
  }
}

impl InterpolatedUniaxialImpl {
  /// Get refractive indices at a specific wavelength
  pub fn get_indices(&self, wavelength: Wavelength, _temperature: Kelvin<f64>) -> Indices {
    let wavelength_nm = *(wavelength / (NANO * M));
    let result = self.interpolator.interpolate(wavelength_nm);
    // For uniaxial: (no, no, ne) where no=nx=ny, ne=nz
    Indices::new(na::Vector3::new(result[0], result[0], result[1]))
  }

  /// Get metadata for this crystal
  pub fn get_meta(&self) -> CrystalMeta {
    CrystalMeta {
      id: "InterpolatedUniaxial",
      name: "Interpolated Uniaxial Crystal",
      reference_url: "",
      axis_type: OpticAxisType::PositiveUniaxial,
      point_group: PointGroup::HM_mm2,
      transmission_range: None,
      temperature_dependence_known: false,
    }
  }
}

/// An interpolated crystal with refractive index data specified as discrete points
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "name")]
pub enum InterpolatedCrystal {
  /// Uniaxial crystal with ordinary (no) and extraordinary (ne) refractive indices
  InterpolatedUniaxial(InterpolatedUniaxialImpl),
  // Future: InterpolatedBiaxial variant
}

impl InterpolatedCrystal {
  /// Get refractive indices at a specific wavelength
  pub fn get_indices(&self, wavelength: Wavelength, temperature: Kelvin<f64>) -> Indices {
    match self {
      InterpolatedCrystal::InterpolatedUniaxial(inner) => {
        inner.get_indices(wavelength, temperature)
      }
    }
  }

  /// Get metadata for this crystal
  pub fn get_meta(&self) -> CrystalMeta {
    match self {
      InterpolatedCrystal::InterpolatedUniaxial(inner) => inner.get_meta(),
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::utils::from_celsius_to_kelvin;
  use float_cmp::approx_eq;
  use serde_json;

  #[test]
  fn test_validate_valid_data() {
    let json = r#"{
      "wavelengths_nm": [400.0, 500.0, 600.0],
      "no": [1.66, 1.65, 1.64],
      "ne": [1.55, 1.54, 1.53]
    }"#;

    let result: Result<InterpolatedUniaxialImpl, _> = serde_json::from_str(json);
    assert!(result.is_ok());
  }

  #[test]
  fn test_validate_length_mismatch() {
    let json = r#"{
      "wavelengths_nm": [400.0, 500.0, 600.0],
      "no": [1.66, 1.65],
      "ne": [1.55, 1.54, 1.53]
    }"#;

    let result: Result<InterpolatedUniaxialImpl, _> = serde_json::from_str(json);
    assert!(result.is_err());
    assert!(result
      .unwrap_err()
      .to_string()
      .contains("Array length mismatch"));
  }

  #[test]
  fn test_validate_empty_data() {
    let json = r#"{
      "wavelengths_nm": [],
      "no": [],
      "ne": []
    }"#;

    let result: Result<InterpolatedUniaxialImpl, _> = serde_json::from_str(json);
    assert!(result.is_err());
    assert!(result
      .unwrap_err()
      .to_string()
      .contains("At least one data point"));
  }

  #[test]
  fn test_validate_no_out_of_range_low() {
    let json = r#"{
      "wavelengths_nm": [400.0, 500.0],
      "no": [0.5, 1.65],
      "ne": [1.55, 1.54]
    }"#;

    let result: Result<InterpolatedUniaxialImpl, _> = serde_json::from_str(json);
    assert!(result.is_err());
    let err = result.unwrap_err().to_string();
    assert!(err.contains("must be greater than 1.0"));
    assert!(err.contains("no"));
  }

  #[test]
  fn test_validate_ne_out_of_range() {
    let json = r#"{
      "wavelengths_nm": [400.0, 500.0],
      "no": [1.66, 1.65],
      "ne": [0.8, 1.54]
    }"#;

    let result: Result<InterpolatedUniaxialImpl, _> = serde_json::from_str(json);
    assert!(result.is_err());
    let err = result.unwrap_err().to_string();
    assert!(err.contains("must be greater than 1.0"));
    assert!(err.contains("ne"));
  }

  #[test]
  fn test_validate_duplicate_wavelengths() {
    let json = r#"{
      "wavelengths_nm": [400.0, 500.0, 500.0],
      "no": [1.66, 1.65, 1.64],
      "ne": [1.55, 1.54, 1.53]
    }"#;

    let result: Result<InterpolatedUniaxialImpl, _> = serde_json::from_str(json);
    assert!(result.is_err());
    assert!(result
      .unwrap_err()
      .to_string()
      .contains("Wavelength validation failed"));
  }

  #[test]
  fn test_get_indices_exact_match() {
    let json = r#"{
      "wavelengths_nm": [400.0, 500.0, 600.0],
      "no": [1.66, 1.65, 1.64],
      "ne": [1.55, 1.54, 1.53]
    }"#;

    let crystal: InterpolatedUniaxialImpl = serde_json::from_str(json).unwrap();
    let indices = crystal.get_indices(500.0 * NANO * M, from_celsius_to_kelvin(20.0));

    assert_eq!(indices.x, 1.65); // no
    assert_eq!(indices.z, 1.54); // ne
  }

  #[test]
  fn test_get_indices_interpolation() {
    let json = r#"{
      "wavelengths_nm": [400.0, 600.0],
      "no": [1.66, 1.64],
      "ne": [1.55, 1.53]
    }"#;

    let crystal: InterpolatedUniaxialImpl = serde_json::from_str(json).unwrap();
    let indices = crystal.get_indices(500.0 * NANO * M, from_celsius_to_kelvin(20.0));

    assert!(approx_eq!(f64, indices.x, 1.65, epsilon = 1e-10)); // no
    assert!(approx_eq!(f64, indices.z, 1.54, epsilon = 1e-10)); // ne
  }

  #[test]
  fn test_get_indices_extrapolation_below() {
    let json = r#"{
      "wavelengths_nm": [500.0, 600.0],
      "no": [1.65, 1.64],
      "ne": [1.54, 1.53]
    }"#;

    let crystal: InterpolatedUniaxialImpl = serde_json::from_str(json).unwrap();
    let indices = crystal.get_indices(400.0 * NANO * M, from_celsius_to_kelvin(20.0));

    assert_eq!(indices.x, 1.65); // Clamped to first value
    assert_eq!(indices.z, 1.54);
  }

  #[test]
  fn test_get_indices_extrapolation_above() {
    let json = r#"{
      "wavelengths_nm": [400.0, 500.0],
      "no": [1.66, 1.65],
      "ne": [1.55, 1.54]
    }"#;

    let crystal: InterpolatedUniaxialImpl = serde_json::from_str(json).unwrap();
    let indices = crystal.get_indices(700.0 * NANO * M, from_celsius_to_kelvin(20.0));

    assert_eq!(indices.x, 1.65); // Clamped to last value
    assert_eq!(indices.z, 1.54);
  }

  // Serde tests

  #[test]
  fn test_deserialize_with_tag() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [400.0, 500.0, 600.0],
      "no": [1.66, 1.65, 1.64],
      "ne": [1.55, 1.54, 1.53]
    }"#;

    let crystal: InterpolatedCrystal = serde_json::from_str(json).unwrap();

    match crystal {
      InterpolatedCrystal::InterpolatedUniaxial(inner) => {
        assert_eq!(inner.no, vec![1.66, 1.65, 1.64]);
        assert_eq!(inner.ne, vec![1.55, 1.54, 1.53]);
      }
    }
  }

  #[test]
  fn test_deserialize_minimal_data() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [500.0, 600.0],
      "no": [1.65, 1.64],
      "ne": [1.54, 1.53]
    }"#;

    let result: Result<InterpolatedCrystal, _> = serde_json::from_str(json);
    assert!(result.is_ok());
  }

  #[test]
  fn test_deserialize_many_points() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [400, 450, 500, 550, 600, 650, 700, 750, 800],
      "no": [1.70, 1.68, 1.66, 1.65, 1.64, 1.63, 1.62, 1.61, 1.60],
      "ne": [1.60, 1.58, 1.56, 1.55, 1.54, 1.53, 1.52, 1.51, 1.50]
    }"#;

    let crystal: InterpolatedCrystal = serde_json::from_str(json).unwrap();
    let indices = crystal.get_indices(575.0 * NANO * M, from_celsius_to_kelvin(20.0));

    // Should interpolate between 550 and 600
    assert!(approx_eq!(f64, indices.x, 1.645, epsilon = 1e-10)); // no
    assert!(approx_eq!(f64, indices.z, 1.545, epsilon = 1e-10)); // ne
  }

  #[test]
  fn test_serialize_roundtrip() {
    let json_in = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [400.0, 500.0, 600.0],
      "no": [1.66, 1.65, 1.64],
      "ne": [1.55, 1.54, 1.53]
    }"#;

    let crystal: InterpolatedCrystal = serde_json::from_str(json_in).unwrap();
    let json_out = serde_json::to_string(&crystal).unwrap();
    let deserialized: InterpolatedCrystal = serde_json::from_str(&json_out).unwrap();

    assert_eq!(crystal, deserialized);
  }

  #[test]
  fn test_deserialize_missing_field() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [400.0, 500.0],
      "no": [1.66, 1.65]
    }"#; // Missing "ne" field

    let result: Result<InterpolatedCrystal, _> = serde_json::from_str(json);
    assert!(result.is_err());
  }

  #[test]
  fn test_deserialize_wrong_type() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": "not an array",
      "no": [1.66, 1.65],
      "ne": [1.55, 1.54]
    }"#;

    let result: Result<InterpolatedCrystal, _> = serde_json::from_str(json);
    assert!(result.is_err());
  }
}
