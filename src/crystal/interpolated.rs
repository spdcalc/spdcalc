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
use crate::utils::vacuum_wavelength_to_frequency;
use crate::SPDCError;
use dim::f64prefixes::{NANO, TERA};
use dim::ucum::{HZ, Hertz, Kelvin, M, RAD};
use serde_derive::{Deserialize, Serialize};
use std::f64::consts::TAU; // 2π — converts angular frequency (rad/s) to ordinary frequency (Hz)

/// Data struct for serialization (plain f64 values).
///
/// Accepts `wavelengths_nm` for backward-compatible input.
/// Serializes as `frequencies_thz` for exact round-trips.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct InterpolatedUniaxialData {
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  wavelengths_nm: Vec<f64>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  frequencies_thz: Vec<f64>,
  no: Vec<f64>,
  ne: Vec<f64>,
}

/// Runtime struct with interpolator indexed by frequency (THz, ascending)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(
  try_from = "InterpolatedUniaxialData",
  into = "InterpolatedUniaxialData"
)]
pub struct InterpolatedUniaxialImpl {
  interpolator: Interpolator<2>, // ascending THz → [no, ne]
}

impl InterpolatedUniaxialImpl {
  /// Create a new uniaxial interpolated crystal from wavelength and refractive index data.
  ///
  /// Wavelengths are expected in ascending order
  pub fn try_new(
    wavelengths: Vec<Wavelength>,
    no: Vec<f64>,
    ne: Vec<f64>,
  ) -> Result<Self, SPDCError> {
    // Convert ascending wavelengths → descending freq → reverse to ascending THz
    let freq_vals_asc: Vec<f64> = wavelengths
      .iter()
      .rev()
      .map(|&wl| *(vacuum_wavelength_to_frequency(wl) / (TAU * TERA * HZ * RAD)))
      .collect();
    let no_asc: Vec<f64> = no.into_iter().rev().collect();
    let ne_asc: Vec<f64> = ne.into_iter().rev().collect();
    Self::try_new_raw(freq_vals_asc, no_asc, ne_asc)
  }

  /// Create a new uniaxial interpolated crystal from angular frequency values (rad/s)
  /// and refractive index data. Frequencies must be provided in ascending order.
  ///
  /// Converts ω to ordinary frequency via f [THz] = ω / (2π × 10¹²).
  pub fn try_new_from_angular_frequencies(
    frequencies: Vec<Frequency>,
    no: Vec<f64>,
    ne: Vec<f64>,
  ) -> Result<Self, SPDCError> {
    let freqs = frequencies
      .iter()
      .map(|&f| *(f / (TAU * TERA * HZ * RAD)))
      .collect();
    Self::try_new_raw(freqs, no, ne)
  }

  /// Create a new uniaxial interpolated crystal from ordinary frequency values (Hz)
  /// and refractive index data. Frequencies must be provided in ascending order.
  pub fn try_new_from_ordinary_frequencies(
    frequencies: Vec<Hertz<f64>>,
    no: Vec<f64>,
    ne: Vec<f64>,
  ) -> Result<Self, SPDCError> {
    let freqs = frequencies
      .iter()
      .map(|&f| *(f / (HZ * TERA)))
      .collect();
    Self::try_new_raw(freqs, no, ne)
  }

  /// Create a new uniaxial interpolated crystal from ordinary frequencies in THz
  /// (f = ω / 2π, must be ascending) and refractive index data.
  pub fn try_new_raw(
    frequencies_thz: Vec<f64>,
    no: Vec<f64>,
    ne: Vec<f64>,
  ) -> Result<Self, SPDCError> {
    if frequencies_thz.len() != no.len() || frequencies_thz.len() != ne.len() {
      return Err(SPDCError::new(format!(
        "Array length mismatch: {} wavelengths, {} no values, {} ne values",
        frequencies_thz.len(),
        no.len(),
        ne.len()
      )));
    }

    if frequencies_thz.is_empty() {
      return Err(SPDCError::new(
        "At least one data point is required".to_string(),
      ));
    }

    if no.iter().any(|&n| n < 1.0) {
      return Err(SPDCError::new(
        "no values must be greater than or equal to 1.0".to_string(),
      ));
    }

    if ne.iter().any(|&n| n < 1.0) {
      return Err(SPDCError::new(
        "ne values must be greater than or equal to 1.0".to_string(),
      ));
    }

    let outputs: Vec<[f64; 2]> = no.iter().zip(&ne).map(|(&n_o, &n_e)| [n_o, n_e]).collect();

    let interpolator = Interpolator::<2>::new(frequencies_thz, outputs)
      .map_err(|e| SPDCError::new(format!("Wavelength validation failed: {}", e)))?;

    Ok(Self { interpolator })
  }
}

impl TryFrom<InterpolatedUniaxialData> for InterpolatedUniaxialImpl {
  type Error = SPDCError;

  fn try_from(data: InterpolatedUniaxialData) -> Result<Self, Self::Error> {
    if !data.frequencies_thz.is_empty() {
      Self::try_new_raw(data.frequencies_thz, data.no, data.ne)
    } else if !data.wavelengths_nm.is_empty() {
      let wavelengths: Vec<Wavelength> = data
        .wavelengths_nm
        .into_iter()
        .map(|wl_nm| wl_nm * NANO * M)
        .collect();
      Self::try_new(wavelengths, data.no, data.ne)
    } else {
      Err(SPDCError::new(
        "Must provide either frequencies_thz or wavelengths_nm",
      ))
    }
  }
}

impl From<&InterpolatedUniaxialImpl> for InterpolatedUniaxialData {
  fn from(value: &InterpolatedUniaxialImpl) -> Self {
    let (frequencies_thz, (no, ne)) = value
      .interpolator
      .iter()
      .map(|(&f, output)| (f, (output[0], output[1])))
      .unzip();
    Self {
      wavelengths_nm: vec![],
      frequencies_thz,
      no,
      ne,
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
    let freq_val = *(vacuum_wavelength_to_frequency(wavelength) / (TAU * TERA * HZ * RAD));
    let result = self.interpolator.interpolate(freq_val);
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
  /// Create a new uniaxial interpolated crystal from wavelength and refractive index data
  ///
  /// # Arguments
  /// * `wavelengths_nm` - Wavelengths in nanometers (must be sorted, will be sorted automatically)
  /// * `no` - Ordinary refractive indices (nx = ny = no for uniaxial)
  /// * `ne` - Extraordinary refractive indices (nz = ne for uniaxial)
  ///
  /// # Example
  /// ```
  /// use spdcalc::crystal::InterpolatedCrystal;
  /// use spdcalc::dim::f64prefixes::NANO;
  /// use spdcalc::dim::ucum::M;
  ///
  /// let crystal = InterpolatedCrystal::new_uniaxial(
  ///   vec![400.0 * NANO * M, 500.0 * NANO * M, 600.0 * NANO * M],
  ///   vec![1.66, 1.65, 1.64],
  ///   vec![1.55, 1.54, 1.53],
  /// ).unwrap();
  /// ```
  pub fn new_uniaxial(
    wavelengths_nm: Vec<Wavelength>,
    no: Vec<f64>,
    ne: Vec<f64>,
  ) -> Result<Self, SPDCError> {
    let inner = InterpolatedUniaxialImpl::try_new(wavelengths_nm, no, ne)?;
    Ok(InterpolatedCrystal::InterpolatedUniaxial(inner))
  }

  /// Create a new uniaxial interpolated crystal from angular frequency (rad/s) and
  /// refractive index data. Frequencies must be in ascending order.
  pub fn new_uniaxial_from_angular_frequencies(
    frequencies: Vec<Frequency>,
    no: Vec<f64>,
    ne: Vec<f64>,
  ) -> Result<Self, SPDCError> {
    let inner =
      InterpolatedUniaxialImpl::try_new_from_angular_frequencies(frequencies, no, ne)?;
    Ok(InterpolatedCrystal::InterpolatedUniaxial(inner))
  }

  /// Create a new uniaxial interpolated crystal from ordinary frequency (Hz) and
  /// refractive index data. Frequencies must be in ascending order.
  pub fn new_uniaxial_from_ordinary_frequencies(
    frequencies: Vec<Hertz<f64>>,
    no: Vec<f64>,
    ne: Vec<f64>,
  ) -> Result<Self, SPDCError> {
    let inner =
      InterpolatedUniaxialImpl::try_new_from_ordinary_frequencies(frequencies, no, ne)?;
    Ok(InterpolatedCrystal::InterpolatedUniaxial(inner))
  }

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
    let result = InterpolatedCrystal::new_uniaxial(
      vec![400.0 * NANO * M, 500.0 * NANO * M, 600.0 * NANO * M],
      vec![1.66, 1.65, 1.64],
      vec![1.55, 1.54, 1.53],
    );
    assert!(result.is_ok());
  }

  #[test]
  fn test_validate_length_mismatch() {
    let result = InterpolatedCrystal::new_uniaxial(
      vec![400.0 * NANO * M, 500.0 * NANO * M, 600.0 * NANO * M],
      vec![1.66, 1.65], // Mismatched length
      vec![1.55, 1.54, 1.53],
    );
    assert!(result.is_err());
    assert!(result
      .unwrap_err()
      .to_string()
      .contains("Array length mismatch"));
  }

  #[test]
  fn test_validate_empty_data() {
    let result = InterpolatedCrystal::new_uniaxial(vec![], vec![], vec![]);
    assert!(result.is_err());
    assert!(result
      .unwrap_err()
      .to_string()
      .contains("At least one data point"));
  }

  #[test]
  fn test_validate_no_out_of_range_low() {
    let result = InterpolatedCrystal::new_uniaxial(
      vec![400.0 * NANO * M, 500.0 * NANO * M],
      vec![0.5, 1.65], // Invalid: 0.5 < 1.0
      vec![1.55, 1.54],
    );
    assert!(result.is_err());
    let err = result.unwrap_err().to_string();
    assert!(err.contains("greater than or equal to 1.0"));
    assert!(err.contains("no"));
  }

  #[test]
  fn test_validate_ne_out_of_range() {
    let result = InterpolatedCrystal::new_uniaxial(
      vec![400.0 * NANO * M, 500.0 * NANO * M],
      vec![1.66, 1.65],
      vec![0.8, 1.54], // Invalid: 0.8 < 1.0
    );
    assert!(result.is_err());
    let err = result.unwrap_err().to_string();
    assert!(err.contains("greater than or equal to 1.0"));
    assert!(err.contains("ne"));
  }

  #[test]
  fn test_validate_duplicate_wavelengths() {
    let result = InterpolatedCrystal::new_uniaxial(
      vec![400.0 * NANO * M, 500.0 * NANO * M, 500.0 * NANO * M], // Duplicate wavelength
      vec![1.66, 1.65, 1.64],
      vec![1.55, 1.54, 1.53],
    );
    assert!(result.is_err());
    assert!(result
      .unwrap_err()
      .to_string()
      .contains("Wavelength validation failed"));
  }

  #[test]
  fn test_get_indices_exact_match() {
    let crystal = InterpolatedCrystal::new_uniaxial(
      vec![400.0 * NANO * M, 500.0 * NANO * M, 600.0 * NANO * M],
      vec![1.66, 1.65, 1.64],
      vec![1.55, 1.54, 1.53],
    )
    .unwrap();

    let indices = crystal.get_indices(500.0 * NANO * M, from_celsius_to_kelvin(20.0));

    assert_eq!(indices.x, 1.65); // no
    assert_eq!(indices.z, 1.54); // ne
  }

  #[test]
  fn test_get_indices_interpolation() {
    let crystal = InterpolatedCrystal::new_uniaxial(
      vec![400.0 * NANO * M, 600.0 * NANO * M],
      vec![1.66, 1.64],
      vec![1.55, 1.53],
    )
    .unwrap();

    let indices = crystal.get_indices(500.0 * NANO * M, from_celsius_to_kelvin(20.0));

    // Interpolation is linear in ordinary frequency f = c/λ (THz).
    // t = (f₅₀₀ - f₆₀₀) / (f₄₀₀ - f₆₀₀) = (1/500 - 1/600) / (1/400 - 1/600) = 2/5 = 0.4
    // no = 1.64 + 0.4 × (1.66 - 1.64) = 1.648
    // ne = 1.53 + 0.4 × (1.55 - 1.53) = 1.538
    assert!(approx_eq!(f64, indices.x, 1.648, epsilon = 1e-6)); // no
    assert!(approx_eq!(f64, indices.z, 1.538, epsilon = 1e-6)); // ne
  }

  #[test]
  fn test_get_indices_extrapolation_below() {
    let crystal = InterpolatedCrystal::new_uniaxial(
      vec![500.0 * NANO * M, 600.0 * NANO * M],
      vec![1.65, 1.64],
      vec![1.54, 1.53],
    )
    .unwrap();

    // 400 nm has higher frequency than both stored points → clamps to last value
    // (last in ascending-freq order = highest freq = shortest wavelength = 500 nm)
    let indices = crystal.get_indices(400.0 * NANO * M, from_celsius_to_kelvin(20.0));

    assert_eq!(indices.x, 1.65); // Clamped to 500 nm value
    assert_eq!(indices.z, 1.54);
  }

  #[test]
  fn test_get_indices_extrapolation_above() {
    let crystal = InterpolatedCrystal::new_uniaxial(
      vec![400.0 * NANO * M, 500.0 * NANO * M],
      vec![1.66, 1.65],
      vec![1.55, 1.54],
    )
    .unwrap();

    // 700 nm has lower frequency than both stored points → clamps to first value
    // (first in ascending-freq order = lowest freq = longest wavelength = 500 nm)
    let indices = crystal.get_indices(700.0 * NANO * M, from_celsius_to_kelvin(20.0));

    assert_eq!(indices.x, 1.65); // Clamped to 500 nm value
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
    let temp = from_celsius_to_kelvin(20.0);

    // Verify the stored data via get_indices at the known exact wavelengths
    let i400 = crystal.get_indices(400.0 * NANO * M, temp);
    assert_eq!(i400.x, 1.66); // no at 400 nm
    assert_eq!(i400.z, 1.55); // ne at 400 nm

    let i500 = crystal.get_indices(500.0 * NANO * M, temp);
    assert_eq!(i500.x, 1.65);
    assert_eq!(i500.z, 1.54);

    let i600 = crystal.get_indices(600.0 * NANO * M, temp);
    assert_eq!(i600.x, 1.64);
    assert_eq!(i600.z, 1.53);
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

    // Interpolation between 550 nm and 600 nm in ordinary frequency space (f = c/λ):
    // t = (f₅₇₅ - f₆₀₀) / (f₅₅₀ - f₆₀₀) = (1/575 - 1/600) / (1/550 - 1/600) = 11/23 ≈ 0.4783
    // no ≈ 1.64 + (11/23) × 0.01 ≈ 1.64478
    // ne ≈ 1.54 + (11/23) × 0.01 ≈ 1.54478
    let t = 11.0_f64 / 23.0;
    let expected_no = 1.64 + t * (1.65 - 1.64);
    let expected_ne = 1.54 + t * (1.55 - 1.54);
    assert!(approx_eq!(f64, indices.x, expected_no, epsilon = 1e-6));
    assert!(approx_eq!(f64, indices.z, expected_ne, epsilon = 1e-6));
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

    // json_out uses frequencies_thz (exact), so the second deserialization
    // produces identical interpolator data and PartialEq holds.
    assert_eq!(crystal, deserialized);
  }

  #[test]
  fn test_frequencies_thz_deserialize_and_roundtrip() {
    // Ordinary frequencies in THz (f = ω/2π), ascending order.
    // ~500 THz ≈ 600 nm, ~600 THz ≈ 500 nm, ~750 THz ≈ 400 nm
    let json_in = r#"{
      "name": "InterpolatedUniaxial",
      "frequencies_thz": [500.0, 600.0, 750.0],
      "no": [1.64, 1.65, 1.66],
      "ne": [1.53, 1.54, 1.55]
    }"#;

    let crystal: InterpolatedCrystal = serde_json::from_str(json_in).unwrap();

    // Serialization must emit frequencies_thz, not wavelengths_nm
    let json_out = serde_json::to_string(&crystal).unwrap();
    assert!(json_out.contains("frequencies_thz"));
    assert!(!json_out.contains("wavelengths_nm"));

    // Round-trip is exact because frequencies are stored and restored without conversion
    let deserialized: InterpolatedCrystal = serde_json::from_str(&json_out).unwrap();
    assert_eq!(crystal, deserialized);
  }

  #[test]
  fn test_frequencies_thz_interpolation() {
    use crate::utils::frequency_to_vacuum_wavelength;

    // Two data points at round ordinary frequencies (THz = cycles/s × 10⁻¹²).
    // 400 THz ≈ 750 nm, 600 THz ≈ 500 nm; midpoint is exactly 500 THz.
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "frequencies_thz": [400.0, 600.0],
      "no": [1.64, 1.66],
      "ne": [1.53, 1.55]
    }"#;

    let crystal: InterpolatedCrystal = serde_json::from_str(json).unwrap();

    // Convert the midpoint ordinary frequency (500 THz) to angular frequency then
    // to wavelength, so get_indices queries at the stored value exactly.
    let mid_wl = frequency_to_vacuum_wavelength(500.0 * TAU * TERA * HZ * RAD);
    let indices = crystal.get_indices(mid_wl, from_celsius_to_kelvin(20.0));

    assert!(approx_eq!(f64, indices.x, 1.65, epsilon = 1e-10));
    assert!(approx_eq!(f64, indices.z, 1.54, epsilon = 1e-10));
  }

  #[test]
  fn test_new_uniaxial_from_ordinary_frequencies() {
    use crate::utils::frequency_to_vacuum_wavelength;

    // 500 THz ≈ 600 nm, 600 THz ≈ 500 nm
    let crystal = InterpolatedCrystal::new_uniaxial_from_ordinary_frequencies(
      vec![500.0 * TERA * HZ, 600.0 * TERA * HZ],
      vec![1.64, 1.65],
      vec![1.53, 1.54],
    )
    .unwrap();

    // Exact match at the 600 THz data point
    let wl = frequency_to_vacuum_wavelength(600.0 * TAU * TERA * HZ * RAD);
    let indices = crystal.get_indices(wl, from_celsius_to_kelvin(20.0));
    assert_eq!(indices.x, 1.65);
    assert_eq!(indices.z, 1.54);

    // Midpoint between 500 THz and 600 THz → 550 THz → no = 1.645, ne = 1.535
    let wl_mid = frequency_to_vacuum_wavelength(550.0 * TAU * TERA * HZ * RAD);
    let indices_mid = crystal.get_indices(wl_mid, from_celsius_to_kelvin(20.0));
    assert!(approx_eq!(f64, indices_mid.x, 1.645, epsilon = 1e-10));
    assert!(approx_eq!(f64, indices_mid.z, 1.535, epsilon = 1e-10));
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
