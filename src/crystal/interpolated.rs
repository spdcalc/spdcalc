//! # Interpolated Crystal Types
//!
//! Support for user-defined crystals specified via data points that are linearly interpolated.
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
use serde_derive::{Deserialize, Serialize};

/// An interpolated crystal with refractive index data specified as discrete points
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "name")]
pub enum InterpolatedCrystal {
    /// Uniaxial crystal with ordinary (no) and extraordinary (ne) refractive indices
    InterpolatedUniaxial {
        /// Wavelengths in nanometers
        wavelengths_nm: Vec<f64>,
        /// Ordinary refractive indices (no)
        no: Vec<f64>,
        /// Extraordinary refractive indices (ne)
        ne: Vec<f64>,
    },
    // Future: InterpolatedBiaxial variant
}

impl InterpolatedCrystal {
    /// Validate the interpolated crystal data
    ///
    /// Checks:
    /// - Array lengths match
    /// - Refractive indices are in valid range (1.0 < n < 4.0)
    /// - Wavelengths are monotonic (via Interpolator construction)
    pub fn validate(&self) -> Result<(), SPDCError> {
        match self {
            InterpolatedCrystal::InterpolatedUniaxial {
                wavelengths_nm,
                no,
                ne,
            } => {
                // Check array lengths
                if wavelengths_nm.len() != no.len() || wavelengths_nm.len() != ne.len() {
                    return Err(SPDCError::new(format!(
                        "Array length mismatch: {} wavelengths, {} no values, {} ne values",
                        wavelengths_nm.len(),
                        no.len(),
                        ne.len()
                    )));
                }

                if wavelengths_nm.is_empty() {
                    return Err(SPDCError::new(
                        "At least one data point is required".to_string(),
                    ));
                }

                // Validate refractive index ranges
                for (i, &n) in no.iter().enumerate() {
                    if n <= 1.0 || n >= 4.0 {
                        return Err(SPDCError::new(format!(
                            "no at wavelength {} nm is out of valid range (1.0, 4.0): {}",
                            wavelengths_nm[i], n
                        )));
                    }
                }

                for (i, &n) in ne.iter().enumerate() {
                    if n <= 1.0 || n >= 4.0 {
                        return Err(SPDCError::new(format!(
                            "ne at wavelength {} nm is out of valid range (1.0, 4.0): {}",
                            wavelengths_nm[i], n
                        )));
                    }
                }

                // Test that we can construct an interpolator (validates monotonicity)
                let outputs: Vec<[f64; 2]> = no
                    .iter()
                    .zip(ne.iter())
                    .map(|(&n_o, &n_e)| [n_o, n_e])
                    .collect();

                Interpolator::<2>::new(wavelengths_nm.clone(), outputs)
                    .map_err(|e| SPDCError::new(format!("Wavelength validation failed: {}", e)))?;

                Ok(())
            }
        }
    }

    /// Get refractive indices at a specific wavelength
    ///
    /// ## Arguments
    /// * `wavelength_nm` - Wavelength in nanometers
    ///
    /// ## Returns
    /// `(no, ne)` tuple of refractive indices
    ///
    /// Uses linear interpolation between data points.
    /// Out-of-range wavelengths are clamped to boundary values.
    pub fn get_indices(&self, wavelength_nm: f64) -> Result<(f64, f64), SPDCError> {
        match self {
            InterpolatedCrystal::InterpolatedUniaxial {
                wavelengths_nm,
                no,
                ne,
            } => {
                // Construct interpolator (cached version would be better for performance)
                let outputs: Vec<[f64; 2]> = no
                    .iter()
                    .zip(ne.iter())
                    .map(|(&n_o, &n_e)| [n_o, n_e])
                    .collect();

                let interpolator = Interpolator::<2>::new(wavelengths_nm.clone(), outputs)
                    .map_err(|e| SPDCError::new(e))?;

                let result = interpolator.interpolate(wavelength_nm);
                Ok((result[0], result[1]))
            }
        }
    }

    /// Get metadata for this crystal
    pub fn get_meta(&self) -> CrystalMeta {
        match self {
            InterpolatedCrystal::InterpolatedUniaxial { .. } => CrystalMeta {
                id: "InterpolatedUniaxial",
                name: "Interpolated Uniaxial Crystal",
                reference_url: "Custom data",
                axis_type: OpticAxisType::PositiveUniaxial,
                point_group: PointGroup::HM_mm2,
                transmission_range: None,
                temperature_dependence_known: false,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::approx_eq;
    use serde_json;

    #[test]
    fn test_validate_valid_data() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0, 600.0],
            no: vec![1.66, 1.65, 1.64],
            ne: vec![1.55, 1.54, 1.53],
        };

        assert!(crystal.validate().is_ok());
    }

    #[test]
    fn test_validate_length_mismatch() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0, 600.0],
            no: vec![1.66, 1.65], // Only 2 values
            ne: vec![1.55, 1.54, 1.53],
        };

        let result = crystal.validate();
        assert!(result.is_err());
        assert!(result.unwrap_err().0.contains("Array length mismatch"));
    }

    #[test]
    fn test_validate_empty_data() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![],
            no: vec![],
            ne: vec![],
        };

        let result = crystal.validate();
        assert!(result.is_err());
        assert!(result.unwrap_err().0.contains("At least one data point"));
    }

    #[test]
    fn test_validate_no_out_of_range_low() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0],
            no: vec![0.5, 1.65], // 0.5 is too low
            ne: vec![1.55, 1.54],
        };

        let result = crystal.validate();
        assert!(result.is_err());
        assert!(result.unwrap_err().0.contains("out of valid range"));
        assert!(result.unwrap_err().0.contains("no"));
    }

    #[test]
    fn test_validate_no_out_of_range_high() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0],
            no: vec![1.66, 5.0], // 5.0 is too high
            ne: vec![1.55, 1.54],
        };

        let result = crystal.validate();
        assert!(result.is_err());
        assert!(result.unwrap_err().0.contains("out of valid range"));
    }

    #[test]
    fn test_validate_ne_out_of_range() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0],
            no: vec![1.66, 1.65],
            ne: vec![0.8, 1.54], // 0.8 is too low
        };

        let result = crystal.validate();
        assert!(result.is_err());
        assert!(result.unwrap_err().0.contains("out of valid range"));
        assert!(result.unwrap_err().0.contains("ne"));
    }

    #[test]
    fn test_validate_duplicate_wavelengths() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0, 500.0], // Duplicate 500.0
            no: vec![1.66, 1.65, 1.64],
            ne: vec![1.55, 1.54, 1.53],
        };

        let result = crystal.validate();
        assert!(result.is_err());
        assert!(result.unwrap_err().0.contains("Wavelength validation failed"));
    }

    #[test]
    fn test_get_indices_exact_match() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0, 600.0],
            no: vec![1.66, 1.65, 1.64],
            ne: vec![1.55, 1.54, 1.53],
        };

        let (no, ne) = crystal.get_indices(500.0).unwrap();
        assert_eq!(no, 1.65);
        assert_eq!(ne, 1.54);
    }

    #[test]
    fn test_get_indices_interpolation() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 600.0],
            no: vec![1.66, 1.64],
            ne: vec![1.55, 1.53],
        };

        let (no, ne) = crystal.get_indices(500.0).unwrap(); // Midpoint
        assert!(approx_eq!(f64, no, 1.65, epsilon = 1e-10));
        assert!(approx_eq!(f64, ne, 1.54, epsilon = 1e-10));
    }

    #[test]
    fn test_get_indices_extrapolation_below() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![500.0, 600.0],
            no: vec![1.65, 1.64],
            ne: vec![1.54, 1.53],
        };

        let (no, ne) = crystal.get_indices(400.0).unwrap(); // Below range
        assert_eq!(no, 1.65); // Should clamp to first value
        assert_eq!(ne, 1.54);
    }

    #[test]
    fn test_get_indices_extrapolation_above() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0],
            no: vec![1.66, 1.65],
            ne: vec![1.55, 1.54],
        };

        let (no, ne) = crystal.get_indices(700.0).unwrap(); // Above range
        assert_eq!(no, 1.65); // Should clamp to last value
        assert_eq!(ne, 1.54);
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
            InterpolatedCrystal::InterpolatedUniaxial {
                wavelengths_nm,
                no,
                ne,
            } => {
                assert_eq!(wavelengths_nm, vec![400.0, 500.0, 600.0]);
                assert_eq!(no, vec![1.66, 1.65, 1.64]);
                assert_eq!(ne, vec![1.55, 1.54, 1.53]);
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

        let crystal: InterpolatedCrystal = serde_json::from_str(json).unwrap();
        assert!(crystal.validate().is_ok());
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
        assert!(crystal.validate().is_ok());

        let (no, ne) = crystal.get_indices(575.0).unwrap();
        // Should interpolate between 550 and 600
        assert!(no > 1.64 && no < 1.65);
        assert!(ne > 1.54 && ne < 1.55);
    }

    #[test]
    fn test_serialize_roundtrip() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0, 600.0],
            no: vec![1.66, 1.65, 1.64],
            ne: vec![1.55, 1.54, 1.53],
        };

        let json = serde_json::to_string(&crystal).unwrap();
        let deserialized: InterpolatedCrystal = serde_json::from_str(&json).unwrap();

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

    #[test]
    fn test_get_meta() {
        let crystal = InterpolatedCrystal::InterpolatedUniaxial {
            wavelengths_nm: vec![400.0, 500.0],
            no: vec![1.66, 1.65],
            ne: vec![1.55, 1.54],
        };

        let meta = crystal.get_meta();
        assert_eq!(meta.id, "InterpolatedUniaxial");
        assert_eq!(meta.temperature_dependence_known, false);
        assert_eq!(meta.transmission_range, None);
    }
}
