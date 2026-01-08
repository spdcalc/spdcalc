use super::*;
use crate::SPDCError;
use dim::f64prefixes::MICRO;
use dim::ucum::{Kelvin, K, M};
use std::fmt;
use std::str::FromStr;
use utils::from_celsius_to_kelvin;

/// A mathematical expression for a crystal's refractive indices
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum CrystalExpr {
  /// A Uniaxial crystal
  Uniaxial {
    /// The ordinary refractive index
    #[serde(skip_serializing)]
    no: meval::Expr,
    /// The extraordinary refractive index
    #[serde(skip_serializing)]
    ne: meval::Expr,
  },
  /// A Biaxial crystal
  Biaxial {
    /// The x refractive index
    #[serde(skip_serializing)]
    nx: meval::Expr,
    /// The y refractive index
    #[serde(skip_serializing)]
    ny: meval::Expr,
    /// The z refractive index
    #[serde(skip_serializing)]
    nz: meval::Expr,
  },
}

/// The type of crystal
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum CrystalType {
  /// BBO (ref 1)
  BBO_1,
  /// KTP
  KTP,

  /// BiBO_1
  BiBO_1,
  /// LiNbO3_1
  LiNbO3_1,
  /// LiNb_MgO
  LiNb_MgO,
  /// KDP_1
  KDP_1,
  /// AgGaSe2_1
  AgGaSe2_1,
  /// AgGaSe2_2
  AgGaSe2_2,

  /// LiIO3_2
  LiIO3_2,
  /// LiIO3_1
  LiIO3_1,
  /// AgGaS2_1
  AgGaS2_1,
  // Sellmeier(sellmeier::SellmeierCrystal<Q, T>),
  /// Crystal Expression
  #[serde(untagged)]
  Expr(CrystalExpr),

  /// Interpolated crystal data
  #[serde(untagged)]
  Interpolated(InterpolatedCrystal),
}

impl fmt::Display for CrystalType {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{}", self.get_meta().id)
  }
}

impl FromStr for CrystalType {
  type Err = SPDCError;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    CrystalType::from_string(s)
  }
}

impl CrystalType {
  /// Get the crystal from its id string or expression
  ///
  /// To create a crystal from an expression, use either a JSON (or HJSON) string,
  /// or an equation string. Expressions for biaxial crystals should have properties
  /// `nx`, `ny`, and `nz`, while uniaxial crystals should have `no` and `ne`.
  ///
  /// The JSON string should be of the form:
  ///
  /// ```
  /// use spdcalc::prelude::*;
  /// // BBO
  /// // no2=2.7359+0.01878/(λ2-0.01822)-0.01354λ2
  /// // ne2=2.3753+0.01224/(λ2-0.01667)-0.01516λ2
  /// // dno/dT = -9.3 x 10-6/°C
  /// // dne/dT = -16.6 x 10-6/°C
  ///
  /// let json = r#"
  /// {
  ///   "no": "sqrt(2.7359+0.01878/(l^2-0.01822)-0.01354*l^2) - 9.3e-6 * T",
  ///   "ne": "sqrt(2.3753+0.01224/(l^2-0.01667)-0.01516*l^2) - 16.6e-6 * T"
  /// }
  /// "#;
  /// let crystal = CrystalType::from_string(json).unwrap();
  /// ```
  ///
  /// The equation string should be of the form:
  ///
  /// ```
  /// use spdcalc::prelude::*;
  /// let expr = r#"
  ///   no = sqrt(2.7359+0.01878/(l^2-0.01822)-0.01354*l^2) - 9.3e-6 * T
  ///   ne = sqrt(2.3753+0.01224/(l^2-0.01667)-0.01516*l^2) - 16.6e-6 * T
  /// "#;
  /// let crystal = CrystalType::from_string(expr).unwrap();
  /// ```
  pub fn from_string(id: &str) -> Result<Self, SPDCError> {
    match id {
      "BBO_1" => Ok(CrystalType::BBO_1),
      "KTP" => Ok(CrystalType::KTP),
      "BiBO_1" => Ok(CrystalType::BiBO_1),
      "LiIO3_1" => Ok(CrystalType::LiIO3_1),
      "LiIO3_2" => Ok(CrystalType::LiIO3_2),
      "LiNbO3_1" => Ok(CrystalType::LiNbO3_1),
      "LiNb_MgO" => Ok(CrystalType::LiNb_MgO),
      "KDP_1" => Ok(CrystalType::KDP_1),
      "AgGaS2_1" => Ok(CrystalType::AgGaS2_1),
      "AgGaSe2_1" => Ok(CrystalType::AgGaSe2_1),
      "AgGaSe2_2" => Ok(CrystalType::AgGaSe2_2),
      _ => {
        // we replace all = with : since we're using hjson as a hack
        let mut s = id.replace('=', ":");
        if !s.trim().starts_with('{') {
          s = format!("{{{}}}", s);
        }
        Ok(CrystalType::Expr(
          deser_hjson::from_str(&s).map_err(|e| SPDCError(e.to_string()))?,
        ))
      }
    }
  }

  /// Get all crystal meta information
  pub fn get_all_meta() -> Vec<CrystalMeta> {
    vec![
      bbo_1::META,
      ktp::META,
      bibo_1::META,
      linbo3_1::META,
      linb_mgo::META,
      kdp_1::META,
      aggase2_1::META,
      aggase2_2::META,
      liio3_2::META,
      liio3_1::LiIO3_1.get_meta(),
      aggas2_1::AgGaS2_1.get_meta(),
    ]
  }

  /// Get the crystal refraction indices for this crystal
  ///
  /// ## Example
  /// ```
  /// use spdcalc::{dim::ucum, na::Vector3, utils::*, CrystalType};
  /// let crystal = CrystalType::BBO_1;
  /// let nm = spdcalc::dim::f64prefixes::NANO * ucum::M;
  /// let indices = crystal.get_indices(720.0 * nm, from_celsius_to_kelvin(30.));
  /// let expected = ucum::Unitless::new(Vector3::new(
  ///   1.6631650519167869,
  ///   1.6631650519167869,
  ///   1.5463903834707935,
  /// ));
  /// assert_eq!(indices, expected)
  /// ```
  pub fn get_indices(&self, vacuum_wavelength: Wavelength, temperature: Kelvin<f64>) -> Indices {
    match &self {
      CrystalType::BBO_1 => bbo_1::get_indices(vacuum_wavelength, temperature),
      CrystalType::KTP => ktp::get_indices(vacuum_wavelength, temperature),
      CrystalType::LiNbO3_1 => linbo3_1::get_indices(vacuum_wavelength, temperature),
      CrystalType::LiNb_MgO => linb_mgo::get_indices(vacuum_wavelength, temperature),
      CrystalType::BiBO_1 => bibo_1::get_indices(vacuum_wavelength, temperature),
      CrystalType::KDP_1 => kdp_1::get_indices(vacuum_wavelength, temperature),
      CrystalType::AgGaSe2_1 => aggase2_1::get_indices(vacuum_wavelength, temperature),
      CrystalType::AgGaSe2_2 => aggase2_2::get_indices(vacuum_wavelength, temperature),
      CrystalType::LiIO3_2 => liio3_2::get_indices(vacuum_wavelength, temperature),
      CrystalType::LiIO3_1 => liio3_1::LiIO3_1.get_indices(vacuum_wavelength, temperature),
      CrystalType::AgGaS2_1 => aggas2_1::AgGaS2_1.get_indices(vacuum_wavelength, temperature),
      // CrystalType::Sellmeier(crystal) => crystal.get_indices(vacuum_wavelength, temperature),
      CrystalType::Expr(expr) => {
        let mut ctx = meval::Context::new();
        ctx.var("T", *((temperature - from_celsius_to_kelvin(20.0)) / K));
        match expr {
          CrystalExpr::Uniaxial { no, ne } => {
            let no = no.clone().bind_with_context(ctx.clone(), "l").unwrap()(
              *(vacuum_wavelength / (MICRO * M)),
            );
            let ne =
              ne.clone().bind_with_context(ctx, "l").unwrap()(*(vacuum_wavelength / (MICRO * M)));
            Indices::new(na::Vector3::new(no, no, ne))
          }
          CrystalExpr::Biaxial { nx, ny, nz } => {
            let nx = nx.clone().bind_with_context(ctx.clone(), "l").unwrap()(
              *(vacuum_wavelength / (MICRO * M)),
            );
            let ny = ny.clone().bind_with_context(ctx.clone(), "l").unwrap()(
              *(vacuum_wavelength / (MICRO * M)),
            );
            let nz =
              nz.clone().bind_with_context(ctx, "l").unwrap()(*(vacuum_wavelength / (MICRO * M)));
            Indices::new(na::Vector3::new(nx, ny, nz))
          }
        }
      }

      CrystalType::Interpolated(interpolated) => {
        use dim::f64prefixes::NANO;
        // Convert wavelength from meters to nanometers
        let wavelength_nm = *(vacuum_wavelength / (NANO * M));

        // Temperature parameter is ignored for interpolated crystals
        let (no, ne) = interpolated.get_indices(wavelength_nm).unwrap_or_else(|_| {
          // Fallback to reasonable default if interpolation fails
          (1.5, 1.5)
        });

        Indices::new(na::Vector3::new(no, no, ne))
      }
    }
  }

  /// Get the crystal meta information for specified crystal type
  pub fn get_meta(&self) -> CrystalMeta {
    match &self {
      CrystalType::BBO_1 => bbo_1::META,
      CrystalType::KTP => ktp::META,
      CrystalType::BiBO_1 => bibo_1::META,
      CrystalType::LiNbO3_1 => linbo3_1::META,
      CrystalType::LiNb_MgO => linb_mgo::META,
      CrystalType::KDP_1 => kdp_1::META,
      CrystalType::AgGaSe2_1 => aggase2_1::META,
      CrystalType::AgGaSe2_2 => aggase2_2::META,
      CrystalType::LiIO3_2 => liio3_2::META,
      CrystalType::LiIO3_1 => liio3_1::LiIO3_1.get_meta(),
      CrystalType::AgGaS2_1 => aggas2_1::AgGaS2_1.get_meta(),
      // CrystalType::Sellmeier(crystal) => crystal.get_meta(),
      CrystalType::Expr(_) => CrystalMeta {
        id: "Expr",
        name: "Expr",
        reference_url: "Expr",
        axis_type: OpticAxisType::PositiveUniaxial,
        point_group: PointGroup::HM_mm2,
        transmission_range: None,
        temperature_dependence_known: false,
      },

      CrystalType::Interpolated(interpolated) => interpolated.get_meta(),
    }
  }
}

#[cfg(feature = "pyo3")]
mod pyo3_impls {
  use super::*;
  use pyo3::{exceptions::PyValueError, prelude::*};

  impl FromPyObject<'_> for CrystalType {
    fn extract_bound(ob: &Bound<'_, pyo3::PyAny>) -> PyResult<Self> {
      let s: &str = ob.extract()?;
      CrystalType::from_str(s).map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
    }
  }

  impl ToPyObject for CrystalType {
    fn to_object(&self, py: Python<'_>) -> PyObject {
      self.to_string().to_object(py)
    }
  }

  impl IntoPy<PyObject> for CrystalType {
    fn into_py(self, py: Python<'_>) -> PyObject {
      self.to_string().into_py(py)
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use dim::f64prefixes::NANO;
  use dim::ucum::M;

  #[test]
  fn test_crystal_ids() {
    for meta in CrystalType::get_all_meta() {
      let crystal = CrystalType::from_string(meta.id).unwrap();
      assert_eq!(meta, crystal.get_meta());
    }
  }

  #[test]
  fn test_crystal_expr() {
    // BBO
    // no2=2.7359+0.01878/(λ2-0.01822)-0.01354λ2
    // ne2=2.3753+0.01224/(λ2-0.01667)-0.01516λ2
    // dno/dT = -9.3 x 10-6/°C
    // dne/dT = -16.6 x 10-6/°C
    let expr = r#"{
      "no": "sqrt(2.7359+0.01878/(l^2-0.01822)-0.01354*l^2) - 9.3e-6 * T",
      "ne": "sqrt(2.3753+0.01224/(l^2-0.01667)-0.01516*l^2) - 16.6e-6 * T"
    }"#;

    let crystal: CrystalType = serde_json::from_str(expr).unwrap();
    let other = CrystalType::from_string(
      r#"
      no = sqrt(2.7359+0.01878/(l^2-0.01822)-0.01354*l^2) - 9.3e-6 * T
      ne = sqrt(2.3753+0.01224/(l^2-0.01667)-0.01516*l^2) - 16.6e-6 * T
    "#,
    )
    .unwrap();

    assert_eq!(crystal, other);

    let indices = crystal.get_indices(1064.0 * NANO * M, from_celsius_to_kelvin(45.0));
    let expected = CrystalType::BBO_1.get_indices(1064.0 * NANO * M, from_celsius_to_kelvin(45.0));
    assert_eq!(indices, expected);
  }

  #[test]
  fn test_interpolated_crystal_deserialization() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [400.0, 500.0, 600.0],
      "no": [1.66, 1.65, 1.64],
      "ne": [1.55, 1.54, 1.53]
    }"#;

    let crystal: CrystalType = serde_json::from_str(json).unwrap();

    // Verify it's the Interpolated variant
    match crystal {
      CrystalType::Interpolated(_) => {},
      _ => panic!("Expected Interpolated variant"),
    }

    // Test get_indices with exact match
    let indices = crystal.get_indices(500.0 * NANO * M, from_celsius_to_kelvin(20.0));
    assert_eq!((*indices).x, 1.65);
    assert_eq!((*indices).y, 1.65);
    assert_eq!((*indices).z, 1.54);
  }

  #[test]
  fn test_interpolated_crystal_interpolation() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [400.0, 600.0],
      "no": [1.66, 1.64],
      "ne": [1.55, 1.53]
    }"#;

    let crystal: CrystalType = serde_json::from_str(json).unwrap();

    // Test interpolation at midpoint (500 nm)
    let indices = crystal.get_indices(500.0 * NANO * M, from_celsius_to_kelvin(20.0));

    use float_cmp::approx_eq;
    assert!(approx_eq!(f64, (*indices).x, 1.65, epsilon = 1e-10));
    assert!(approx_eq!(f64, (*indices).y, 1.65, epsilon = 1e-10));
    assert!(approx_eq!(f64, (*indices).z, 1.54, epsilon = 1e-10));
  }

  #[test]
  fn test_interpolated_crystal_extrapolation() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [500.0, 600.0],
      "no": [1.65, 1.64],
      "ne": [1.54, 1.53]
    }"#;

    let crystal: CrystalType = serde_json::from_str(json).unwrap();

    // Test below range
    let indices = crystal.get_indices(400.0 * NANO * M, from_celsius_to_kelvin(20.0));
    assert_eq!((*indices).x, 1.65);
    assert_eq!((*indices).z, 1.54);

    // Test above range
    let indices = crystal.get_indices(700.0 * NANO * M, from_celsius_to_kelvin(20.0));
    assert_eq!((*indices).x, 1.64);
    assert_eq!((*indices).z, 1.53);
  }

  #[test]
  fn test_interpolated_crystal_roundtrip() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [400.0, 500.0, 600.0],
      "no": [1.66, 1.65, 1.64],
      "ne": [1.55, 1.54, 1.53]
    }"#;

    let crystal: CrystalType = serde_json::from_str(json).unwrap();
    let serialized = serde_json::to_string(&crystal).unwrap();
    let deserialized: CrystalType = serde_json::from_str(&serialized).unwrap();

    assert_eq!(crystal, deserialized);
  }

  #[test]
  fn test_interpolated_crystal_get_meta() {
    let json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [400.0, 500.0],
      "no": [1.66, 1.65],
      "ne": [1.55, 1.54]
    }"#;

    let crystal: CrystalType = serde_json::from_str(json).unwrap();
    let meta = crystal.get_meta();

    assert_eq!(meta.id, "InterpolatedUniaxial");
    assert_eq!(meta.temperature_dependence_known, false);
    assert_eq!(meta.transmission_range, None);
  }

  #[test]
  fn test_interpolated_vs_expr() {
    // Create an Expr-based crystal with linear relationship
    let expr_json = r#"{
      "no": "1.66 - 0.0001 * (l - 0.4)",
      "ne": "1.55 - 0.0001 * (l - 0.4)"
    }"#;
    let expr_crystal: CrystalType = serde_json::from_str(expr_json).unwrap();

    // Create equivalent interpolated crystal
    let interp_json = r#"{
      "name": "InterpolatedUniaxial",
      "wavelengths_nm": [400.0, 500.0, 600.0],
      "no": [1.66, 1.65, 1.64],
      "ne": [1.55, 1.54, 1.53]
    }"#;
    let interp_crystal: CrystalType = serde_json::from_str(interp_json).unwrap();

    // Compare at test wavelength
    let nm = NANO * M;
    let expr_indices = expr_crystal.get_indices(500.0 * nm, from_celsius_to_kelvin(20.0));
    let interp_indices = interp_crystal.get_indices(500.0 * nm, from_celsius_to_kelvin(20.0));

    use float_cmp::approx_eq;
    assert!(approx_eq!(
      f64,
      (*expr_indices).x,
      (*interp_indices).x,
      epsilon = 1e-6
    ));
    assert!(approx_eq!(
      f64,
      (*expr_indices).z,
      (*interp_indices).z,
      epsilon = 1e-6
    ));
  }
}
