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
}
