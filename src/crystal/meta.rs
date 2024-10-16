use strum_macros::EnumString;

/// Type of Optic Axis the crystal has
#[derive(
  Debug, Serialize, PartialEq, Deserialize, Copy, Clone, EnumString, strum_macros::Display,
)]
pub enum OpticAxisType {
  /// Uniaxial crystal with positive birifringence
  PositiveUniaxial,
  /// Uniaxial crystal with negative birifringence
  NegativeUniaxial,
  /// Biaxial crystal with positive birifringence
  PositiveBiaxial,
  /// Biaxial crystal with negative birifringence
  NegativeBiaxial,
}

/// Crystallographic Point Group
///
/// Represented with Hermann-Mauguin (Universal) short notation
///
/// **Note**: Due to the character limitations of rust enum names,
/// `s` represents `/` and `i` represents rotary-inversion
///
/// [General Information](https://en.wikipedia.org/wiki/Crystallographic_point_group)
#[allow(non_camel_case_types)]
#[derive(
  Debug, Serialize, PartialEq, Deserialize, Copy, Clone, EnumString, strum_macros::Display,
)]
pub enum PointGroup {
  /// Triclinic 1
  HM_1,
  /// Triclinic -1
  HM_i1,

  /// Monoclinic 2
  HM_2,
  /// Monoclinic m
  HM_m,
  /// Monoclinic 2/m
  HM_2sm,

  /// Orthorhombic 222
  HM_222,
  /// Orthorhombic mm2
  HM_mm2,
  /// Orthorhombic mmm
  HM_mmm,

  /// Tetragonal 4
  HM_4,
  /// Tetragonal -4
  HM_i4,
  /// Tetragonal 4/m
  HM_4sm,
  /// Tetragonal 422
  HM_422,
  /// Tetragonal 4mm
  HM_4mm,
  /// Tetragonal -42m
  HM_i42m,
  /// Tetragonal 4/mmm
  HM_4smmm,

  /// Trigonal 3
  HM_3,
  /// Trigonal -3
  HM_i3,
  /// Trigonal 32
  HM_32,
  /// Trigonal 3m
  HM_3m,
  /// Trigonal -3m
  HM_i3m,

  /// Hexagonal 6
  HM_6,
  /// Hexagonal -6
  HM_i6,
  /// Hexagonal 6/m
  HM_6sm,
  /// Hexagonal 622
  HM_622,
  /// Hexagonal 6mm
  HM_6mm,
  /// Hexagonal -62m
  HM_i62m,
  /// Hexagonal 6/mmm
  HM_6smmm,

  /// Cubic 23
  HM_23,
  /// Cubic m-3
  HM_mi3,
  /// Cubic 432
  HM_432,
  /// Cubic -43m
  HM_i43m,
  /// Cubic m-3m
  HM_mi3m,
}

/// A valid wavelength range for a crystal
#[derive(Debug, Serialize, PartialEq, Deserialize, Copy, Clone)]
pub struct ValidWavelengthRange(pub f64, pub f64);

/// Meta information about the crystal
#[derive(Debug, Serialize, PartialEq, Deserialize, Copy, Clone)]
pub struct CrystalMeta {
  /// the internal enum name
  pub id: &'static str,
  /// The name of the crystal
  pub name: &'static str,
  /// A url to a reference for the crystal
  pub reference_url: &'static str,
  /// Optic axis type of the crystal
  pub axis_type: OpticAxisType,
  /// Point Group (class) of the crystal
  pub point_group: PointGroup,
  /// The wavelength range within which the crystal transmits light
  pub transmission_range: Option<ValidWavelengthRange>,
  /// Whether or not temperature dependence is known
  pub temperature_dependence_known: bool,
}

#[cfg(feature = "pyo3")]
mod pyo3_impls {
  use super::*;
  use pyo3::{
    prelude::*,
    types::{PyDict, PyTuple},
  };

  impl ToPyObject for ValidWavelengthRange {
    fn to_object(&self, py: Python<'_>) -> PyObject {
      let tuple = PyTuple::new_bound(py, &[self.0, self.1]);
      tuple.into()
    }
  }

  impl ToPyObject for CrystalMeta {
    fn to_object(&self, py: Python<'_>) -> PyObject {
      let dict = PyDict::new_bound(py);
      dict.set_item("id", self.id).unwrap();
      dict.set_item("name", self.name).unwrap();
      dict.set_item("reference_url", self.reference_url).unwrap();
      dict
        .set_item("axis_type", self.axis_type.to_string())
        .unwrap();
      dict
        .set_item("point_group", self.point_group.to_string())
        .unwrap();
      dict
        .set_item("transmission_range", self.transmission_range)
        .unwrap();
      dict
        .set_item(
          "temperature_dependence_known",
          self.temperature_dependence_known,
        )
        .unwrap();

      dict.into()
    }
  }

  impl IntoPy<PyObject> for CrystalMeta {
    fn into_py(self, py: Python<'_>) -> PyObject {
      self.to_object(py)
    }
  }
}
