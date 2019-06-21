use super::*;
use photon::PhotonType;
use dim::ucum::*;
use na::{Vector3, Rotation3};

/// The phasematch type
#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone)]
pub enum PMType {
  /// Type 0:   o -> o + o
  Type0_o_oo,
  /// Type 0:   e -> e + e
  Type0_e_ee,
  /// Type 1:   e -> o + o
  Type1_e_oo,

  /// Type 2:   e -> e + o
  Type2_e_eo,
  /// Type 2:   e -> o + e
  Type2_e_oe,
}

#[derive(Debug, Copy, Clone)]
pub struct CrystalSetup {
  pub crystal :Crystal,
  pub pm_type : PMType,
  pub theta :Angle,
  pub phi :Angle,
  pub length : Meter<f64>,
  pub temperature : Kelvin<f64>,
}

// internal helper to solve Equation (11)
// sign = (slow: 1, fast: -1)
fn solve_for_n(b :f64, c :f64, sign : i16) -> f64 {
  let d = b * b - 4.0 * c;
  (2.0 / (b + (sign as f64) * d.sqrt())).sqrt()
}

impl CrystalSetup {
  pub fn get_index_along(&self, wavelength : Wavelength, direction : Direction, photon_type : &PhotonType) -> RIndex {
    // Calculation follows https://physics.nist.gov/Divisions/Div844/publications/migdall/phasematch.pdf
    let indices = self.crystal.get_indices(wavelength, self.temperature);
    let n_inv2 = indices.map(|i| i.value_unsafe.powi(-2));
    let crystal_rotation = Rotation3::from_euler_angles(0., self.theta.value_unsafe, self.phi.value_unsafe);
    let s = crystal_rotation * direction;
    let s_squared = s.map(|i| i * i);

    let sum_recip = Vector3::new(
      n_inv2.y + n_inv2.z,
      n_inv2.x + n_inv2.z,
      n_inv2.x + n_inv2.y
    );
    let prod_recip = Vector3::new(
      n_inv2.y * n_inv2.z,
      n_inv2.x * n_inv2.z,
      n_inv2.x * n_inv2.y
    );

    // Equation (11)
    // x² + bx + c = 0
    let b = s_squared.dot(&sum_recip);
    let c = s_squared.dot(&prod_recip);

    let slow = 1;
    let fast = -1;

    dim::ucum::ONE * match &self.pm_type {
      PMType::Type0_o_oo =>
        solve_for_n(b, c, fast),
      PMType::Type0_e_ee =>
        solve_for_n(b, c, slow),
      PMType::Type1_e_oo => {
        let sign = match photon_type { PhotonType::Pump => slow, _ => fast };
        solve_for_n(b, c, sign)
      },
      PMType::Type2_e_eo => {
        let sign = match photon_type { PhotonType::Idler => fast, _ => slow };
        solve_for_n(b, c, sign)
      },
      PMType::Type2_e_oe => {
        let sign = match photon_type { PhotonType::Signal => fast, _ => slow };
        solve_for_n(b, c, sign)
      },
    }
  }
}