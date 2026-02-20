use super::*;
use crate::beam::*;
use crate::math::nelder_mead_1d;
use dim::ucum::*;
use na::{Rotation3, Vector3};
use std::f64::consts::FRAC_PI_2;

/// Crystal setup
///
/// This struct contains all the information about the crystal.
#[derive(Debug, Clone, PartialEq)]
pub struct CrystalSetup {
  /// The type of the crystal that influences the refractive indices
  pub crystal: CrystalType,
  /// The phasematching type through the crystal
  pub pm_type: PMType,
  /// The polar angle of the crystal optic axis
  pub phi: Angle,
  /// The azimuthal angle of the crystal optic axis
  pub theta: Angle,
  /// The length of the crystal
  pub length: Meter<f64>,
  /// The temperature of the crystal
  pub temperature: Kelvin<f64>,
  /// Whether the signal beam is traveling in the opposite z direction as the idler
  pub counter_propagation: bool,
}

impl AsRef<CrystalSetup> for CrystalSetup {
  fn as_ref(&self) -> &Self {
    self
  }
}

impl CrystalSetup {
  /// Convert a direction relative to the pump beam to a direction relative to the crystal optic axes
  pub fn to_crystal_frame(&self, direction: Direction) -> Direction {
    let crystal_rotation = Rotation3::from_euler_angles(0., *(self.theta / RAD), *(self.phi / RAD));
    crystal_rotation * direction
  }

  /// get the refractive index along specified direction, with wavelength and polarization type
  pub fn index_along(
    &self,
    vacuum_wavelength: Wavelength,
    direction: Direction,
    polarization: PolarizationType,
  ) -> RIndex {
    // Calculation follows https://physics.nist.gov/Divisions/Div844/publications/migdall/phasematch.pdf
    let indices = *self
      .crystal
      .get_indices(vacuum_wavelength, self.temperature);
    let n_inv2 = indices.map(|i| i.powi(-2));
    let s = self.to_crystal_frame(direction);
    let s_squared = s.map(|i| i * i);

    let sum_recip = Vector3::new(
      n_inv2.y + n_inv2.z,
      n_inv2.x + n_inv2.z,
      n_inv2.x + n_inv2.y,
    );
    let prod_recip = Vector3::new(
      n_inv2.y * n_inv2.z,
      n_inv2.x * n_inv2.z,
      n_inv2.x * n_inv2.y,
    );

    // Equation (11)
    // x² - bx + c = 0
    let b = -s_squared.dot(&sum_recip);
    let c = s_squared.dot(&prod_recip);

    let invxsq = match roots::find_roots_quadratic(1., b, c) {
      roots::Roots::One([x]) => -x,
      // x1 is < x2
      roots::Roots::Two([x1, x2]) => {
        match polarization {
          // fast
          PolarizationType::Ordinary => x2,
          // slow
          PolarizationType::Extraordinary => x1,
        }
      }
      _ => return RIndex::new(0.), // imaginary index
    };

    if invxsq < 0. {
      RIndex::new(0.)
    } else {
      RIndex::new(1. / invxsq.sqrt())
    }
  }

  /// calculate the optimal crystal theta
  /// by minimizing delta k
  pub fn optimum_theta(&self, signal: &SignalBeam, pump: &PumpBeam) -> Angle {
    let theta_s_e = signal.theta_external(self);

    let delta_k = move |theta| {
      let mut crystal_setup = self.clone();
      let mut signal = signal.clone();

      crystal_setup.theta = theta * RAD;
      signal.set_theta_external(theta_s_e, &crystal_setup);

      let idler =
        IdlerBeam::try_new_optimum(&signal, pump, &crystal_setup, PeriodicPoling::Off).unwrap();
      let del_k = delta_k(
        signal.frequency(),
        idler.frequency(),
        &signal,
        &idler,
        pump,
        &crystal_setup,
        PeriodicPoling::Off,
      );

      (del_k * M / RAD).z.abs()
      // let (ws, wi) = (signal.frequency(), idler.frequency());
      // let zs = self.optimal_waist_position(signal.vacuum_wavelength(), signal.polarization());
      // let zi = self.optimal_waist_position(idler.vacuum_wavelength(), idler.polarization());

      // let spdc = SPDC::new(
      //   crystal_setup,
      //   signal,
      //   idler,
      //   pump.clone(),
      //   1. * M,
      //   1. * dim::ucum::W,
      //   0.1,
      //   PeriodicPoling::Off,
      //   zs,
      //   zi,
      //   MetersPerMilliVolt::new(1.)
      // );
      // - crate::jsa::jsa_raw(ws, wi, &spdc, None).norm_sqr()
    };

    let guess = PI / 6.;
    let theta = nelder_mead_1d(delta_k, (guess, guess + 1.), 1000, 0., FRAC_PI_2, 1e-6);

    theta * RAD
  }

  /// Assign the optimum crystal theta
  pub fn assign_optimum_theta(&mut self, signal: &SignalBeam, pump: &PumpBeam) {
    self.theta = self.optimum_theta(signal, pump);
  }

  /// Calculate the optimal waist position inside the crystal.
  ///
  /// The position is the distance from the crystal exit surface along the z-axis.
  /// z_{s,i} = -\frac{1}{2}\frac{L}{n_z(\lambda_{s,i})}
  pub fn optimal_waist_position(
    &self,
    wavelength: Wavelength,
    polarization: PolarizationType,
  ) -> Distance {
    -0.5 * self.length
      / self.index_along(
        wavelength,
        na::Unit::new_normalize(na::Vector3::z()),
        polarization,
      )
  }
}

#[cfg(test)]
mod test {
  use super::*;
  use crate::utils::testing::*;

  #[test]
  fn index_along_test() {
    let mut spdc = SPDC::default();
    spdc.crystal_setup.phi = Angle::new(0.);
    spdc.crystal_setup.theta = Angle::new(PI / 2.);
    spdc.crystal_setup.crystal = CrystalType::BBO_1;
    spdc.signal.set_angles(0. * DEG, 53. * DEG);
    let n = spdc.crystal_setup.index_along(
      spdc.signal.vacuum_wavelength(),
      spdc.signal.direction(),
      spdc.signal.polarization(),
    );

    assert_eq!(n, Unitless::new(1.6017685463810718));
  }

  #[test]
  fn index_along_angle_test() {
    let mut spdc = SPDC::default();
    spdc.crystal_setup.phi = Angle::new(0.);
    spdc.crystal_setup.theta = Angle::new(0.1);
    spdc.crystal_setup.crystal = CrystalType::KTP;
    spdc.signal.set_angles(0. * DEG, 5. * DEG);
    let no = spdc.crystal_setup.index_along(
      spdc.signal.vacuum_wavelength(),
      spdc.signal.direction(),
      PolarizationType::Ordinary,
    );
    let ne = spdc.crystal_setup.index_along(
      spdc.signal.vacuum_wavelength(),
      spdc.signal.direction(),
      PolarizationType::Extraordinary,
    );

    let actual = (*no, *ne);
    let expected = (1.7340844750789677, 1.7319709938216157);

    assert_nearly_equal!("index_along_angle no", actual.0, expected.0, 1e-12);
    assert_nearly_equal!("index_along_angle ne", actual.1, expected.1, 1e-12);
  }

  #[test]
  fn index_along_expr_crystal_test() {
    // Crystal with constant indices: no = 1, ne = 2
    let crystal = CrystalType::from_string(r#"no = 1, ne = 2"#).unwrap();
    let crystal_setup = CrystalSetup {
      crystal,
      pm_type: PMType::Type2_e_eo,
      phi: Angle::new(0.),
      theta: 90. * DEG,
      length: 1e-3 * M,
      temperature: crate::utils::from_celsius_to_kelvin(20.),
      counter_propagation: false,
    };

    let wavelength = 1550e-9 * M;
    // Direction of propagation along z
    let direction = Beam::new(
      PolarizationType::Extraordinary,
      Angle::new(0.),
      Angle::new(0.),
      wavelength,
      1. * M,
    )
    .direction();

    // For e->eo: signal is Extraordinary, idler is Ordinary
    let signal_pol = crystal_setup.pm_type.signal_polarization();
    let idler_pol = crystal_setup.pm_type.idler_polarization();

    let n_signal = crystal_setup.index_along(wavelength, direction, signal_pol);
    let n_idler = crystal_setup.index_along(wavelength, direction, idler_pol);

    assert_eq!(signal_pol, PolarizationType::Extraordinary);
    assert_eq!(idler_pol, PolarizationType::Ordinary);
    dbg!("n_signal: {}, n_idler: {}", *n_signal, *n_idler);
    assert_eq!(
      *n_signal, 2.0,
      "Expected signal index to be 2.0, got {}",
      *n_signal
    );
    assert_eq!(
      *n_idler, 1.0,
      "Expected idler index to be 1.0, got {}",
      *n_idler
    );
  }

  #[test]
  fn optimum_theta_test() {
    let mut spdc = testing_props(false);
    spdc.assign_optimum_crystal_theta();
    assert_nearly_equal!(
      "crystal_theta",
      spdc.crystal_setup.theta.value_unsafe,
      0.5123773602350467,
      0.001
    );
  }
}
