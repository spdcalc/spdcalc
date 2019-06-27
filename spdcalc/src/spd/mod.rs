use crate::*;
use dim::ucum;
use photon::Photon;
use crystal::CrystalSetup;
use dim::f64prefixes::MILLI;

pub struct PeriodicPolling {
  pub period : ucum::Meter<f64>,
  pub sign: Sign,
}

pub fn calc_delta_k(
  signal :&Photon,
  idler :&Photon,
  pump :&Photon,
  crystal_setup :&CrystalSetup,
  pp: Option<PeriodicPolling>
) -> Momentum3 {

  let r_s = signal.get_direction();
  let r_i = idler.get_direction();

  let ns_by_ls = *(ucum::M * signal.get_index(&crystal_setup) / signal.get_wavelength());
  let ni_by_li = *(ucum::M * idler.get_index(&crystal_setup) / idler.get_wavelength());
  let np_by_lp = *(ucum::M * pump.get_index(&crystal_setup) / pump.get_wavelength());

  let mut dk =
      r_s.as_ref() * ns_by_ls
    + r_i.as_ref() * ni_by_li;

  dk.z = np_by_lp - dk.z;

  // put into milliJoule seconds
  (PI2 / MILLI) * match pp {
    Some(poling) => {
      dk.z -= 1.0 / (poling.sign * (*(poling.period/ucum::M)));
      Momentum3::new(dk)
    },
    None => Momentum3::new(dk),
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::crystal::CrystalSetup;
  extern crate float_cmp;
  use float_cmp::*;
  use crate::utils::*;
  use photon::PhotonType;
  use ucum::*;
  use dim::f64prefixes::*;

  fn init() -> (CrystalSetup, Photon, Photon, Photon) {
    let wavelength = 1550. * NANO * M;
    let waist = WaistSize::new(100.0 * MICRO * M, 100.0 * MICRO * M);
    let crystal_setup = CrystalSetup{
      crystal: Crystal::BBO_1,
      pm_type : crystal::PMType::Type2_e_eo,
      theta : -3.0 * DEG,
      phi : 1.0 * DEG,
      length : 2_000.0 * MICRO * M,
      temperature : from_celsius_to_kelvin(20.0),
    };

    let signal = Photon::new(PhotonType::Signal, 2.0 * DEG, 3.0 * DEG, wavelength, waist);
    let idler = Photon::new(PhotonType::Idler, 182. * DEG, 0.05132276669556576 * RAD, wavelength, waist);
    let pump = Photon::new(PhotonType::Pump, 0. * DEG, 0. * DEG, 775. * NANO * M, waist);

    (crystal_setup, signal, idler, pump)
  }

  #[test]
  fn calc_delta_k_test() {
    let (crystal_setup, signal, idler, pump) = init();
    let expected = na::Vector3::new(6908.816094920548, 241.26117431161828, -358.7018265926745);
    let pp = PeriodicPolling{
      period: 0.00004656366863331685 * ucum::M,
      sign: Sign::POSITIVE,
    };

    let del_k = *(calc_delta_k(&signal, &idler, &pump, &crystal_setup, Option::Some(pp)) / ucum::J / ucum::S);
    assert!(approx_eq!(f64, del_k.x, expected.x, ulps = 2, epsilon = 1e-6), "actual: {}, expected: {}", del_k.x, expected.x);
    assert!(approx_eq!(f64, del_k.y, expected.y, ulps = 2, epsilon = 1e-6), "actual: {}, expected: {}", del_k.y, expected.y);
    assert!(approx_eq!(f64, del_k.z, expected.z, ulps = 2, epsilon = 1e-6), "actual: {}, expected: {}", del_k.z, expected.z);
  }
}
