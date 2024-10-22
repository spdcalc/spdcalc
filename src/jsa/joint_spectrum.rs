use crate::{
  math::Integrator, phasematch::*, Complex, Frequency, FrequencySpace, IntoSignalIdlerIterator,
  JSIUnits, JsiNorm, JsiSinglesNorm, PerMeter3, PerMeter4, SPDCError, SPDC,
};
use rayon::prelude::*;

// This defines a more than reasonable box around frequency ranges to...
// 1. speed up calculations
// 2. avoid non-sensical values
fn invalid_frequencies(omega_s: Frequency, omega_i: Frequency, spdc: &SPDC) -> bool {
  let omega_p = spdc.pump.frequency();
  use crate::dim::Abs;
  omega_s <= Frequency::new(0.)
    || omega_i <= Frequency::new(0.)
    || omega_s > omega_p
    || omega_i > omega_p
    || (omega_s - omega_i).abs() > 0.75 * omega_p
}

/// The raw joint spectrum amplitude
///
/// This is the JSA for coincidences that does not include count rate constants.
pub fn jsa_raw(
  omega_s: Frequency,
  omega_i: Frequency,
  spdc: &SPDC,
  integrator: Integrator,
) -> Complex<f64> {
  if invalid_frequencies(omega_s, omega_i, spdc) {
    return Complex::new(0., 0.);
  }
  let alpha = pump_spectral_amplitude(omega_s + omega_i, spdc);
  // check the threshold
  if alpha < spdc.pump_spectrum_threshold {
    Complex::new(0., 0.)
  } else {
    let f = phasematch_fiber_coupling(omega_s, omega_i, spdc, integrator) / PerMeter4::new(1.);
    *(alpha * f)
  }
}

/// The raw joint spectrum intensity
///
/// This is the JSI for singles that does not include count rate constants.
pub fn jsi_singles_raw(
  omega_s: Frequency,
  omega_i: Frequency,
  spdc: &SPDC,
  integrator: Integrator,
) -> f64 {
  if invalid_frequencies(omega_s, omega_i, spdc) {
    return 0.;
  }
  let alpha = pump_spectral_amplitude(omega_s + omega_i, spdc);
  // check the threshold
  if alpha < spdc.pump_spectrum_threshold {
    0.
  } else {
    let fs =
      phasematch_singles_fiber_coupling(omega_s, omega_i, spdc, integrator) / PerMeter3::new(1.);
    // use crate::utils::frequency_to_vacuum_wavelength;
    // let mut setup : SPDCSetup = spdc.clone().into();
    // setup.signal.set_wavelength(frequency_to_vacuum_wavelength(omega_s));
    // setup.idler.set_wavelength(frequency_to_vacuum_wavelength(omega_s));
    // let fs = calc_singles_phasematch_fiber_coupling(&setup).0;
    *(alpha.powi(2) * fs)
  }
}

/// Joint Spectrum Calculation Helper
///
/// This is the primary way of calculating aspects of the joint spectrum for a given setup.
/// It provides some optimization by caching the JSA and JSI at the center frequency used to calculate normalized values.
#[derive(Clone, Debug)]
pub struct JointSpectrum {
  spdc: SPDC,
  integrator: Integrator,
  jsa_center: f64,
  jsi_singles_center: f64,
}

impl JointSpectrum {
  /// Create a new instance
  pub fn new(spdc: SPDC, integrator: Integrator) -> Self {
    let spdc_optimal = spdc.clone().try_as_optimum().unwrap();
    let jsa_center_norm = jsi_normalization(
      spdc_optimal.signal.frequency(),
      spdc_optimal.idler.frequency(),
      &spdc_optimal,
    ) / JsiNorm::new(1.);
    let jsa_center = jsa_center_norm.sqrt()
      * jsa_raw(
        spdc_optimal.signal.frequency(),
        spdc_optimal.idler.frequency(),
        &spdc_optimal,
        integrator,
      )
      .norm();
    let jsi_singles_center_norm = *(jsi_singles_normalization(
      spdc_optimal.signal.frequency(),
      spdc_optimal.idler.frequency(),
      &spdc_optimal,
    ) / JsiSinglesNorm::new(1.));
    let jsi_singles_center = jsi_singles_center_norm
      * jsi_singles_raw(
        spdc_optimal.signal.frequency(),
        spdc_optimal.idler.frequency(),
        &spdc_optimal,
        integrator,
      );

    Self {
      spdc,
      integrator,
      jsa_center,
      jsi_singles_center,
    }
  }

  /// Get the value of the JSA at specified signal/idler frequencies
  ///
  /// Technically the units should be 1/sqrt(s)/(rad/s)
  pub fn jsa(&self, omega_s: Frequency, omega_i: Frequency) -> Complex<f64> {
    let jsa = jsa_raw(omega_s, omega_i, &self.spdc, self.integrator);
    use num::Zero;
    if jsa == Complex::zero() {
      Complex::zero()
    } else {
      let n = jsi_normalization(omega_s, omega_i, &self.spdc) / JsiNorm::new(1.);
      n.sqrt() * jsa
    }
  }

  /// Get the JSA over a specified range of signal/idler frequencies
  pub fn jsa_range<T: IntoSignalIdlerIterator>(&self, range: T) -> Vec<Complex<f64>> {
    range
      .into_signal_idler_par_iterator()
      .map(|(ws, wi)| self.jsa(ws, wi))
      .collect()
  }

  /// Get the normalized value of the JSA at specified signal/idler frequencies
  ///
  /// This is unitless and normalized to the optimal setup
  pub fn jsa_normalized(&self, omega_s: Frequency, omega_i: Frequency) -> Complex<f64> {
    self.jsa(omega_s, omega_i) / self.jsa_center
  }

  /// Get the normalized value of the JSA at specified signal/idler frequencies
  pub fn jsa_normalized_range<T: IntoSignalIdlerIterator>(&self, range: T) -> Vec<Complex<f64>> {
    range
      .into_signal_idler_par_iterator()
      .map(|(ws, wi)| self.jsa_normalized(ws, wi))
      .collect()
  }

  /// Get the value of the JSI at specified signal/idler frequencies
  ///
  /// Units are: per second per (rad/s)^2, which when integrated over
  /// gives a value proportional to the count rate (counts/s)
  pub fn jsi(&self, omega_s: Frequency, omega_i: Frequency) -> JSIUnits<f64> {
    let jsa = jsa_raw(omega_s, omega_i, &self.spdc, self.integrator);
    use num::Zero;
    if jsa == Complex::zero() {
      JSIUnits::new(0.)
    } else {
      let n = jsi_normalization(omega_s, omega_i, &self.spdc) / JsiNorm::new(1.);
      JSIUnits::new(*(n * jsa.norm_sqr()))
    }
  }

  /// Get the JSI over a specified range of signal/idler frequencies
  pub fn jsi_range<T: IntoSignalIdlerIterator>(&self, range: T) -> Vec<JSIUnits<f64>> {
    range
      .into_signal_idler_par_iterator()
      .map(|(ws, wi)| self.jsi(ws, wi))
      .collect()
  }

  /// Get the normalized value of the JSI at specified signal/idler frequencies
  ///
  /// This is unitless and normalized to the optimal setup
  pub fn jsi_normalized(&self, omega_s: Frequency, omega_i: Frequency) -> f64 {
    *(self.jsi(omega_s, omega_i) / JSIUnits::new(1.)) / self.jsa_center.powi(2)
  }

  /// Get the normalized value of the JSI at specified signal/idler frequencies
  pub fn jsi_normalized_range<T: IntoSignalIdlerIterator>(&self, range: T) -> Vec<f64> {
    range
      .into_signal_idler_par_iterator()
      .map(|(ws, wi)| self.jsi_normalized(ws, wi))
      .collect()
  }

  // /// Get the value of the JSA Singles at specified signal/idler frequencies
  // ///
  // /// Technically the units should be 1/sqrt(s)/(rad/s)
  // pub fn jsa_singles(&self, omega_s: Frequency, omega_i: Frequency) -> f64 {
  //   let jsi = jsi_singles_raw(omega_s, omega_i, &self.spdc, self.integrator);
  //   if jsi == 0. {
  //     0.
  //   } else {
  //     let n = jsi_singles_normalization(omega_s, omega_i, &self.spdc) / JsiSinglesNorm::new(1.);
  //     (n * jsi).sqrt()
  //   }
  // }

  // /// Get the JSA Singles over a specified range of signal/idler frequencies
  // pub fn jsa_singles_range<T: IntoSignalIdlerIterator>(&self, range: T) -> Vec<f64> {
  //   range
  //     .into_signal_idler_par_iterator()
  //     .map(|(ws, wi)| self.jsa_singles(ws, wi))
  //     .collect()
  // }

  // /// Get the normalized value of the JSA Singles at specified signal/idler frequencies
  // ///
  // /// This is unitless and normalized to the optimal setup
  // pub fn jsa_singles_normalized(&self, omega_s: Frequency, omega_i: Frequency) -> f64 {
  //   self.jsi_singles_normalized(omega_s, omega_i).sqrt()
  // }

  // /// Get the normalized value of the JSA Singles at specified signal/idler frequencies
  // pub fn jsa_singles_normalized_range<T: IntoSignalIdlerIterator>(&self, range: T) -> Vec<f64> {
  //   range
  //     .into_signal_idler_par_iterator()
  //     .map(|(ws, wi)| self.jsa_singles_normalized(ws, wi))
  //     .collect()
  // }

  /// Get the value of the JSI Singles at specified signal/idler frequencies
  ///
  /// Units are: per second per (rad/s)^2, which when integrated over
  /// gives a value proportional to the count rate (counts/s)
  pub fn jsi_singles(&self, omega_s: Frequency, omega_i: Frequency) -> JSIUnits<f64> {
    let jsi = jsi_singles_raw(omega_s, omega_i, &self.spdc, self.integrator);
    if jsi == 0. {
      JSIUnits::new(0.)
    } else {
      let n = jsi_singles_normalization(omega_s, omega_i, &self.spdc) / JsiSinglesNorm::new(1.);
      JSIUnits::new(*n * jsi)
    }
  }

  /// Get the JSI Singles over a specified range of signal/idler frequencies
  pub fn jsi_singles_range<T: IntoSignalIdlerIterator>(&self, range: T) -> Vec<JSIUnits<f64>> {
    range
      .into_signal_idler_par_iterator()
      .map(|(ws, wi)| self.jsi_singles(ws, wi))
      .collect()
  }

  /// Get the JSI Singles for the idler over a specified range of signal/idler frequencies
  pub fn jsi_singles_idler_range<T: IntoSignalIdlerIterator>(
    &self,
    range: T,
  ) -> Vec<JSIUnits<f64>> {
    let swapped = self.spdc.clone().with_swapped_signal_idler();
    let idler_spectrum = Self::new(swapped, self.integrator);
    range
      .into_signal_idler_par_iterator()
      .map(|(ws, wi)| idler_spectrum.jsi_singles(wi, ws))
      .collect()
  }

  /// Get the normalized value of the JSI Singles at specified signal/idler frequencies
  ///
  /// This is unitless and normalized to the optimal setup
  pub fn jsi_singles_normalized(&self, omega_s: Frequency, omega_i: Frequency) -> f64 {
    *(self.jsi_singles(omega_s, omega_i) / JSIUnits::new(1.)) / self.jsi_singles_center
  }

  /// Get the normalized value of the JSI Singles at specified signal/idler frequencies
  pub fn jsi_singles_normalized_range<T: IntoSignalIdlerIterator>(&self, range: T) -> Vec<f64> {
    range
      .into_signal_idler_par_iterator()
      .map(|(ws, wi)| self.jsi_singles_normalized(ws, wi))
      .collect()
  }

  /// Get the normalized value of the JSI Singles for the idler at specified signal/idler frequencies
  pub fn jsi_singles_idler_normalized_range<T: IntoSignalIdlerIterator>(
    &self,
    range: T,
  ) -> Vec<f64> {
    let swapped = self.spdc.clone().with_swapped_signal_idler();
    let idler_spectrum = Self::new(swapped, self.integrator);
    range
      .into_signal_idler_par_iterator()
      .map(|(ws, wi)| idler_spectrum.jsi_singles_normalized(wi, ws))
      .collect()
  }

  /// Calculate the schmidt number for this SPDC configuration over a specified range of signal/idler frequencies
  pub fn schmidt_number<R: Into<FrequencySpace>>(&self, range: R) -> Result<f64, SPDCError> {
    crate::math::schmidt_number(self.jsa_range(range.into()))
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{
    utils::{frequency_to_vacuum_wavelength, vacuum_wavelength_to_frequency, Steps},
    PeriodicPoling, SPDCConfig,
  };
  use dim::{f64prefixes::*, ucum::*};

  fn get_spdc() -> SPDC {
    let json = serde_json::json!({
      "crystal": {
        "kind": "KTP",
        "pm_type": "e->eo",
        "phi_deg": 0,
        "theta_deg": 90,
        "length_um": 14_000,
        "temperature_c": 20
      },
      "pump": {
        "wavelength_nm": 775,
        "waist_um": 200,
        "bandwidth_nm": 0.5,
        "average_power_mw": 300
      },
      "signal": {
        "wavelength_nm": 1550,
        "phi_deg": 0,
        "theta_external_deg": 0,
        "waist_um": 100,
        "waist_position_um": "auto"
      },
      "idler": "auto",
      "periodic_poling": {
        "poling_period_um": "auto"
      },
      "deff_pm_per_volt": 7.6
    });

    // let json = serde_json::json!({
    //   "crystal": {
    //     "kind": "KTP",
    //     "pm_type": "e->eo",
    //     "phi_deg": 0,
    //     "theta_deg": 0,
    //     "length_um": 30_000,
    //     "temperature_c": 20
    //   },
    //   "pump": {
    //     "wavelength_nm": 775,
    //     "waist_um": 50,
    //     "bandwidth_nm": 5.35,
    //     "average_power_mw": 500
    //   },
    //   "signal": {
    //     "wavelength_nm": 1550,
    //     "phi_deg": 0,
    //     "theta_external_deg": 0,
    //     "waist_um": 50,
    //     "waist_position_um": "auto"
    //   },
    //   "idler": "auto",
    //   "periodic_poling": {
    //     "poling_period_um": "auto"
    //   }
    // });

    let config: SPDCConfig = serde_json::from_value(json).expect("Could not unwrap json");
    
    config
      .try_as_spdc()
      .expect("Could not convert to SPDC instance")
  }

  #[test]
  fn test_efficiency() {
    let spdc = get_spdc();

    let spectrum = JointSpectrum::new(spdc.clone(), Integrator::default());
    let frequencies = FrequencySpace::new(
      (
        vacuum_wavelength_to_frequency(1541.54 * NANO * M),
        vacuum_wavelength_to_frequency(1558.46 * NANO * M),
        20,
      ),
      (
        vacuum_wavelength_to_frequency(1541.63 * NANO * M),
        vacuum_wavelength_to_frequency(1558.56 * NANO * M),
        20,
      ),
    );

    let jsi = spectrum.jsi_range(frequencies);
    let jsi_singles = spectrum.jsi_singles_range(frequencies);
    let jsi_singles_idler = spectrum.jsi_singles_idler_range(frequencies);

    let steps = frequencies.as_steps();
    let dxdy = Steps::from(steps.0).division_width() * Steps::from(steps.1).division_width();
    let coinc_rate: Hertz<f64> = jsi.into_iter().sum::<JSIUnits<f64>>() * dxdy;
    let singles_signal_rate: Hertz<f64> = jsi_singles.into_iter().sum::<JSIUnits<f64>>() * dxdy;
    let singles_idler_rate: Hertz<f64> =
      jsi_singles_idler.into_iter().sum::<JSIUnits<f64>>() * dxdy;

    assert!(coinc_rate < singles_signal_rate);
    assert!(coinc_rate < singles_idler_rate);
  }

  #[test]
  fn test_normalized_jsa() {
    let json = serde_json::json!({
      "crystal": {
        "kind": "KTP",
        "pm_type": "e oo",
        "phi_deg": 0.0,
        "theta_deg": "auto",
        "length_um": 20000.0,
        "temperature_c": 20.0
      },
      "pump": {
        "wavelength_nm": 775,
        "waist_um": 100.0,
        "bandwidth_nm": 5.35,
        "average_power_mw": 1.0,
        "spectrum_threshold": 0.01
      },
      "signal": {
        "wavelength_nm": 1550,
        "phi_deg": 0.0,
        "theta_deg": 0.0,
        "theta_external_deg": null,
        "waist_um": 100.0,
        "waist_position_um": -5766.731750218876
      },
      "idler": {
        "wavelength_nm": 1550,
        "phi_deg": 180.0,
        "theta_deg": 0.0,
        "theta_external_deg": null,
        "waist_um": 100.0,
        "waist_position_um": -5506.780644729153
      },
      // "periodic_poling": {
      //   "poling_period_um": 46.52032850062398,
      //   "apodization": null
      // },
      "deff_pm_per_volt": 1.0
    });
    let integrator = Integrator::default();
    let config: SPDCConfig = serde_json::from_value(json).expect("Could not unwrap json");
    let spdc = config
      .try_as_spdc()
      .expect("Could not convert to SPDC instance");
    let optimal = spdc.clone().try_as_optimum().unwrap();
    dbg!(&spdc);
    dbg!(&optimal);
    dbg!(spdc
      .joint_spectrum(integrator)
      .jsa(spdc.signal.frequency(), spdc.idler.frequency()));
    dbg!(optimal
      .joint_spectrum(integrator)
      .jsa(optimal.signal.frequency(), optimal.idler.frequency()));
    let sp = optimal.joint_spectrum(integrator);
    let jsa = sp.jsa_normalized(spdc.signal.frequency(), spdc.idler.frequency());
    dbg!(jsa.norm());
    // assert!(float_cmp::approx_eq!(f64, jsa.norm(), 1.0));
    dbg!(delta_k(
      optimal.signal.frequency(),
      optimal.idler.frequency(),
      &optimal.signal,
      &optimal.idler,
      &optimal.pump,
      &optimal.crystal_setup,
      PeriodicPoling::Off
    ));
    let range = optimal.optimum_range(100);
    let jsi = optimal
      .joint_spectrum(integrator)
      .jsi_normalized_range(range);
    // dbg!(&jsi);
    // check the max value isn't > 1
    let steps = range.as_steps();
    let max = jsi
      .iter()
      .enumerate()
      .fold((0.0_f64, 0. * RAD / S, 0. * RAD / S), |a, (i, &b)| {
        let (x, y) = steps.value(i);
        if b > a.0 {
          (b, x, y)
        } else {
          a
        }
      });
    dbg!(
      max.0,
      frequency_to_vacuum_wavelength(max.1),
      frequency_to_vacuum_wavelength(max.2)
    );
    assert!(max.0 <= 1.0);
  }
}
