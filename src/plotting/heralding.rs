use crate::{jsa::*, utils::transpose_vec};
use super::*;
use spdc_setup::*;
use utils::Steps2D;
use std::f64::consts::{SQRT_2};
use math::{fwhm_to_sigma, sq};
use na::Vector2;
use dim::{
  ucum::{self, UCUM, M, S, C_, Meter, Hertz, RAD, EPS_0},
};

derived!(ucum, UCUM: RateConstUnits = Meter * Meter * Meter * Meter * Second * Second * Second );
/// Calculates \eta = \frac{2 \sqrt{2}}{\sqrt{\pi}} \frac{\mathbb{K}^2 d_{eff}^2 A_m^2}{\epsilon_0 c^3} (W_{sf}^2\sec\theta_{sf} ) \frac{W_{0x}W_{0y} L^2 P_{av}}{\sigma}.
/// with units m^6 s^3
#[allow(non_snake_case)]
fn calc_rate_constant(spdc_setup : &SPDCSetup) -> RateConstUnits<f64> {
  let PI2c = PI2 * C_;
  // K degeneracy factor (nonlinear polarization).
  // K = 1 for non-degenerate emission frequencies.
  // K = 1/2 for degenerate emission frequencies.
  let degeneracy_factor = 1.0_f64;
  // A_m = e^{i\frac{\pi}{2}m} \mathrm{sinc}\left ( {\frac{\pi}{2}m} \right ): periodic poling coefficient. Am = 1 for m = 0 (no periodic poling);
  // and A_m = 2i/\pi for m = 1 (sinusoidal variation)
  let periodic_poling_coeff = 1.0_f64;
  let constants = (2. * SQRT_2 / PI.sqrt())
    * degeneracy_factor.powi(2)
    * periodic_poling_coeff.powi(2) / (EPS_0 * C_ * C_ * C_);

  let p_bw = spdc_setup.pump_bandwidth;
  let lamda_p = spdc_setup.pump.get_wavelength();
  // TODO: ask krister why there's a factor of 2 here again. and why we're doing this for sigma
  let sigma = 2. * fwhm_to_sigma(PI2c * (
    1. / (lamda_p - 0.5 * p_bw)
    - 1. / (lamda_p + 0.5 * p_bw)
  ));

  let Wp = *(spdc_setup.pump.waist / M);
  let Wp_SQ = (Wp.x * Wp.y) * M * M;

  let L = spdc_setup.crystal_setup.length;
  let Pav = spdc_setup.pump_average_power;
  let deff = spdc_setup.crystal_setup.crystal.get_effective_nonlinear_coefficient();
  constants * sq(deff) * Wp_SQ * sq(L) * Pav / sigma
}

derived!(ucum, UCUM: CoincRateConstUnits = Meter * Meter * Meter * Meter * Meter * Meter * Meter * Meter * Second * Second * Second );
/// Calculates \eta = \frac{2 \sqrt{2}}{\sqrt{\pi}} \frac{\mathbb{K}^2 d_{eff}^2 A_m^2}{\epsilon_0 c^3} (W_{sf}^2\sec\theta_{sf} ) (W_{if}^2\sec\theta_{if}) \frac{W_{0x}W_{0y} L^2 P_{av}}{\sigma}.
/// with units (m^8 s^3)
#[allow(non_snake_case)]
fn calc_coincidence_rate_constant(spdc_setup : &SPDCSetup) -> CoincRateConstUnits<f64> {
  let Ws = *(spdc_setup.signal.waist / M);
  let Ws_SQ = (Ws.x * Ws.y) * M * M;
  let Wi = *(spdc_setup.idler.waist / M);
  let Wi_SQ = (Wi.x * Wi.y) * M * M;
  let theta_i_e = *(spdc_setup.idler.get_external_theta(&spdc_setup.crystal_setup) / RAD);
  let sec_i = 1. / f64::cos(theta_i_e);
  let theta_s_e = *(spdc_setup.signal.get_external_theta(&spdc_setup.crystal_setup) / RAD);
  let sec_s = 1. / f64::cos(theta_s_e);

  calc_rate_constant(&spdc_setup) * (Ws_SQ * sec_s) * (Wi_SQ * sec_i)
}

derived!(ucum, UCUM: JacobianDetUnits = Unitless / Second / Second );
// Used for coordinate transformation in double integrals
// l(\omega_s, \omega_i) = \frac{\omega_s \omega_i}{n_s^2(\omega_s) n_i^2(\omega_i)}
// but we input wavelengths
fn calc_jacobian_det_lambda_to_omega(l_s : Wavelength, l_i : Wavelength, spdc_setup : &SPDCSetup) -> JacobianDetUnits<f64> {
  let pi2c = PI2 * C_;
  let n_s = spdc_setup.crystal_setup.get_index_along(l_s, spdc_setup.signal.get_direction(), &spdc_setup.signal.get_type());
  let n_i = spdc_setup.crystal_setup.get_index_along(l_i, spdc_setup.idler.get_direction(), &spdc_setup.idler.get_type());
  let denominator = l_s * l_i * sq(n_s * n_i);
  pi2c * pi2c / denominator
}

/// Calculate the coincidence rates for each wavelength^2.
/// Integrating over all wavelengths (all array items) will give total coincidence rate.
#[allow(non_snake_case)]
pub fn calc_coincidences_rate_distribution(spdc_setup : &SPDCSetup, wavelength_range : &Steps2D<Wavelength>) -> Vec<Hertz<f64>> {

  let PI2c = PI2 * C_;

  let Ws = *(spdc_setup.signal.waist / M);
  let Wi = *(spdc_setup.idler.waist / M);

  if Ws.x == 0. || Ws.y == 0. || Wi.x == 0. || Wi.y == 0. {
    return vec![0. / S; wavelength_range.len()];
  }

  // apply angle adjustments for offset fiber theta
  let spdc_setup = spdc_setup.with_fiber_theta_offsets_applied();

  // let n_p = spdc_setup.pump.get_index(&spdc_setup.crystal_setup);

  // NOTE: original version incorrectly computed the step size.
  // originally it was (x_f - x_i) / n_{points}, which should be (x_f - x_i / (n-1))
  let (dlamda_s, dlamda_i) = wavelength_range.division_widths();

  let eta = calc_coincidence_rate_constant(&spdc_setup);

  let jsa_units = JSAUnits::new(1.);

  // NOTE: moving this inside the integral and using running lamda values
  // results in efficiency change of ~0.01%

  // TODO: check with krister
  // NOTE: this was the old way
  let l_s = spdc_setup.signal.get_wavelength();
  let l_i = spdc_setup.idler.get_wavelength();
  let d_omega_s = PI2c * dlamda_s / sq(l_s);
  let d_omega_i = PI2c * dlamda_i / sq(l_i);

  wavelength_range
    .into_iter()
    .map(|(l_s, l_i)| {
      // let d_omega_s = PI2c * dlamda_s / sq(l_s);
      // let d_omega_i = PI2c * dlamda_i / sq(l_i);
      // dbg!(d_omega_s, d_omega_i);
      let lomega = calc_jacobian_det_lambda_to_omega(l_s, l_i, &spdc_setup);
      let factor = eta * lomega * d_omega_s * d_omega_i;
      let amplitude = calc_jsa(&spdc_setup, l_s, l_i) / jsa_units;

      amplitude.norm_sqr() * sq(jsa_units) * factor
    })
    .collect()
}

/// Calculate the singles rate per unit wavelength^2.
/// Integrating over all wavelengths (all array items) will give total singles rate.
/// technically this distribution has units of 1/(s m^2).. but we return 1/s
#[allow(non_snake_case)]
pub fn calc_singles_rate_distribution_signal(spdc_setup : &SPDCSetup, wavelength_range : &Steps2D<Wavelength>) -> Vec<Hertz<f64>> {

  let PI2c = PI2 * C_;

  let Ws = *(spdc_setup.signal.waist / M);
  let Wi = *(spdc_setup.idler.waist / M);

  // TODO: ask krister, why are we squaring waist like this?
  let Ws_SQ = (Ws.x * Ws.y) * M * M;
  let Wi_SQ = (Wi.x * Wi.y) * M * M;

  // early out...
  if *(Ws_SQ / M / M) == 0. || *(Wi_SQ / M / M) == 0. {
    return vec![0. / S; wavelength_range.len()];
  }

  // apply fiber offset theta for the signal only
  let mut spdc_setup = spdc_setup.with_signal_fiber_theta_offsets_applied();
  // place the idler at the correct angle for this adjusted signal position
  // because we want to have a "bucket collector" at the idler channel
  // to collect all modes
  spdc_setup.assign_optimum_idler();

  let theta_s_e = *(spdc_setup.signal.get_external_theta(&spdc_setup.crystal_setup) / RAD);

  // let n_p = spdc_setup.pump.get_index(&spdc_setup.crystal_setup);
  let PHI_s = 1. / f64::cos(theta_s_e);

  let scale_s = Ws_SQ * PHI_s;
  let (dlamda_s, dlamda_i) = wavelength_range.division_widths();

  let eta_s = calc_rate_constant(&spdc_setup);

  let jsa_units = JSAUnits::new(1.);

  wavelength_range
    .into_iter()
    .map(|(l_s, l_i)| {
      let d_omega_s = PI2c * dlamda_s / sq(l_s);
      let d_omega_i = PI2c * dlamda_i / sq(l_i);

      let lomega = calc_jacobian_det_lambda_to_omega(l_s, l_i, &spdc_setup);
      let f = eta_s * d_omega_s * d_omega_i * lomega;
      let factor_s = scale_s * f;

      let amplitude_s = *(calc_jsa_singles(&spdc_setup, l_s, l_i) / jsa_units);
      // technically this distribution has units of 1/(s m^2).. but we return 1/s
      amplitude_s.norm() * factor_s * jsa_units / M / M
    })
    .collect()
}

/// Calculate the singles rate for both signal and idler per unit wavelength^2.
/// Integrating over all wavelengths (all array items) will give total singles rate.
/// technically this distribution has units of 1/(s m^2).. but we return 1/s
#[allow(non_snake_case)]
pub fn calc_singles_rate_distributions(spdc_setup : &SPDCSetup, wavelength_range : &Steps2D<Wavelength>) -> (Vec<Hertz<f64>>, Vec<Hertz<f64>>) {

  let PI2c = PI2 * C_;

  let Ws = *(spdc_setup.signal.waist / M);
  let Wi = *(spdc_setup.idler.waist / M);

  // TODO: ask krister, why are we squaring waist like this?
  let Ws_SQ = (Ws.x * Ws.y) * M * M;
  let Wi_SQ = (Wi.x * Wi.y) * M * M;

  // early out...
  if *(Ws_SQ / M / M) == 0. || *(Wi_SQ / M / M) == 0. {
    return (vec![0. / S; wavelength_range.len()], vec![0. / S; wavelength_range.len()]);
  }

  // spdc_setup copy with signal and idler swapped
  let mut spd_swap = spdc_setup.with_swapped_signal_idler();

  // apply fiber offset theta for the signal only
  let mut spdc_setup = spdc_setup.with_signal_fiber_theta_offsets_applied();
  // place the idler at the correct angle for this adjusted signal position
  // because we want to have a "bucket collector" at the idler channel
  // to collect all modes
  spdc_setup.assign_optimum_idler();

  // remember... spd_swap has signal and idler swapped...
  spd_swap = spd_swap.with_signal_fiber_theta_offsets_applied();
  // place the "idler" (signal) at the correct angle for this adjusted signal position
  // because we want to have a "bucket collector" at the "idler" channel
  // to collect all modes
  spd_swap.assign_optimum_idler();

  // constants that work for both channels...
  // let n_p = spdc_setup.pump.get_index(&spdc_setup.crystal_setup);
  let (dlamda_s, dlamda_i) = wavelength_range.division_widths();
  let eta_s = calc_rate_constant(&spdc_setup);

  let jsa_units = JSAUnits::new(1.);

  // constants for signal or idler channel
  let theta_s_e = *(spdc_setup.signal.get_external_theta(&spdc_setup.crystal_setup) / RAD);
  let PHI_s = 1. / f64::cos(theta_s_e);

  let theta_i_e = *(spdc_setup.idler.get_external_theta(&spdc_setup.crystal_setup) / RAD);
  let PHI_i = 1. / f64::cos(theta_i_e);

  let scale_s = Ws_SQ * PHI_s;
  let scale_i = Wi_SQ * PHI_i;

  let (signal, idler) = wavelength_range
    .into_iter()
    .map(|(l_s, l_i)| {
      let d_omega_s = PI2c * dlamda_s / sq(l_s);
      let d_omega_i = PI2c * dlamda_i / sq(l_i);

      let f = eta_s * d_omega_s * d_omega_i;

      // it's symmetric so no need to swap
      let lomega = calc_jacobian_det_lambda_to_omega(l_s, l_i, &spdc_setup);
      let factor_s = scale_s * lomega * f;
      let factor_i = scale_i * lomega * f;

      let amplitude_s = *(calc_jsa_singles(&spdc_setup, l_s, l_i) / jsa_units);
      // swap signal/idler wavelengths for spd_swap...
      let amplitude_i = *(calc_jsa_singles(&spd_swap, l_i, l_s) / jsa_units);
      // technically this distribution has units of 1/(s m^2).. but we return 1/s
      (
        amplitude_s.norm() * factor_s * jsa_units / M / M,
        amplitude_i.norm() * factor_i * jsa_units / M / M,
      )
    })
    .unzip();

  (signal, idler)
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub struct HeraldingResults {
  pub signal_singles_rate : f64,
  pub idler_singles_rate : f64,
  pub coincidences_rate : f64,
  pub signal_efficiency : f64,
  pub idler_efficiency : f64,
  pub symmetric_efficiency : f64,
}

impl HeraldingResults {
  pub fn from_distributions(
    coincidences_rate_distribution : Vec<Hertz<f64>>,
    signal_singles_rate_distribution : Vec<Hertz<f64>>,
    idler_singles_rate_distribution: Vec<Hertz<f64>>
  ) -> Self {

    let coincidences_rate = coincidences_rate_distribution.iter()
      .map(|&r| *(r * S)).sum();

    let signal_singles_rate = signal_singles_rate_distribution.iter()
      .map(|&r| *(r * S)).sum();

    let idler_singles_rate = idler_singles_rate_distribution.iter()
      .map(|&r| *(r * S)).sum();

    let signal_efficiency = if idler_singles_rate == 0. { 0. } else { coincidences_rate / idler_singles_rate };
    let idler_efficiency = if signal_singles_rate == 0. { 0. } else { coincidences_rate / signal_singles_rate };
    let symmetric_efficiency = if signal_singles_rate == 0. || idler_singles_rate == 0. {
      0.
    } else {
      let denom : f64 = signal_singles_rate * idler_singles_rate;
      coincidences_rate / denom.sqrt()
    };

    HeraldingResults {
      signal_singles_rate,
      idler_singles_rate,
      coincidences_rate,
      signal_efficiency,
      idler_efficiency,
      symmetric_efficiency,
    }
  }
}

/// Get the heralding results for a given spdc_setup setup and signal/idler wavelength range
pub fn calc_heralding_results(spdc_setup : &SPDCSetup, wavelength_range : &Steps2D<Wavelength>) -> HeraldingResults {
  let coincidences_rate_distribution = calc_coincidences_rate_distribution(&spdc_setup, &wavelength_range);
  let (signal, idler) = calc_singles_rate_distributions(&spdc_setup, &wavelength_range);

  HeraldingResults::from_distributions(
    coincidences_rate_distribution,
    signal,
    idler
  )
}

/// Calculate the count rates, and efficiencies for signal, idler singles and coincidences
/// as well as the efficiencies over a range of signal/idler waist sizes.
pub fn plot_heralding_results_by_signal_idler_waist(
  spdc_setup : &SPDCSetup,
  si_waists : &Steps2D<Meter<f64>>,
  wavelength_range : &Steps2D<Wavelength>
) -> Vec<HeraldingResults> {
  si_waists
    .into_iter()
    .map(|(ws, wi)| {
      let mut spdc_setup = spdc_setup.clone();
      spdc_setup.signal.waist = Meter::new(Vector2::new(*(ws / M), *(ws / M)));
      spdc_setup.idler.waist = Meter::new(Vector2::new(*(wi / M), *(wi / M)));

      calc_heralding_results(&spdc_setup, &wavelength_range)
    })
    .collect()
}

/// Calculate the count rates, and efficiencies for signal, idler singles and coincidences
/// as well as the efficiencies over a range of pump vs signal/idler waist sizes.
pub fn plot_heralding_results_by_pump_signal_idler_waist(
  spdc_setup : &SPDCSetup,
  ps_waists : &Steps2D<Meter<f64>>,
  wavelength_range : &Steps2D<Wavelength>
) -> Vec<HeraldingResults> {
  ps_waists
    .into_iter()
    .map(|(wp, ws)| {
      let mut spdc_setup = spdc_setup.clone();
      spdc_setup.pump.waist = Meter::new(Vector2::new(*(wp / M), *(wp / M)));
      spdc_setup.signal.waist = Meter::new(Vector2::new(*(ws / M), *(ws / M)));
      spdc_setup.idler.waist = spdc_setup.signal.waist.clone();

      calc_heralding_results(&spdc_setup, &wavelength_range)
    })
    .collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  // extern crate float_cmp;
  // use float_cmp::*;
  use dim::f64prefixes::{NANO, MICRO};
  use ucum::DEG;

  fn percent_diff(actual : f64, expected : f64) -> f64 {
    100. * ((expected - actual) / expected).abs()
  }

  #[test]
  fn zero_rates_test() {
    let mut spdc_setup = SPDCSetup::default();
    spdc_setup.fiber_coupling = true;
    spdc_setup.crystal_setup.crystal = Crystal::KTP;
    spdc_setup.signal.waist = Meter::new(Vector2::new(0., 0.));
    spdc_setup.idler.waist = Meter::new(Vector2::new(0., 0.));
    spdc_setup.assign_optimum_periodic_poling();

    let wavelength_range = Steps2D(
      (1490.86 * NANO * M, 1609.14 * NANO * M, 10),
      (1495.05 * NANO * M, 1614.03 * NANO * M, 10)
    );

    let results = calc_heralding_results(&spdc_setup, &wavelength_range);
    assert_eq!(results.coincidences_rate, 0.);
  }

  #[test]
  fn coincidence_rates_test() {
    let mut spdc_setup = SPDCSetup::default();
    spdc_setup.fiber_coupling = true;
    spdc_setup.crystal_setup.crystal = Crystal::KTP;
    spdc_setup.crystal_setup.theta = 90. * DEG;
    spdc_setup.crystal_setup.length = 6000. * MICRO * M;
    spdc_setup.pump.waist = WaistSize::new(Vector2::new(200. * MICRO, 200. * MICRO));
    spdc_setup.assign_optimum_periodic_poling();

    let wavelength_range = Steps2D(
      (1546.02 * NANO * M, 1553.99 * NANO * M, 30),
      (1542.05 * NANO * M, 1550.00 * NANO * M, 30)
    );

    let rates = calc_coincidences_rate_distribution(&spdc_setup, &wavelength_range);
    let actual = rates.iter().map(|&r| *(r * S)).sum();
    let expected = 3800.0;

    let accept_diff = 1e-9; // percent difference accpetable
    let pdiff = percent_diff(actual, expected);

    assert!(
      pdiff < accept_diff,
      "norm percent difference: {}. (actual: {}, expected: {})",
      pdiff,
      actual,
      expected
    );
  }

  #[test]
  fn singles_rates_test() {
    let mut spdc_setup = SPDCSetup::default();
    spdc_setup.fiber_coupling = true;
    spdc_setup.crystal_setup.crystal = Crystal::KTP;
    spdc_setup.crystal_setup.theta = 90. * DEG;
    spdc_setup.crystal_setup.length = 6000. * MICRO * M;
    spdc_setup.pump.waist = WaistSize::new(Vector2::new(200. * MICRO, 200. * MICRO));
    spdc_setup.assign_optimum_periodic_poling();

    let wavelength_range = Steps2D(
      (1546.02 * NANO * M, 1553.99 * NANO * M, 30),
      (1542.05 * NANO * M, 1550.00 * NANO * M, 30)
    );

    let (signal, idler) = calc_singles_rate_distributions(&spdc_setup, &wavelength_range);
    let actual = signal.iter().zip(idler.iter())
      .map(|(&r_s, &r_i)|
        ( *(r_s * S), *(r_i * S) )
      )
      // sum tuples
      .fold((0., 0.), |col, (r_s, r_i)|
        (col.0 + r_s, col.1 + r_i)
      );

    let expected = (3901.133069046596, 3901.3937911182834);

    let accept_diff = 1e-9; // percent difference accpetable
    let pdiff_s = percent_diff(actual.0, expected.0);

    assert!(
      pdiff_s < accept_diff,
      "norm percent difference: {}. (actual: {}, expected: {})",
      pdiff_s,
      actual.0,
      expected.0
    );

    let pdiff_i = percent_diff(actual.1, expected.1);

    assert!(
      pdiff_i < accept_diff,
      "norm percent difference: {}. (actual: {}, expected: {})",
      pdiff_i,
      actual.1,
      expected.1
    );
  }

  #[test]
  fn efficiency_apodization_test() {
    let mut spdc_setup = SPDCSetup::default();
    spdc_setup.fiber_coupling = true;
    spdc_setup.crystal_setup.crystal = Crystal::KTP;
    spdc_setup.crystal_setup.length = 1000. * MICRO * M;
    spdc_setup.assign_optimum_periodic_poling();
    spdc_setup.pp = spdc_setup.pp.map(|pp| {
      PeriodicPoling {
        apodization: Some(Apodization {
          fwhm: 1000. * MICRO * M
        }),
        ..pp
      }
    });
    // spdc_setup.assign_optimum_periodic_poling();

    let wavelength_range = Steps2D(
      (1490.86 * NANO * M, 1609.14 * NANO * M, 30),
      (1495.05 * NANO * M, 1614.03 * NANO * M, 30)
    );

    let coinc_rate_distr = calc_coincidences_rate_distribution(&spdc_setup, &wavelength_range);
    let (signal, idler) = calc_singles_rate_distributions(&spdc_setup, &wavelength_range);

    let results = HeraldingResults::from_distributions(coinc_rate_distr, signal, idler);

    let accept_diff = 1e-1;

    // old bugged code would have given
    // coinc rate sum 3005.068611324783
    // singles s rate sum 3448.594543844433
    // singles i rate sum 3448.71420649685
    // idler efficiency 0.8713893654702548
    // signal efficiency 0.8713591302125567

    let expected_idler = 0.8890414777632809;
    let pdiff_i = percent_diff(results.idler_efficiency, expected_idler);
    assert!(
      pdiff_i < accept_diff,
      "norm percent difference: {}. (actual: {}, expected: {})",
      pdiff_i,
      results.idler_efficiency,
      expected_idler
    );

    let expected_signal = 0.8890107154534013;
    let pdiff_s = percent_diff(results.signal_efficiency, expected_signal);

    assert!(
      pdiff_s < accept_diff,
      "norm percent difference: {}. (actual: {}, expected: {})",
      pdiff_s,
      results.signal_efficiency,
      expected_signal
    );
  }

  #[test]
  fn efficiency_test() {
    let mut spdc_setup = SPDCSetup::default();
    spdc_setup.fiber_coupling = true;
    spdc_setup.crystal_setup.crystal = Crystal::KTP;
    spdc_setup.crystal_setup.theta = 90. * DEG;
    spdc_setup.crystal_setup.length = 6000. * MICRO * M;
    spdc_setup.pump.waist = WaistSize::new(Vector2::new(200. * MICRO, 200. * MICRO));
    spdc_setup.assign_optimum_periodic_poling();

    let wavelength_range = Steps2D(
      (1546.02 * NANO * M, 1553.99 * NANO * M, 30),
      (1542.05 * NANO * M, 1550.00 * NANO * M, 30)
    );

    let coinc_rate_distr = calc_coincidences_rate_distribution(&spdc_setup, &wavelength_range);
    let (signal, idler) = calc_singles_rate_distributions(&spdc_setup, &wavelength_range);

    let results = HeraldingResults::from_distributions(coinc_rate_distr, signal, idler);
    dbg!(results);
    let accept_diff = 1e-1;

    // old bugged code would give
    // coinc rate sum 8767.90113100421
    // singles s rate sum 10192.880932347976
    // singles i rate sum 10193.119015740884
    // idler efficiency 0.8601985237734435
    // signal efficiency 0.860178431887653

    let expected_signal = 0.9783834540624885;
    let pdiff_s = percent_diff(results.signal_efficiency, expected_signal);

    assert!(
      pdiff_s < accept_diff,
      "norm percent difference: {}. (actual: {}, expected: {})",
      pdiff_s,
      results.signal_efficiency,
      expected_signal
    );

    let expected_idler = 0.9784488417733235;
    let pdiff_i = percent_diff(results.idler_efficiency, expected_idler);
    assert!(
      pdiff_i < accept_diff,
      "norm percent difference: {}. (actual: {}, expected: {})",
      pdiff_i,
      results.idler_efficiency,
      expected_idler
    );
  }
}
