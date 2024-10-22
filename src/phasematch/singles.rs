use super::*;
use crate::utils::frequency_to_wavenumber;
use crate::Complex;
use dim::ucum::{M, RAD};

/// Evaluate the fiber coupled singles phasematching function for a given set of frequencies
///
/// This places a "bucket collector" at the idler
#[allow(non_snake_case)]
pub fn phasematch_singles_fiber_coupling(
  omega_s: Frequency,
  omega_i: Frequency,
  spdc: &SPDC,
  integrator: Integrator,
) -> PerMeter3<f64> {
  let M2 = M * M; // meters squared
                  // crystal length
  let L = spdc.crystal_setup.length;

  let theta_s = spdc.signal.theta_internal();
  let phi_s = spdc.signal.phi();
  let theta_s_e = spdc.signal.theta_external(&spdc.crystal_setup);

  let Ws_SQ = spdc.signal.waist().x_by_y();

  let Wx_SQ = sq(spdc.pump.waist().x);
  let Wy_SQ = sq(spdc.pump.waist().y);

  // Counter-propagation requires a sign change in the wavenumber but
  // not for the free propagation constants ks_f and ki_f.
  let sign_ks = spdc.signal.direction().z.signum();
  let sign_ki = spdc.idler.direction().z.signum();

  let omega_p = omega_s + omega_i; // spdc.pump.frequency();
  let n_p = spdc.pump.refractive_index(omega_p, &spdc.crystal_setup);
  let k_p = frequency_to_wavenumber(omega_p, n_p);
  let n_s = spdc.signal.refractive_index(omega_s, &spdc.crystal_setup);
  let n_i = spdc.idler.refractive_index(omega_i, &spdc.crystal_setup);
  let k_s = sign_ks * frequency_to_wavenumber(omega_s, n_s);
  let k_i = sign_ki * frequency_to_wavenumber(omega_i, n_i);
  // let k_s = (spdc.signal.wavevector(omega_s, &spdc.crystal_setup) * M / RAD).z * RAD / M;
  // let k_i = (spdc.idler.wavevector(omega_i, &spdc.crystal_setup) * M / RAD).z * RAD / M;

  let PHI_s = cos(theta_s_e).powi(-2);

  let z0 = 0. * M; //put pump in middle of the crystal
  let z0s = spdc.signal_waist_position;

  // Height of the collected spots from the z axis.
  let hs = L * 0.5 * tan(theta_s) * cos(phi_s);
  // let hi = L * 0.5 * tan(theta_i) * cos(phi_i);

  let RHOpx = tan(spdc.pump.walkoff_angle(&spdc.crystal_setup));

  // Now calculate the the coeficients that get repeatedly used. This is from
  // Karina's code. Assume a symmetric pump waist (Wx = Wy)
  use dim::Abs;
  let ks_f = k_s.abs() / n_s; // exact
  let SIN_THETA_s_e = sin(theta_s_e); // 1e-9
  let COS_PHI_s = cos(phi_s);
  let GAM2s = -0.25 * Ws_SQ; // exact
  let GAM1s = GAM2s * PHI_s; // 1e-10
  let GAM3s = -2. * ks_f * GAM1s * SIN_THETA_s_e * COS_PHI_s; // 1e-10
  let GAM4s = -0.5 * ks_f * SIN_THETA_s_e * COS_PHI_s * GAM3s; // 1e-5
  let zhs = z0s + hs * SIN_THETA_s_e * COS_PHI_s; // 1e-13
  let DEL2s = (0.5 / ks_f) * zhs; // 1e-9
  let DEL1s = DEL2s * PHI_s; // 1e-9
  let DEL3s = -hs - zhs * PHI_s * SIN_THETA_s_e * COS_PHI_s; // 1e-11
  let KpKs = *(k_p * k_s * M2 / RAD / RAD); // exact

  let dksi = k_s + k_i + spdc.pp.k_eff();
  let C7 = k_p - dksi; // 1e-7
  let C3 = L * C7; // 1e-10
  let C4 = L * (1. / k_i - 1. / k_p); // 1e-13
  let C5 = k_s / k_p; // exact
  let C9 = Complex::new(*(k_p * Wx_SQ / M / RAD), 0.); // exact
  let C10 = Complex::new(*(k_p * Wy_SQ / M / RAD), 0.); // exact
  let LRho = L * RHOpx; // DIFFERENT SIGN and 1e-5
  let LRho_sq = LRho * LRho;

  let alpha1 = 4. * KpKs * Complex::new(*(GAM1s / M2), -*(DEL1s / M2 * RAD));
  let alpha2 = 4. * KpKs * Complex::new(*(GAM2s / M2), -*(DEL2s / M2 * RAD));
  let alpha3 = Complex::new(*(GAM3s / RAD / M), -*(DEL3s / M));

  let k_p_L = k_p * L;
  let KpKs4inv = 1. / (4. * KpKs);
  let imag = Complex::i();

  let fn_z = |z1: f64, z2: f64| {
    let B0 = z1 - z2;

    // krister broke this out of the integral so that repeat calculations didn't happen
    // over z1. Might not be necessary though.
    let A1 = 2. * z0 - L * z1;
    let B1 = 1. - z1;
    let B3 = 1. + z1;

    let A2 = 2. * z0 - L * z2;
    let B2 = 1. - z2;
    let B4 = 1. + z2;

    let B6a = *(C4 * B0 * RAD / M2);
    let gamma1 = *(-k_p_L * B1 / RAD + k_s * A1 / RAD) * imag; // exact
    let gamma2 = *(-k_p_L * B2 / RAD + k_s * A2 / RAD) * imag; // exact
    let Ha = alpha1 + gamma1;
    let Hb = alpha2 + gamma1;
    let Hc = alpha1.conj() - gamma2;
    let Hd = alpha2.conj() - gamma2;

    let ks = *(k_s * M / RAD);

    let AA1 = (Ha - C9 * ks) * KpKs4inv;
    let AA2 = (Hc - C9 * ks) * KpKs4inv;
    let BB1 = (Hb - C10 * ks) * KpKs4inv;
    let BB2 = (Hd - C10 * ks) * KpKs4inv;

    // TODO: verify with krister that this is correct in the original version
    let X11 = C9 * ks - Ha;
    let X12 = (Hc - C9 * ks) * imag;
    let Y21 = C10 * ks - Hb;
    let Y22 = (Hd - C10 * ks) * imag;

    // Now to calculate the term EE
    // EE = 1/4*(-  2*Wx^2 + I B6a + C5/X11*(C9 - I A1)^2 - I C5/X12*(C9 + I A2)^2  )
    let EE = 0.25
      * (-Complex::new(2. * (*(Wx_SQ / M2)), 0.)
        + imag * B6a
        + (*C5) / X11 * sq(C9 - imag * (*(A1 / M)))
        - imag * (*C5) / X12 * sq(C9 + imag * (*(A2 / M))));

    // Now to calculate the term FF
    // FF = 1/4*(-2*Wy^2 + I B6a - C5/Y21 *(I C10 + A1)^2 + I C5/Y22 *(-I C10 + A2)^2)
    let FF = 0.25
      * (-Complex::new(2. * (*(Wy_SQ / M2)), 0.) + imag * B6a
        - (*C5) / Y21 * sq(imag * C10 + (*(A1 / M)))
        + imag * (*C5) / Y22 * sq(-imag * C10 + (*(A2 / M))));

    // Now to calculate the term GG
    // GG = ks*( \[Alpha]3c/X12 *(I C9 - A2)  +  \[Alpha]3/X11 *(-C9 + I A1));
    let GG = ks
      * (alpha3.conj() / X12 * (imag * C9 - (*(A2 / M)))
        + alpha3 / X11 * (-C9 + imag * (*(A1 / M))));

    // Now to calculate the term HH
    // HH = L * \[Rho]/2 *(I B0 + ks*(B3/Y21 *(-I C10 - A1)  +  B4/Y22 *(C10 + I A2)));
    let HH = 0.5
      * (*(LRho / M))
      * (imag * B0
        + ks * (B3 / Y21 * (-imag * C10 - (*(A1 / M))) + B4 / Y22 * (C10 + imag * (*(A2 / M)))));

    // Now to calculate the term II
    // IIrho = 1/4* ks*kp*L^2*\[Rho]^2 ( -B3^2/Y21 +I B4^2/Y22)
    // IIgam = kp*ks*(\[Alpha]3^2/X11 - I \[Alpha]3c^2/X12)
    // IIdelk = 2 \[CapitalGamma]4s + 0.5 I (C3*B0)
    // II = IIrho + IIgam + IIdelk
    let IIrho = 0.25 * KpKs * (*(LRho_sq / M2)) * (-B3.powi(2) / Y21 + imag * B4.powi(2) / Y22);
    let IIgam = KpKs * (sq(alpha3) / X11 - imag * sq(alpha3.conj()) / X12);
    let IIdelk = 2. * *(GAM4s / RAD / RAD) + 0.5 * imag * *(C3 / RAD) * B0;
    let II = IIrho + IIgam + IIdelk;

    // Now calculate terms in the numerator
    // Exp(-(GG^2/(4 EE)) - HH^2/(4 FF) + II)
    let numerator = (-sq(GG) / (4. * EE) - sq(HH) / (4. * FF) + II).exp();

    // Now calculate terms in the Denominator
    // 8 * Sqrt[AA1 BB1 AA2 BB2 EE FF]
    let denominator = 8. * (AA1 * BB1 * AA2 * BB2 * EE * FF).sqrt();

    // Take into account apodized crystals
    // Apodization 1/e^2
    let pmzcoeff = spdc.pp.integration_constant(z1, L) * spdc.pp.integration_constant(z2, L);

    // Now calculate the full term in the integral.
    pmzcoeff * numerator / denominator
  };

  // let integrator = SimpsonIntegration2D::new(|z1, z2, _| fn_z(z1, z2));

  // // h(\omega_s, \omega_i) = \frac{1}{4} \int_{-1}^{1} d\xi_1 \int_{-1}^{1} d\xi_2 \psi(\xi_1, \xi_2).
  // let result = 0.25
  //   * integrator
  //     .integrate(
  //       (-1., 1.),
  //       (-1., 1.),
  //       steps.unwrap_or_else(|| integration_steps_best_guess(L)),
  //     )
  //     .norm();

  // let integrator = Integrator::AdaptiveSimpson {
  //   tolerance: 1e-5,
  //   max_depth: 1000,
  // };
  //
  // TODO: Where does this factor of 0.25 come from?
  let result = 0.25 * integrator.integrate2d(fn_z, -1., 1., -1., 1.).norm();

  PerMeter3::new(result)
}

#[cfg(test)]
mod tests {
  use super::*;

  fn percent_diff(actual: f64, expected: f64) -> f64 {
    100. * ((expected - actual) / expected).abs()
  }

  #[test]
  fn phasematch_singles_test() {
    let mut spdc = SPDC::default();
    spdc.crystal_setup.theta = 0.5515891191131287 * RAD;
    // spdc.signal.set_from_external_theta(0.0523598775598298 * RAD, &spdc.crystal_setup);
    spdc.signal.set_angles(0. * RAD, 0.03253866877817829 * RAD);
    // spdc.assign_optimum_idler();
    // spdc.assign_optimum_theta();
    spdc.idler.set_angles(PI * RAD, 0.03178987094602039 * RAD);

    spdc.signal_waist_position = -0.0007348996031796276 * M;
    spdc.idler_waist_position = -0.0007348996031796276 * M;

    let amp = *(phasematch_singles_fiber_coupling(
      spdc.signal.frequency(),
      spdc.idler.frequency(),
      &spdc,
      Integrator::default(),
    ) / PerMeter3::new(1.));

    let actual = amp;
    let expected = Complex::new(9.518_572_188_658_382e23, 95667755.72451791).norm();

    let accept_diff = 1e-4;

    let normdiff = percent_diff(actual, expected);
    assert!(
      normdiff < accept_diff,
      "norm percent difference: {}",
      normdiff
    );
  }

  #[test]
  fn phasematch_singles_pp_test() {
    let mut spdc = SPDC::default();
    spdc.pp = PeriodicPoling::On {
      sign: Sign::NEGATIVE,
      period: 0.000018041674656364844 * M,
      apodization: Apodization::Off,
    };
    // spdc.signal.set_from_external_theta(3. * DEG, &spdc.crystal_setup);
    spdc.signal.set_angles(0. * RAD, 0.0341877166429185 * RAD);
    // spdc.assign_optimum_idler();
    // spdc.assign_optimum_theta();

    // FIXME This isn't matching.
    spdc.idler.set_angles(PI * RAD, 0.031789820056487665 * RAD);
    spdc.crystal_setup.theta = 1.5707963267948966 * RAD;
    spdc.signal_waist_position = -0.0006311635856188344 * M;
    spdc.idler_waist_position = -0.0006311635856188344 * M;

    let amp = *(phasematch_singles_fiber_coupling(
      spdc.signal.frequency(),
      spdc.idler.frequency(),
      &spdc,
      Integrator::default(),
    ) / PerMeter3::new(1.));

    let actual = amp;
    let expected = Complex::new(1.667_581_141_397_712_8e24, -126659122.3067034).norm();

    let accept_diff = 1e-4;

    let normdiff = percent_diff(actual, expected);
    assert!(
      normdiff < accept_diff,
      "norm percent difference: {}",
      normdiff
    );
  }
}
