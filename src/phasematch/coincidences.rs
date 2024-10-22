use super::*;
use crate::utils::frequency_to_wavenumber;
use dim::ucum::{M, RAD};
use dim::{ucum::UCUM, Dimensioned};

/// Evaluate the phasematching function using a gaussian approximation
#[allow(non_snake_case)]
pub fn phasematch_gaussian(
  omega_s: Frequency,
  omega_i: Frequency,
  spdc: &SPDC,
) -> PerMeter4<Complex<f64>> {
  let L = spdc.crystal_setup.length;
  let delk = *(spdc.delta_k(omega_s, omega_i) / Wavenumber::new(1.));
  let delta_k_z = Wavenumber::new(1.) * delk.z;
  let arg = L * 0.5 * delta_k_z;
  let pmz = Complex::new(gaussian_pm(*(arg / RAD)), 0.);

  PerMeter4::new(pmz)
}

/// Evaluate the phasematching function using a sinc approximation
#[allow(non_snake_case)]
pub fn phasematch_sinc(
  omega_s: Frequency,
  omega_i: Frequency,
  spdc: &SPDC,
) -> PerMeter4<Complex<f64>> {
  let L = spdc.crystal_setup.length;
  let delk = *(spdc.delta_k(omega_s, omega_i) / Wavenumber::new(1.));
  let delta_k_z = Wavenumber::new(1.) * delk.z;
  let arg = L * 0.5 * delta_k_z;
  // no fiber coupling
  let pmz = Complex::new(sinc(arg), 0.);
  let waist = spdc.pump.waist();
  let pmt = (-0.5 * (sq(delk.x * waist.x / M) + sq(delk.y * waist.y / M))).exp();

  PerMeter4::new(pmz * pmt)
}

/// Get the phasematching function integrand for a given setup at given frequencies
#[allow(non_snake_case)]
pub fn get_pm_integrand(
  omega_s: Frequency,
  omega_i: Frequency,
  spdc: &SPDC,
) -> impl Fn(f64) -> Complex<f64> + '_ {
  let L = spdc.crystal_setup.length;

  let phi_s = spdc.signal.phi();
  let phi_i = spdc.idler.phi();
  let theta_s = spdc.signal.theta_internal();
  let theta_i = spdc.idler.theta_internal();
  let theta_s_e = spdc.signal.theta_external(&spdc.crystal_setup);
  let theta_i_e = spdc.idler.theta_external(&spdc.crystal_setup);

  let Ws_SQ = spdc.signal.waist().x_by_y();
  let Wi_SQ = spdc.idler.waist().x_by_y();

  let Wx_SQ = sq(spdc.pump.waist().x);
  let Wy_SQ = sq(spdc.pump.waist().y);

  let SEC_2_THETA_s = cos(theta_s_e).powi(-2);
  let SEC_2_THETA_i = cos(theta_i_e).powi(-2);

  let z0 = 0. * M; //put pump in middle of the crystal
  let z0s = spdc.signal_waist_position;
  let z0i = spdc.idler_waist_position;

  // Height of the collected spots from the z axis.
  let hs = L * 0.5 * tan(theta_s) * cos(phi_s);
  let hi = L * 0.5 * tan(theta_i) * cos(phi_i);

  let SIN_THETA_s_e = sin(theta_s_e);
  let SIN_THETA_i_e = sin(theta_i_e);
  let TAN_THETA_s_e = tan(theta_s_e);
  let TAN_THETA_i_e = tan(theta_i_e);
  let COS_PHI_s = cos(phi_s);
  let COS_PHI_i = cos(phi_i);

  let omega_p = omega_s + omega_i;

  // Counter-propagation requires a sign change in the wavenumber but
  // not for the free propagation constants ks_f and ki_f.
  let sign_ks = spdc.signal.direction().z.signum();
  let sign_ki = spdc.idler.direction().z.signum();

  let n_p = spdc.pump.refractive_index(omega_p, &spdc.crystal_setup);
  let k_p = frequency_to_wavenumber(omega_p, n_p);
  let n_s = spdc.signal.refractive_index(omega_s, &spdc.crystal_setup);
  let n_i = spdc.idler.refractive_index(omega_i, &spdc.crystal_setup);
  let k_s = sign_ks * frequency_to_wavenumber(omega_s, n_s);
  let k_i = sign_ki * frequency_to_wavenumber(omega_i, n_i);

  // Now calculate the the coeficients that get repeatedly used.
  // This is from Karina's code.

  use dim::Abs;
  let ks_f = k_s.abs() / n_s;
  let ki_f = k_i.abs() / n_i;

  let GAM2s = -0.25 * Ws_SQ;
  let GAM2i = -0.25 * Wi_SQ;
  let GAM1s = GAM2s * SEC_2_THETA_s;
  let GAM1i = GAM2i * SEC_2_THETA_i;

  let GAM3s = -2. * ks_f * GAM1s * SIN_THETA_s_e * COS_PHI_s;
  let GAM3i = -2. * ki_f * GAM1i * SIN_THETA_i_e * COS_PHI_i;
  let GAM4s = -0.5 * ks_f * SIN_THETA_s_e * COS_PHI_s * GAM3s;
  let GAM4i = -0.5 * ki_f * SIN_THETA_i_e * COS_PHI_i * GAM3i;
  let zhs = z0s + hs * SIN_THETA_s_e * COS_PHI_s;
  let zhi = z0i + hi * SIN_THETA_i_e * COS_PHI_i;
  let DEL2s = (0.5 / ks_f) * zhs * RAD;
  let DEL2i = (0.5 / ki_f) * zhi * RAD;
  let DEL1s = DEL2s * SEC_2_THETA_s;
  let DEL1i = DEL2i * SEC_2_THETA_i;
  let DEL3s = -hs - zhs * SEC_2_THETA_s * SIN_THETA_s_e * COS_PHI_s;
  let DEL3i = -hi - zhi * SEC_2_THETA_i * SIN_THETA_i_e * COS_PHI_i;

  let DEL4s = 0.5 * ks_f * zhs * TAN_THETA_s_e.powi(2) - ks_f * z0s;
  let DEL4i = 0.5 * ki_f * zhi * TAN_THETA_i_e.powi(2) - ki_f * z0i;

  // dbg!(
  //   zhs, zhi, GAM1s, GAM1i, GAM2s, GAM2i, GAM3s, GAM3i, GAM4s, GAM4i, DEL1s, DEL1i, DEL2s, DEL2i,
  //   DEL3s, DEL3i, DEL4s, DEL4i
  // );

  let M2 = M * M; // meters squared
                  // let As_r = -0.25 * Wx_SQ + GAM1s;
                  // let As_i = -DEL1s;
  let As = Complex::new(*((-0.25 * Wx_SQ + GAM1s) / M2), *(-DEL1s / M2));
  // let Ai_r = -0.25 * Wx_SQ + GAM1i;
  // let Ai_i = -DEL1i;
  let Ai = Complex::new(*((-0.25 * Wx_SQ + GAM1i) / M2), *(-DEL1i / M2));
  // let Bs_r = -0.25 * Wy_SQ + GAM2s;
  // let Bs_i = -DEL2s;
  let Bs = Complex::new(*((-0.25 * Wy_SQ + GAM2s) / M2), *(-DEL2s / M2));
  // let Bi_r = -0.25 * Wy_SQ + GAM2i;
  // let Bi_i = -DEL2i;
  let Bi = Complex::new(*((-0.25 * Wy_SQ + GAM2i) / M2), *(-DEL2i / M2));
  let Cs = -0.25 * (L / k_s - 2. * z0 / k_p);
  let Ci = -0.25 * (L / k_i - 2. * z0 / k_p);
  let Ds = 0.25 * L * (1. / k_s - 1. / k_p);
  let Di = 0.25 * L * (1. / k_i - 1. / k_p);
  // let mx_real = -0.50 * Wx_SQ;
  // let mx_imag = z0/k_p;
  let mx = Complex::new(*((-0.50 * Wx_SQ) / M2), *((z0 / k_p) * RAD / M2));
  // let my_real = -0.50 * Wy_SQ;
  // let my_imag = mx_imag;
  let my = Complex::new(*((-0.50 * Wy_SQ) / M2), *((z0 / k_p) * RAD / M2));
  let m = L / (2. * k_p);
  let n = 0.5 * L * tan(spdc.pump.walkoff_angle(&spdc.crystal_setup));

  let hh = Complex::new(
    *(GAM4s / RAD / RAD + GAM4i / RAD / RAD),
    -*(DEL4s / RAD + DEL4i / RAD),
  );

  // let A5R = GAM3s;
  // let A5I = -DEL3s;
  let A5 = Complex::new(*(GAM3s / RAD / M), -*(DEL3s / M));
  let A5sq = A5 * A5;
  // let A7R = GAM3i;
  // let A7I = -DEL3i;
  let A7 = Complex::new(*(GAM3i / RAD / M), -*(DEL3i / M));

  let dksi = k_s + k_i + spdc.pp.k_eff();
  let ee = 0.5 * L * (k_p + dksi);
  let ff = 0.5 * L * (k_p - dksi);

  // dbg!(Cs, Ci, Ds, Di);

  // dbg!(As, Ai, Bs, Bi, Cs, Ci, Ds, Di);

  

  move |z: f64| {
    let Ds_z = Ds * z;
    let Di_z = Di * z;
    let CsDs = Complex::new(0., *((Cs + Ds_z) * RAD / M2));
    let CiDi = Complex::new(0., *((Ci + Di_z) * RAD / M2));

    // let A1R = As_r;
    // let A1I = As_i + Cs + Ds_z;
    let A1 = As + CsDs;

    // let A2R = Bs_r;
    // let A2I = Bs_i + Cs + Ds_z;
    let A2 = Bs + CsDs;

    // let A3R = Ai_r;
    // let A3I = Ai_i + Ci + Di_z;
    let A3 = Ai + CiDi;

    // let A4R = Bi_r;
    // let A4I = Bi_i + Ci + Di_z;
    let A4 = Bi + CiDi;

    // let A6R = 0.;
    // let A6I = n * (1 + z);
    let A6 = Complex::new(0., *(n / M) * (1. + z));

    // let A8R = mx_real;
    // let A8I = mx_imag - m * z;
    let mz = Complex::new(0., *(m * z * RAD / M2));
    let A8 = mx - mz;

    // let A9R = my_real;
    // let A9I = my_imag - m * z;
    let A9 = my - mz;

    // let A10R = hh_r;
    // let A10I = hh_i + ee + ff * z;
    let A10 = hh + Complex::new(0., *((ee + ff * z) / RAD));

    // First calculate terms in the exponential of the integral
    // exp(1/4 (4 A10 - A5^2/A1 - A6^2/A2 - (-2 A1 A7 + A5 A8)^2/(A1 (4 A1 A3 - A8^2)) - (A6^2 (-2 A2 + A9)^2)/(A2 (4 A2 A4 - A9^2))))
    let A6sq = A6 * A6;
    let A8sq = A8 * A8;
    let A9sq = A9 * A9;
    let invA1 = A1.inv();
    let invA2 = A2.inv();
    let term4 = -2. * A1 * A7 + A5 * A8;
    let term5 = -2. * A2 + A9;
    let denom1 = 4. * A1 * A3 - A8sq;
    let denom2 = 4. * A2 * A4 - A9sq;

    // dbg!(A1, A2);

    let numerator = ((4. * A10
      - invA1 * (A5sq + (term4 * term4) / denom1)
      - invA2 * A6sq * (1. + (term5 * term5) / denom2))
      / 4.)
      .exp();

    // Now deal with the denominator in the integral:
    // Sqrt[A1 A2 (-4 A3 + A8^2/A1) (-4 A4 + A9^2/A2)]
    let denominator = (denom1 * denom2).sqrt();

    // dbg!(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10);

    // Take into account apodized crystals
    let pmzcoeff = spdc.pp.integration_constant(z, L);

    // dbg!(pmzcoeff, z, numerator, denominator);
    // Now calculate the full term in the integral.
    pmzcoeff * numerator / denominator
  }
}

/// Evaluate the phasematching function for fiber coupling
///
/// This is the secret sauce of spdcalc.
#[allow(non_snake_case)]
pub fn phasematch_fiber_coupling(
  omega_s: Frequency,
  omega_i: Frequency,
  spdc: &SPDC,
  integrator: Integrator,
) -> PerMeter4<Complex<f64>> {
  // return phasematch_fiber_coupling2(omega_s, omega_i, spdc, integrator);
  // return phasematch_fiber_coupling_v3(omega_s, omega_i, spdc, integrator);
  // return phasematch_sinc(omega_s, omega_i, spdc);

  let fn_z = get_pm_integrand(omega_s, omega_i, spdc);
  // TODO: Where does this factor of 0.5 come from?
  let result = 0.5 * integrator.integrate(&fn_z, -1., 1.);
  PerMeter4::new(result)
}

/// Get a dimensioned complex number
pub fn complex_dim<U>(re: UCUM<f64, U>, im: UCUM<f64, U>) -> UCUM<Complex<f64>, U> {
  UCUM::<Complex<f64>, U>::new(Complex::new(*re.value_unsafe(), *im.value_unsafe()))
}

/// Get a dimensioned complex number with zero real part
pub fn imag_dim<U>(im: UCUM<f64, U>) -> UCUM<Complex<f64>, U> {
  UCUM::<Complex<f64>, U>::new(Complex::new(0., *im.value_unsafe()))
}

/// Get a dimensioned complex number with zero imaginary part
pub fn real_dim<U>(re: UCUM<f64, U>) -> UCUM<Complex<f64>, U> {
  UCUM::<Complex<f64>, U>::new(Complex::new(*re.value_unsafe(), 0.))
}

/// Evaluate the phasematching function for fiber coupling
///
/// This is the secret sauce of spdcalc.
#[allow(non_snake_case)]
pub fn phasematch_fiber_coupling2(
  omega_s: Frequency,
  omega_i: Frequency,
  spdc: &SPDC,
  integrator: Integrator,
) -> PerMeter4<Complex<f64>> {
  // crystal length
  let L = spdc.crystal_setup.length;
  let half_L = L * 0.5;

  let phi_s = spdc.signal.phi();
  let phi_i = spdc.idler.phi();

  // NOTE: This implementation assumes you are in the Phi = 0 plane
  // But we need to account for the directionality of the signal and idler
  // in that plane.
  // This accounts for the directionality of signal and idler
  // by negating theta if the angle is greater than 90 degrees.
  let sign_s = sin(phi_s).signum();
  let sign_i = sin(phi_i).signum();

  let theta_s = sign_s * spdc.signal.theta_internal();
  let theta_i = sign_i * spdc.idler.theta_internal();
  let theta_s_e = sign_s * spdc.signal.theta_external(&spdc.crystal_setup);
  let theta_i_e = sign_i * spdc.idler.theta_external(&spdc.crystal_setup);

  let Ws_SQ = spdc.signal.waist().x_by_y();
  let Wi_SQ = spdc.idler.waist().x_by_y();

  let Wx_SQ = sq(spdc.pump.waist().x);
  let Wy_SQ = sq(spdc.pump.waist().y);

  let SEC_2_THETA_s = cos(theta_s_e).powi(-2);
  let SEC_2_THETA_i = cos(theta_i_e).powi(-2);

  let ζ = 0. * M; //put pump in middle of the crystal
  let ζs = spdc.signal_waist_position;
  let ζi = spdc.idler_waist_position;

  // Height of the collected spots from the z axis.
  let hs = half_L * tan(theta_s);
  let hi = half_L * tan(theta_i);

  let SIN_THETA_s_e = sin(theta_s_e);
  let SIN_THETA_i_e = sin(theta_i_e);
  let TAN_THETA_s_e = tan(theta_s_e);
  let TAN_THETA_i_e = tan(theta_i_e);

  let omega_p = omega_s + omega_i; // spdc.pump.frequency();

  let n_p = spdc.pump.refractive_index(omega_p, &spdc.crystal_setup);
  let k_p = frequency_to_wavenumber(omega_p, n_p);
  let n_s = spdc.signal.refractive_index(omega_s, &spdc.crystal_setup);
  let n_i = spdc.idler.refractive_index(omega_i, &spdc.crystal_setup);

  // Counter-propagation requires a sign change in the wavenumber but
  // not for the free propagation constants ks_f and ki_f.
  let sign_ks = spdc.signal.direction().z.signum();
  let sign_ki = spdc.idler.direction().z.signum();
  let k_s = sign_ks * frequency_to_wavenumber(omega_s, n_s);
  let k_i = sign_ki * frequency_to_wavenumber(omega_i, n_i);

  let rho = tan(spdc.pump.walkoff_angle(&spdc.crystal_setup));

  // Now calculate the the coeficients that get repeatedly used.
  // This is from Karina's code.
  use dim::Abs;
  let ks_free = k_s.abs() / n_s;
  let ki_free = k_i.abs() / n_i;

  let Γ_s2 = -0.25 * Ws_SQ;
  let Γ_i2 = -0.25 * Wi_SQ;
  let Γ_s1 = Γ_s2 * SEC_2_THETA_s;
  let Γ_i1 = Γ_i2 * SEC_2_THETA_i;

  let Γ_s3 = -2. * ks_free * Γ_s1 * SIN_THETA_s_e / RAD;
  let Γ_i3 = -2. * ki_free * Γ_i1 * SIN_THETA_i_e / RAD;
  let Γ_s4 = -0.5 * ks_free * SIN_THETA_s_e * Γ_s3 / RAD;
  let Γ_i4 = -0.5 * ki_free * SIN_THETA_i_e * Γ_i3 / RAD;

  let zhs = ζs + hs * SIN_THETA_s_e;
  let zhi = ζi + hi * SIN_THETA_i_e;

  // zR is rayleigh range = 0.5 * k * W^2
  let zR_s = 0.5 * ks_free * Ws_SQ;
  let zR_i = 0.5 * ki_free * Wi_SQ;
  //
  // NOTE: v3 and v2 differ by a minus sign in the lambda terms
  // seems to be the version 3 used the non-complex-conjugate version of the mode
  // let Λ_s2 = (0.5 / ks_free) * zhs * RAD;
  // let Λ_i2 = (0.5 / ki_free) * zhi * RAD;
  let Λ_s2 = 0.25 * Ws_SQ * zhs / zR_s * RAD;
  let Λ_i2 = 0.25 * Wi_SQ * zhi / zR_i * RAD;
  let Λ_s1 = Λ_s2 * SEC_2_THETA_s;
  let Λ_i1 = Λ_i2 * SEC_2_THETA_i;
  // sec(a) tan(a) = sec^2(a) sin(a)
  // let Λ_s3 = -hs - zhs * SEC_2_THETA_s * SIN_THETA_s_e;
  // let Λ_i3 = -hi - zhi * SEC_2_THETA_i * SIN_THETA_i_e;
  // let Λ_s4 = 0.5 * ks_free * zhs * TAN_THETA_s_e.powi(2) / RAD - ks_free * ζs / RAD;
  // let Λ_i4 = 0.5 * ki_free * zhi * TAN_THETA_i_e.powi(2) / RAD - ki_free * ζi / RAD;
  let Λ_s3 = -hs - (0.5 * ks_free * Ws_SQ / zR_s) * zhs * SEC_2_THETA_s * SIN_THETA_s_e;
  let Λ_i3 = -hi - (0.5 * ki_free * Wi_SQ / zR_i) * zhi * SEC_2_THETA_i * SIN_THETA_i_e;
  let Λ_s4 =
    0.25 * (sq(ks_free) * Ws_SQ / zR_s) * zhs * TAN_THETA_s_e.powi(2) / RAD - ks_free * ζs / RAD;
  let Λ_i4 =
    0.25 * (sq(ki_free) * Wi_SQ / zR_i) * zhi * TAN_THETA_i_e.powi(2) / RAD - ki_free * ζi / RAD;
  // dbg!(
  //   zhs, zhi, Γ_s1, Γ_i1, Γ_s2, Γ_i2, Γ_s3, Γ_i3, Γ_s4, Γ_i4, Λ_s1, Λ_i1, Λ_s2, Λ_i2, Λ_s3, Λ_i3,
  //   Λ_s4, Λ_i4
  // );

  let zlks = 0.5 * (ζ / k_p - half_L / k_s) * RAD;
  let zlki = 0.5 * (ζ / k_p - half_L / k_i) * RAD;

  // The underscore terms lack the terms that depend on the integration parameter z
  let A1_ = complex_dim(-0.25 * Wx_SQ + Γ_s1, -Λ_s1);
  let A2_ = complex_dim(-0.25 * Wy_SQ + Γ_s2, -Λ_s2);
  let A3_ = complex_dim(-0.25 * Wx_SQ + Γ_i1, -Λ_i1);
  let A4_ = complex_dim(-0.25 * Wy_SQ + Γ_i2, -Λ_i2);
  let A5 = complex_dim(Γ_s3, -Λ_s3);
  let A7 = complex_dim(Γ_i3, -Λ_i3);
  let A8_ = real_dim(-0.5 * Wx_SQ);
  let A9_ = real_dim(-0.5 * Wy_SQ);
  let ksip = k_s + k_i + spdc.pp.k_eff();
  let A10_ = complex_dim(Γ_s4 + Γ_i4, -Λ_s4 - Λ_i4 + half_L * (k_p + ksip) / RAD);

  let Kps = 0.25 * L * (1. / k_p - 1. / k_s) * RAD;
  let Kpi = 0.25 * L * (1. / k_p - 1. / k_i) * RAD;
  let Lrho_by_2 = 0.5 * L * rho;

  // dbg!(A1_, A2_, A3_, A4_, A5, A7, A8_, A9_, A10_);

  let M2 = M * M;
  let A5sq = A5 * A5;

  // dbg!(zlks, zlki, Kps, Kpi);

  let fn_z = |z: f64| {
    let a = imag_dim(Kps * z - zlks);
    let b = imag_dim(Kpi * z - zlki);
    let A1 = A1_ - a;
    let A2 = A2_ - a;
    let A3 = A3_ - b;
    let A4 = A4_ - b;
    let c = imag_dim(Lrho_by_2 * (z + 1.));
    let A6 = c;
    let d = imag_dim((ζ - half_L * z) / k_p * RAD);
    let A8 = A8_ + d;
    let A9 = A9_ + d;
    let A10 = A10_ + imag_dim(half_L * z * (k_p - ksip) / RAD);

    let A6sq = A6 * A6;
    let A8sq = A8 * A8;
    let A9sq = A9 * A9;
    let invA1 = 1. / A1;
    let invA2 = 1. / A2;
    let denom1 = 4. * A1 * A3 - A8sq;
    let denom2 = 4. * A2 * A4 - A9sq;

    // dbg!(A1, A2);

    let numerator = (A10
      - 0.25
        * (invA1 * (A5sq + sq(-2. * A1 * A7 + A5 * A8) / denom1)
          + invA2 * A6sq * (1. + sq(-2. * A2 + A9) / denom2)))
      .exp();

    // Now deal with the denominator in the integral:
    // Sqrt[A1 A2 (-4 A3 + A8^2/A1) (-4 A4 + A9^2/A2)]
    let denom_units = M2 * M2 * M2 * M2;
    let denominator = (denom1 * denom2 / denom_units).sqrt();

    // dbg!(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10);

    // Take into account apodized crystals
    let pmzcoeff = spdc.pp.integration_constant(z, L);

    // dbg!(pmzcoeff, z, numerator, denominator);
    // Now calculate the full term in the integral.
    pmzcoeff * numerator / denominator
  };

  let result = 0.5 * integrator.integrate(fn_z, -1., 1.);

  PerMeter4::new(result)
}

// This implementation is based on Karina's FinalExpressions_Coincidence 2015-04-30.nb
#[allow(dead_code)]
#[allow(non_snake_case)]
fn phasematch_fiber_coupling_v3(
  omega_s: Frequency,
  omega_i: Frequency,
  spdc: &SPDC,
  integrator: Integrator,
) -> PerMeter4<Complex<f64>> {
  // crystal length
  let L = spdc.crystal_setup.length;

  let theta_s = spdc.signal.theta_internal();
  let phi_s = spdc.signal.phi();
  let theta_i = spdc.idler.theta_internal();
  let phi_i = spdc.idler.phi();
  let theta_s_e = spdc.signal.theta_external(&spdc.crystal_setup);
  let theta_i_e = spdc.idler.theta_external(&spdc.crystal_setup);

  let SIN_THETA_s_e = sin(theta_s_e);
  let SIN_THETA_i_e = sin(theta_i_e);
  let TAN_THETA_s_e = tan(theta_s_e);
  let TAN_THETA_i_e = tan(theta_i_e);
  let SIN_PHI_s = sin(phi_s);
  let SIN_PHI_i = sin(phi_i);
  let COS_PHI_s = cos(phi_s);
  let COS_PHI_i = cos(phi_i);
  let SEC_THETA_s_e = sec(theta_s_e);
  let SEC_THETA_i_e = sec(theta_i_e);

  let z0 = 0. * M; //put pump in middle of the crystal
  let ζ_s = spdc.signal_waist_position;
  let ζ_i = spdc.idler_waist_position;
  // Height of the collected spots from the z axis.
  let half_L = L / 2.;
  let h_s = half_L * tan(theta_s) * cos(phi_s);
  let h_i = half_L * tan(theta_i) * cos(phi_i);
  let d_s = half_L * tan(theta_s) * sin(phi_s);
  let d_i = half_L * tan(theta_i) * sin(phi_i);

  let Ws_SQ = spdc.signal.waist().x_by_y();
  let Wi_SQ = spdc.idler.waist().x_by_y();

  let Wx_SQ = sq(spdc.pump.waist().x);
  let Wy_SQ = sq(spdc.pump.waist().y);

  let COS_2_PHI_s = cos(phi_s).powi(2);
  let COS_2_PHI_i = cos(phi_i).powi(2);
  let SIN_2_PHI_s = sin(phi_s).powi(2);
  let SIN_2_PHI_i = sin(phi_i).powi(2);
  let SEC_2_THETA_s = sec(theta_s_e).powi(2);
  let SEC_2_THETA_i = sec(theta_i_e).powi(2);
  let TAN_2_THETA_s = tan(theta_s_e).powi(2);
  let TAN_2_THETA_i = tan(theta_i_e).powi(2);

  let omega_p = omega_s + omega_i; // spdc.pump.frequency();

  let n_p = spdc.pump.refractive_index(omega_p, &spdc.crystal_setup);
  let n_s = spdc.signal.refractive_index(omega_s, &spdc.crystal_setup);
  let n_i = spdc.idler.refractive_index(omega_i, &spdc.crystal_setup);

  // Counter-propagation requires a sign change in the wavenumber but
  // not for the free propagation constants ks_f and ki_f.
  let sign_ks = spdc.signal.direction().z.signum();
  let sign_ki = spdc.idler.direction().z.signum();

  let k_p = frequency_to_wavenumber(omega_p, n_p);
  let k_s = sign_ks * frequency_to_wavenumber(omega_s, n_s);
  let k_i = sign_ki * frequency_to_wavenumber(omega_i, n_i);

  let rho = tan(spdc.pump.walkoff_angle(&spdc.crystal_setup));

  let M2 = M * M; // meters squared

  use dim::Abs;
  let ks_free = k_s.abs() / n_s;
  let ki_free = k_i.abs() / n_i;

  let Γ_s1 = -0.25 * Ws_SQ * (COS_2_PHI_s * SEC_2_THETA_s + SIN_2_PHI_s);
  let Γ_i1 = -0.25 * Wi_SQ * (COS_2_PHI_i * SEC_2_THETA_i + SIN_2_PHI_i);

  let Γ_s2 = -0.25 * Ws_SQ * (SIN_2_PHI_s * SEC_2_THETA_s + COS_2_PHI_s);
  let Γ_i2 = -0.25 * Wi_SQ * (SIN_2_PHI_i * SEC_2_THETA_i + COS_2_PHI_i);

  let Γ_s3 = 0.5 * ks_free * Ws_SQ * COS_PHI_s * SEC_THETA_s_e * TAN_THETA_s_e / RAD;
  let Γ_i3 = 0.5 * ki_free * Wi_SQ * COS_PHI_i * SEC_THETA_i_e * TAN_THETA_i_e / RAD;

  let Γ_s4 = 0.5 * ks_free * Ws_SQ * SIN_PHI_s * SEC_THETA_s_e * TAN_THETA_s_e / RAD;
  let Γ_i4 = 0.5 * ki_free * Wi_SQ * SIN_PHI_i * SEC_THETA_i_e * TAN_THETA_i_e / RAD;

  let Γ_s5 = -0.25 * Ws_SQ * sin(2. * phi_s) * TAN_2_THETA_s;
  let Γ_i5 = -0.25 * Wi_SQ * sin(2. * phi_i) * TAN_2_THETA_i;

  let Γ_s6 = -0.25 * sq(ks_free) * Ws_SQ * TAN_2_THETA_s / RAD / RAD;
  let Γ_i6 = -0.25 * sq(ki_free) * Wi_SQ * TAN_2_THETA_i / RAD / RAD;

  // The lamda terms simplify quite a bit when you choose good intermediate variables:
  let H_s = SIN_THETA_s_e * h_s * COS_PHI_s;
  let H_i = SIN_THETA_i_e * h_i * COS_PHI_i;
  let D_s = SIN_THETA_s_e * d_s * SIN_PHI_s;
  let D_i = SIN_THETA_i_e * d_i * SIN_PHI_i;

  let HDZ_s = H_s + D_s + ζ_s;
  let HDZ_i = H_i + D_i + ζ_i;

  let SIN_COS_SEC_s = SIN_2_PHI_s + COS_2_PHI_s * SEC_2_THETA_s;
  let SIN_COS_SEC_i = SIN_2_PHI_i + COS_2_PHI_i * SEC_2_THETA_i;
  let Λ_s1 = (0.5 / ks_free) * SIN_COS_SEC_s * HDZ_s * RAD;
  let Λ_i1 = (0.5 / ki_free) * SIN_COS_SEC_i * HDZ_i * RAD;

  let COS_SIN_SEC_s = COS_2_PHI_s + SIN_2_PHI_s * SEC_2_THETA_s;
  let COS_SIN_SEC_i = COS_2_PHI_i + SIN_2_PHI_i * SEC_2_THETA_i;
  let Λ_s2 = (0.5 / ks_free) * COS_SIN_SEC_s * HDZ_s * RAD;
  let Λ_i2 = (0.5 / ki_free) * COS_SIN_SEC_i * HDZ_i * RAD;

  let Λ_s3 = -SIN_COS_SEC_s * h_s - COS_PHI_s * TAN_THETA_s_e * SEC_THETA_s_e * (D_s + ζ_s);
  let Λ_i3 = -SIN_COS_SEC_i * h_i - COS_PHI_i * TAN_THETA_i_e * SEC_THETA_s_e * (D_i + ζ_i);

  let Λ_s4 = -SIN_PHI_s * TAN_THETA_s_e * SEC_THETA_s_e * (H_s + ζ_s) - COS_SIN_SEC_s * d_s;
  let Λ_i4 = -SIN_PHI_i * TAN_THETA_i_e * SEC_THETA_s_e * (H_i + ζ_i) - COS_SIN_SEC_i * d_i;

  let SIN_COS_TAN_s = SIN_PHI_s * COS_PHI_s * TAN_2_THETA_s;
  let SIN_COS_TAN_i = SIN_PHI_i * COS_PHI_i * TAN_2_THETA_i;
  let Λ_s5 = (1. / ks_free) * SIN_COS_TAN_s * HDZ_s * RAD;
  let Λ_i5 = (1. / ki_free) * SIN_COS_TAN_i * HDZ_i * RAD;

  let Λ_s6 = (0.5 * ks_free) * (TAN_2_THETA_s * HDZ_s - 2. * ζ_s) / RAD;
  let Λ_i6 = (0.5 * ki_free) * (TAN_2_THETA_i * HDZ_i - 2. * ζ_i) / RAD;

  let zlks = 0.5 * (z0 / k_p - half_L / k_s) * RAD;
  let zlki = 0.5 * (z0 / k_p - half_L / k_i) * RAD;

  // The omega terms lack the terms that depend on the integration parameter z
  let Ω1_ = complex_dim(-0.25 * Wx_SQ + Γ_s1, zlks + Λ_s1);
  let Ω2_ = complex_dim(-0.25 * Wy_SQ + Γ_s2, zlks + Λ_s2);
  let Ω3_ = complex_dim(-0.25 * Wx_SQ + Γ_i1, zlki + Λ_i1);
  let Ω4_ = complex_dim(-0.25 * Wy_SQ + Γ_i2, zlki + Λ_i2);
  let Ω5 = complex_dim(Γ_s3, Λ_s3);
  let Ω6 = complex_dim(Γ_i3, Λ_i3);
  let Ω7_ = complex_dim(Γ_s4, Λ_s4);
  let Ω8_ = complex_dim(Γ_i4, Λ_i4);
  let Ω9 = complex_dim(Γ_s5, Λ_s5);
  let Ω10 = complex_dim(Γ_i5, Λ_i5);
  let Ω11_ = real_dim(-0.5 * Wx_SQ);
  let Ω12_ = real_dim(-0.5 * Wy_SQ);
  let ksip = k_s + k_i + spdc.pp.k_eff();
  let Ω13_ = complex_dim(Γ_s6 + Γ_i6, Λ_s6 + Λ_i6 + half_L * (k_p + ksip) / RAD);

  let Kps = 0.25 * L * (1. / k_p - 1. / k_s) * RAD;
  let Kpi = 0.25 * L * (1. / k_p - 1. / k_i) * RAD;
  let Lrho_by_2 = half_L * rho;

  let fn_z = |z: f64| {
    let a = imag_dim(Kps * z);
    let b = imag_dim(Kpi * z);
    let Ω1 = Ω1_ - a;
    let Ω2 = Ω2_ - a;
    let Ω3 = Ω3_ - b;
    let Ω4 = Ω4_ - b;
    let c = imag_dim(Lrho_by_2 * (z + 1.));
    let Ω7 = Ω7_ + c;
    let Ω8 = Ω8_ + c;
    let d = imag_dim((z0 - half_L * z) / k_p * RAD);
    let Ω11 = Ω11_ + d;
    let Ω12 = Ω12_ + d;
    let Ω13 = Ω13_ + imag_dim(half_L * z * (k_p - ksip) / RAD);

    // -4 Ω1 Ω102 Ω2 + Ω112 - 4 Ω1 Ω3 Ω122 - 4 Ω2 Ω4 - 2 Ω10 Ω11 Ω12 Ω9 + Ω102 - 4 Ω3 Ω4 Ω92
    let Δ = -4. * Ω1 * sq(Ω10) * Ω2 + (sq(Ω11) - 4. * Ω1 * Ω3) * (sq(Ω12) - 4. * Ω2 * Ω4)
      - 2. * Ω10 * Ω11 * Ω12 * Ω9
      + (sq(Ω10) - 4. * Ω3 * Ω4) * sq(Ω9);
    // Ω112 -4Ω1Ω3(Ω12Ω7-2Ω2Ω8)+Ω12(2Ω3Ω5-Ω11Ω6)Ω9-2Ω3Ω8Ω92 +Ω102Ω11Ω2Ω5-4Ω1Ω2Ω6-Ω11Ω7Ω9+Ω6Ω922
    let N = sq(
      (sq(Ω11) - 4. * Ω1 * Ω3) * (Ω12 * Ω7 - 2. * Ω2 * Ω8) + Ω12 * (2. * Ω3 * Ω5 - Ω11 * Ω6) * Ω9
        - 2. * Ω3 * Ω8 * sq(Ω9)
        + Ω10 * (2. * Ω11 * Ω2 * Ω5 - 4. * Ω1 * Ω2 * Ω6 - Ω11 * Ω7 * Ω9 + Ω6 * sq(Ω9)),
    );
    // 4 Ω112 Ω13 Ω2 - 16 Ω1 Ω13 Ω2 Ω3 + 4 Ω2 Ω3 Ω52 - 4 Ω11 Ω2 Ω5 Ω6 + 4 Ω1 Ω2 Ω62 - Ω112 Ω72 +4Ω1Ω3Ω72 -4Ω3Ω5Ω7Ω9+2Ω11Ω6Ω7Ω9+4Ω13Ω3Ω92 -Ω62 Ω92
    let m = 4. * sq(Ω11) * Ω13 * Ω2 - 16. * Ω1 * Ω13 * Ω2 * Ω3 + 4. * Ω2 * Ω3 * sq(Ω5)
      - 4. * Ω11 * Ω2 * Ω5 * Ω6
      + 4. * Ω1 * Ω2 * sq(Ω6)
      - sq(Ω11) * sq(Ω7)
      + 4. * Ω1 * Ω3 * sq(Ω7)
      - 4. * Ω3 * Ω5 * Ω7 * Ω9
      + 2. * Ω11 * Ω6 * Ω7 * Ω9
      + 4. * Ω13 * Ω3 * sq(Ω9)
      - sq(Ω6) * sq(Ω9);
    // 4Ω112 Ω2+Ω3-4Ω1Ω2+Ω92
    let P = 4. * (sq(Ω11) * Ω2 + Ω3 * (-4. * Ω1 * Ω2 + sq(Ω9)));

    let denom_units = M2 * M2 * M2 * M2;
    let ψ = ((m + N / Δ) / P).exp() / (Δ / denom_units).sqrt();

    // Take into account apodized crystals
    let pmzcoeff = spdc.pp.integration_constant(z, L);

    // dbg!(pmzcoeff, z, numerator, denominator);
    // Now calculate the full term in the integral.
    ψ * pmzcoeff
  };

  // let integrator = SimpsonIntegration::new(fn_z);
  // let result = integrator.integrate(
  //   -1.,
  //   1.,
  //   steps.unwrap_or_else(|| integration_steps_best_guess(L)),
  // );

  let result = integrator.integrate(fn_z, -1., 1.);
  PerMeter4::new(result)
}

#[cfg(test)]
mod tests {
  use super::*;
  extern crate float_cmp;
  use crate::utils::testing::assert_nearly_equal;
  use dim::{
    f64prefixes::{MICRO, NANO},
    ucum::DEG,
    Dimensioned,
  };
  use utils::{frequency_to_vacuum_wavelength, testing::testing_props};

  #[test]
  fn phasematch_test() {
    let spdc = crate::utils::testing::testing_props(false);

    let expected = 132056158795971.14;
    let actual = phasematch_fiber_coupling(
      spdc.signal.frequency(),
      spdc.idler.frequency(),
      &spdc,
      Integrator::Simpson { divs: 50 },
    );

    assert_nearly_equal!(
      "coincidences magnitude",
      actual.value_unsafe().norm(),
      expected,
      0.1
    );
  }

  #[test]
  fn phasematch_fiber_coupling_test() {
    let mut spdc = SPDC::default();
    // spdc.signal.set_from_external_theta(3. * DEG, &spdc.crystal_setup);
    // spdc.signal.set_angles(0. * RAD, 0.03253866877817829 * RAD);
    spdc = spdc.try_as_optimum().unwrap();

    // println!("spdc: {:#?}", spdc);
    let jsa_units = JSAUnits::new(1.);
    let amp = *(phasematch_fiber_coupling(
      spdc.signal.frequency(),
      spdc.idler.frequency(),
      &spdc,
      Integrator::Simpson { divs: 50 },
    ) / jsa_units);

    let actual = amp;
    let expected =
      *(phasematch_gaussian(spdc.signal.frequency(), spdc.idler.frequency(), &spdc) / jsa_units);

    // NOTE: this is not a great test anymore
    let accept_diff = 1e-4;

    assert_nearly_equal!(actual.norm(), expected.norm(), accept_diff);
  }

  #[test]
  fn phasematch_fiber_coupling_pp_test() {
    let spdc = testing_props(true);

    dbg!(&spdc);

    let jsa_units = JSAUnits::new(1.);
    let amp = *(phasematch_fiber_coupling(
      spdc.signal.frequency(),
      spdc.idler.frequency(),
      &spdc,
      Integrator::Simpson { divs: 50 },
    ) / jsa_units);

    let actual = amp;
    // let expected = Complex::new(-243675412686457.94, 411264607672255.2);
    // Before refactor
    let expected = Complex::new(-507267866308567.25, -590079455430265.8);
    dbg!(actual);

    let accept_diff = 1e-16;

    assert_nearly_equal!("norm", actual.norm(), expected.norm(), accept_diff);

    assert_nearly_equal!("arg", actual.arg(), expected.arg(), accept_diff);
  }

  #[test]
  fn compare_version1_version2() {
    let json = serde_json::json!({
      "crystal": {
        "kind": "KTP",
        "pm_type": "e->eo",
        "phi_deg": 3,
        "theta_deg": 0,
        "length_um": 2000,
        "temperature_c": 20
      },
      "pump": {
        "wavelength_nm": 775,
        "waist_um": 90,
        "bandwidth_nm": 5.35,
        "average_power_mw": 1
      },
      "signal": {
        "wavelength_nm": 1550,
        "phi_deg": 0,
        "theta_external_deg": 1,
        "waist_um": 10,
        "waist_position_um": "auto"
      },
      "idler": "auto",
      "pp": "auto",
      "deff_pm_per_volt": 1.
    });

    let config: SPDCConfig = serde_json::from_value(json).expect("Could not unwrap json");
    let mut spdc = config
      .try_as_spdc()
      .expect("Could not convert to SPDC instance");

    fn compare(spdc: &SPDC) {
      let wavelengths = WavelengthSpace::new(
        (1510. * NANO * M, 1590. * NANO * M, 10),
        (1510. * NANO * M, 1590. * NANO * M, 10),
      );
      for (ws, wi) in wavelengths.into_signal_idler_iterator() {
        let old = *(phasematch_fiber_coupling(ws, wi, spdc, Integrator::Simpson { divs: 5 })
          / JSAUnits::new(1.));
        let new = *(phasematch_fiber_coupling2(ws, wi, spdc, Integrator::Simpson { divs: 5 })
          / JSAUnits::new(1.));
        let ls = frequency_to_vacuum_wavelength(ws);
        let li = frequency_to_vacuum_wavelength(wi);
        assert_nearly_equal!(
          format!("norm at {}, {}, {:#?}", ls, li, spdc),
          new.norm(),
          old.norm(),
          1e-6
        );
      }
    }

    for len in 1..10 {
      dbg!(len);
      spdc.crystal_setup.length = len as f64 * 1000. * MICRO * M;
      for theta in 0..10 {
        dbg!(theta);
        spdc.crystal_setup.counter_propagation = false;
        spdc.signal.set_angles(0. * DEG, theta as f64 * DEG);
        spdc.assign_optimum_idler().unwrap();
        spdc.assign_optimum_crystal_theta();
        compare(&spdc);
        spdc.signal.set_angles(180. * DEG, theta as f64 * DEG);
        spdc.assign_optimum_idler().unwrap();
        spdc.assign_optimum_crystal_theta();
        compare(&spdc);
        spdc.crystal_setup.counter_propagation = true;
        spdc.signal.set_angles(0. * DEG, theta as f64 * DEG);
        spdc.assign_optimum_idler().unwrap();
        spdc.assign_optimum_crystal_theta();
        compare(&spdc);
        spdc
          .signal
          .set_angles(180. * DEG, (180. - theta as f64) * DEG);
        spdc.assign_optimum_idler().unwrap();
        spdc.assign_optimum_crystal_theta();
        compare(&spdc);
      }
    }

    spdc.crystal_setup.theta = 0. * DEG;

    for len in 1..10 {
      dbg!(len);
      spdc.crystal_setup.length = len as f64 * 1000. * MICRO * M;
      for theta in 0..10 {
        dbg!(theta);
        spdc.crystal_setup.counter_propagation = false;
        spdc.signal.set_angles(0. * DEG, theta as f64 * DEG);
        spdc.assign_optimum_idler().unwrap();
        spdc.assign_optimum_periodic_poling().unwrap();
        compare(&spdc);
        spdc.signal.set_angles(180. * DEG, theta as f64 * DEG);
        spdc.assign_optimum_idler().unwrap();
        spdc.assign_optimum_periodic_poling().unwrap();
        compare(&spdc);
        spdc.crystal_setup.counter_propagation = true;
        spdc.signal.set_angles(0. * DEG, theta as f64 * DEG);
        spdc.assign_optimum_idler().unwrap();
        spdc.assign_optimum_periodic_poling().unwrap();
        compare(&spdc);
        spdc
          .signal
          .set_angles(180. * DEG, (180. - theta as f64) * DEG);
        spdc.assign_optimum_idler().unwrap();
        spdc.assign_optimum_periodic_poling().unwrap();
        compare(&spdc);
      }
    }
  }

  // #[test]
  // fn compare_version1_version3() {
  //   let json = serde_json::json!({
  //     "crystal": {
  //       "kind": "KTP",
  //       "pm_type": "e->eo",
  //       "phi_deg": 3,
  //       "theta_deg": 0,
  //       "length_um": 2000,
  //       "temperature_c": 20
  //     },
  //     "pump": {
  //       "wavelength_nm": 775,
  //       "waist_um": 90,
  //       "bandwidth_nm": 5.35,
  //       "average_power_mw": 1
  //     },
  //     "signal": {
  //       "wavelength_nm": 1550,
  //       "phi_deg": 0,
  //       "theta_external_deg": 1,
  //       "waist_um": 10,
  //       "waist_position_um": "auto"
  //     },
  //     "idler": "auto",
  //     "pp": "auto",
  //     "deff_pm_per_volt": 1.
  //   });

  //   let config: SPDCConfig = serde_json::from_value(json).expect("Could not unwrap json");
  //   let mut spdc = config
  //     .try_as_spdc()
  //     .expect("Could not convert to SPDC instance");

  //   fn compare(spdc: &SPDC) {
  //     let wavelengths = WavelengthSpace::new(
  //       (1510. * NANO * M, 1590. * NANO * M, 10),
  //       (1510. * NANO * M, 1590. * NANO * M, 10),
  //     );
  //     for (ws, wi) in wavelengths.into_signal_idler_iterator() {
  //       let old = *(phasematch_fiber_coupling(ws, wi, &spdc, Integrator::Simpson { divs: 5 })
  //         / JSAUnits::new(1.));
  //       let new = *(phasematch_fiber_coupling_v3(ws, wi, &spdc, Integrator::Simpson { divs: 5 })
  //         / JSAUnits::new(1.));
  //       let ls = frequency_to_vacuum_wavelength(ws);
  //       let li = frequency_to_vacuum_wavelength(wi);
  //       assert_nearly_equal!(
  //         format!("norm at {}, {}, {:#?}", ls, li, spdc),
  //         new.norm(),
  //         old.norm(),
  //         50.
  //       );
  //     }
  //   }

  //   for len in 1..10 {
  //     dbg!(len);
  //     spdc.crystal_setup.length = len as f64 * 1000. * MICRO * M;
  //     for theta in 0..10 {
  //       dbg!(theta);
  //       spdc.crystal_setup.counter_propagation = false;
  //       spdc.signal.set_angles(0. * DEG, theta as f64 * DEG);
  //       spdc.assign_optimum_idler().unwrap();
  //       spdc.assign_optimum_crystal_theta();
  //       compare(&spdc);
  //       spdc.signal.set_angles(180. * DEG, theta as f64 * DEG);
  //       spdc.assign_optimum_idler().unwrap();
  //       spdc.assign_optimum_crystal_theta();
  //       compare(&spdc);
  //       spdc.crystal_setup.counter_propagation = true;
  //       spdc.signal.set_angles(0. * DEG, theta as f64 * DEG);
  //       spdc.assign_optimum_idler().unwrap();
  //       spdc.assign_optimum_crystal_theta();
  //       compare(&spdc);
  //       spdc
  //         .signal
  //         .set_angles(180. * DEG, (180. - theta as f64) * DEG);
  //       spdc.assign_optimum_idler().unwrap();
  //       spdc.assign_optimum_crystal_theta();
  //       compare(&spdc);
  //     }
  //   }

  //   spdc.crystal_setup.theta = 0. * DEG;

  //   for len in 1..10 {
  //     dbg!(len);
  //     spdc.crystal_setup.length = len as f64 * 1000. * MICRO * M;
  //     for theta in 0..10 {
  //       dbg!(theta);
  //       spdc.crystal_setup.counter_propagation = false;
  //       spdc.signal.set_angles(0. * DEG, theta as f64 * DEG);
  //       spdc.assign_optimum_idler().unwrap();
  //       spdc.assign_optimum_periodic_poling().unwrap();
  //       compare(&spdc);
  //       spdc.signal.set_angles(180. * DEG, theta as f64 * DEG);
  //       spdc.assign_optimum_idler().unwrap();
  //       spdc.assign_optimum_periodic_poling().unwrap();
  //       compare(&spdc);
  //       spdc.crystal_setup.counter_propagation = true;
  //       spdc.signal.set_angles(0. * DEG, theta as f64 * DEG);
  //       spdc.assign_optimum_idler().unwrap();
  //       spdc.assign_optimum_periodic_poling().unwrap();
  //       compare(&spdc);
  //       spdc
  //         .signal
  //         .set_angles(180. * DEG, (180. - theta as f64) * DEG);
  //       spdc.assign_optimum_idler().unwrap();
  //       spdc.assign_optimum_periodic_poling().unwrap();
  //       compare(&spdc);
  //     }
  //   }
  // }
}
