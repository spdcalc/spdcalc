//! Mathematical helpers

mod differentiation;
use crate::constants::{PI, TWO_PI};
use crate::{dim::ucum::RAD, Angle};
pub use differentiation::*;

mod integration;
pub use integration::*;

mod nelder_mead;
pub use self::nelder_mead::*;

mod schmidt;
pub use self::schmidt::*;

mod interpolation;
pub use interpolation::*;

lazy_static::lazy_static! {
  static ref FWHM_OVER_WAIST :f64 = f64::sqrt(2. * f64::ln(2.));
}

/// The sin function for dimensioned angles
pub fn sin(a: Angle) -> f64 {
  (a / RAD).sin()
}

/// The cos function for dimensioned angles
pub fn cos(a: Angle) -> f64 {
  (a / RAD).cos()
}
/// The tan function for dimensioned angles
pub fn tan(a: Angle) -> f64 {
  (a / RAD).tan()
}
/// The sec function for dimensioned angles
pub fn sec(a: Angle) -> f64 {
  1. / cos(a)
}
/// The csc function for dimensioned angles
pub fn csc(a: Angle) -> f64 {
  1. / sin(a)
}
/// The cot function for dimensioned angles
pub fn cot(a: Angle) -> f64 {
  1. / tan(a)
}

/// Standard sinc function `sinc(x) = sin(x) / x`
pub fn sinc(x: Angle) -> f64 {
  if x == 0. * RAD {
    1.
  } else {
    sin(x) / *(x / RAD)
  }
}

/// Square a value
pub fn sq<T>(n: T) -> <T as std::ops::Mul>::Output
where
  T: std::ops::Mul + Copy,
{
  n * n
}

/// Normalize an angle to [0, 2π)
pub fn normalize_angle(ang: Angle) -> Angle {
  let val = ang / RAD;
  val.rem_euclid(TWO_PI) * RAD
}

/// Normalize an angle to (-π, π]
pub fn normalize_angle_signed(ang: Angle) -> Angle {
  let val = ang / RAD;
  let rem = val.rem_euclid(TWO_PI);
  if rem > PI {
    (rem - TWO_PI) * RAD
  } else {
    rem * RAD
  }
}

/// Simple implementation of linear interpolation
pub fn lerp<T>(
  first: T,
  second: T,
  t: f64,
) -> <<T as std::ops::Mul<f64>>::Output as std::ops::Add>::Output
where
  T: std::ops::Mul<f64>,
  <T as std::ops::Mul<f64>>::Output: std::ops::Add<<T as std::ops::Mul<f64>>::Output>,
{
  first * (1. - t) + second * t
}

// http://mathworld.wolfram.com/GaussianFunction.html
/// FWHM / sigma = 2 * sqrt(2 * ln(2))
pub fn fwhm_to_sigma<T>(fwhm: T) -> <T as std::ops::Div<f64>>::Output
where
  T: std::ops::Div<f64>,
{
  fwhm / (2. * *FWHM_OVER_WAIST)
}

/// FWHM to 1/e^2 width
pub fn fwhm_to_waist<T>(fwhm: T) -> <T as std::ops::Div<f64>>::Output
where
  T: std::ops::Div<f64>,
{
  fwhm / *FWHM_OVER_WAIST
}

/// FWHM to 1/e^2 width
pub fn waist_to_fwhm<T>(w: T) -> <T as std::ops::Mul<f64>>::Output
where
  T: std::ops::Mul<f64>,
{
  w * *FWHM_OVER_WAIST
}

/// Round a number to a certain number of significant figures
pub fn sigfigs(x: f64, decimals: u8) -> f64 {
  let y = 10u32.pow(decimals as u32) as f64;
  (x * y).round() / y
}

#[cfg(test)]
mod test {
  use super::*;

  #[test]
  fn test_angle_normalizer() {
    assert_approx_eq!(f64, *(normalize_angle(2. * RAD) / RAD), 2.);
    assert_approx_eq!(f64, *(normalize_angle(2. * RAD + TWO_PI * RAD) / RAD), 2.);
    assert_approx_eq!(f64, *(normalize_angle(-2. * RAD) / RAD), TWO_PI - 2.);
    assert_approx_eq!(
      f64,
      *(normalize_angle(-2. * RAD - TWO_PI * RAD) / RAD),
      TWO_PI - 2.
    );
  }

  #[test]
  fn test_angle_normalizer_signed() {
    assert_approx_eq!(f64, *(normalize_angle_signed(1. * RAD) / RAD), 1.);
    assert_approx_eq!(
      f64,
      *(normalize_angle_signed(1. * RAD + TWO_PI * RAD) / RAD),
      1.
    );
    assert_approx_eq!(f64, *(normalize_angle_signed(-1. * RAD) / RAD), -1.);
    assert_approx_eq!(
      f64,
      *(normalize_angle_signed(-1. * RAD - TWO_PI * RAD) / RAD),
      -1.
    );
    assert_approx_eq!(
      f64,
      *(normalize_angle_signed(1. * RAD + PI * RAD) / RAD),
      1. - PI
    );
    assert_approx_eq!(
      f64,
      *(normalize_angle_signed(-1. * RAD - PI * RAD) / RAD),
      -1. + PI
    );
  }
}
