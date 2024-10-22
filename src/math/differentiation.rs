//! Modified from <https://github.com/b52/optimization-rust/>

/// Compute the gradient (or derrivative) of function
///
/// Uses simple one step forward finite difference with step width `h = √εx`.
///
/// # Examples
///
/// ```
/// # use spdcalc::math::*;
/// let square = |x : &[f64]| x[0] * x[0];
///
/// assert!(gradient_at(square, &[0.0])[0] < 1.0e-3);
/// assert!(gradient_at(square, &[1.0])[0] > 1.0);
/// assert!(gradient_at(square, &[-1.0])[0] < 1.0);
/// ```
pub fn gradient_at<T>(func: impl Fn(&[f64]) -> f64, position: T) -> Vec<f64>
where
  T: AsRef<[f64]>,
{
  let mut x: Vec<_> = position.as_ref().to_vec();

  position
    .as_ref()
    .iter()
    .cloned()
    .enumerate()
    .map(|(i, x_i)| {
      let h = if x_i == 0.0 {
        f64::EPSILON.powf(1. / 3.)
        // f64::EPSILON * 1.0e10
      } else {
        f64::EPSILON.powf(1. / 3.) * x_i.abs()
      };

      assert!(h.is_finite(), "Derivative 'h' is infinite!");

      x[i] = x_i + h;

      let forward = func(&x);

      x[i] = x_i - h;

      let backward = func(&x);

      x[i] = x_i;

      let d_i = 0.5 * (forward - backward) / h;

      assert!(d_i.is_finite(), "Derivative is infinite!");

      d_i
    })
    .collect()
}

/// Compute the derrivative of function at a given position
pub fn derivative_at(func: impl Fn(f64) -> f64, position: f64) -> f64 {
  gradient_at(|x: &[f64]| func(x[0]), [position])[0]
}

#[cfg(test)]
mod tests {
  use super::*;
  extern crate float_cmp;
  use float_cmp::*;

  #[test]
  fn derrivative_test() {
    let func = |x: f64| x.sin();
    let func_prime = |x: f64| x.cos();

    let x = 0.4;
    let actual = derivative_at(func, x);
    let expected = func_prime(x);

    assert!(
      approx_eq!(f64, actual, expected, ulps = 2, epsilon = 1e-8),
      "actual: {}, expected: {}",
      actual,
      expected
    );
  }
}
