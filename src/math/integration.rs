use crate::utils::Steps;
use crate::Complex;
use num::Integer;
use quad_rs::Integrable;
use rayon::prelude::*;

struct Problem<F, Y>(pub F)
where
  F: Fn(Complex<f64>) -> Y;

impl<F, Y> Integrable for Problem<F, Y>
where
  Y: Into<Complex<f64>>,
  F: Fn(Complex<f64>) -> Y,
{
  type Input = Complex<f64>;
  type Output = Complex<f64>;
  fn integrand(
    &self,
    input: &Self::Input,
  ) -> Result<Self::Output, quad_rs::EvaluationError<Self::Input>> {
    Ok((self.0)(*input).into())
  }
}

/// Various integration methods
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(tag = "method")]
pub enum Integrator {
  /// Simpson integration with specified divisions
  Simpson {
    /// Divisions (slices) to use for integration
    divs: usize,
  },
  /// Adaptive Simpson integration
  AdaptiveSimpson {
    /// Relative tolerance
    tolerance: f64,
    /// Maximum depth
    max_depth: usize,
  },
  /// Gauss-Konrod integration
  GaussKonrod {
    /// Relative tolerance
    tolerance: f64,
    /// Maximum depth
    max_depth: usize,
  },
  /// Gauss-Legendre Integration
  GaussLegendre {
    /// Degree parameter
    degree: usize,
  },
  /// Clenshaw-Curtis integration
  ClenshawCurtis {
    /// Relative tolerance
    tolerance: f64,
  },
}

impl Default for Integrator {
  fn default() -> Self {
    Integrator::Simpson { divs: 50 }
  }
}

impl Integrator {
  /// Calculate the integral of function from x -> [a, b]
  pub fn integrate<F, Y>(&self, func: F, a: f64, b: f64) -> Complex<f64>
  where
    F: Send + Sync + Fn(f64) -> Y,
    Y: Into<Complex<f64>>,
  {
    match self {
      Integrator::Simpson { divs } => simpson(func, a, b, *divs),
      Integrator::AdaptiveSimpson {
        tolerance,
        max_depth,
      } => simpson_adaptive(&func, a, b, *tolerance, *max_depth),
      Integrator::GaussKonrod {
        tolerance,
        max_depth,
      } => {
        let a = Complex::new(a, 0.);
        let b = Complex::new(b, 0.);
        let integrator = quad_rs::Integrator::default()
          .relative_tolerance(*tolerance)
          .with_maximum_iter(*max_depth);
        integrator
          .integrate(Problem(|z| func(z.re)), a..b)
          .unwrap()
          .result
          .result
          .unwrap()
      }
      Integrator::ClenshawCurtis { tolerance } => {
        let re =
          quadrature::clenshaw_curtis::integrate(|z| func(z).into().re, a, b, *tolerance).integral;
        let im =
          quadrature::clenshaw_curtis::integrate(|z| func(z).into().im, a, b, *tolerance).integral;
        Complex::new(re, im)
      }
      Integrator::GaussLegendre { degree } => {
        let quad = gauss_quad::GaussLegendre::new(*degree.max(&2)).unwrap();
        let re = quad.integrate(a, b, |z| func(z).into().re);
        let im = quad.integrate(a, b, |z| func(z).into().im);
        Complex::new(re, im)
      }
    }
  }

  /// Perform a 2d integration of specified function over x -> [a, b] and y -> [c, d]
  pub fn integrate2d<F, Y>(&self, func: F, a: f64, b: f64, c: f64, d: f64) -> Complex<f64>
  where
    F: Send + Sync + Fn(f64, f64) -> Y,
    Y: Into<Complex<f64>>,
  {
    match self {
      Integrator::Simpson { divs } => simpson2d(func, a, b, c, d, *divs),
      Integrator::AdaptiveSimpson {
        tolerance,
        max_depth,
      } => simpson_adaptive_2d(&func, a, b, c, d, *tolerance, *max_depth),
      Integrator::GaussKonrod {
        tolerance,
        max_depth,
      } => {
        let a = Complex::new(a, 0.);
        let b = Complex::new(b, 0.);
        let c = Complex::new(c, 0.);
        let d = Complex::new(d, 0.);
        let integrator = quad_rs::Integrator::default()
          .relative_tolerance(*tolerance)
          .with_maximum_iter(*max_depth);
        integrator
          .integrate(
            Problem(|z: Complex<f64>| {
              quad_rs::Integrator::default()
                .relative_tolerance(*tolerance)
                .with_maximum_iter(*max_depth)
                .integrate(Problem(|w: Complex<f64>| func(z.re, w.re)), c..d)
                .unwrap()
                .result
                .result
                .unwrap()
            }),
            a..b,
          )
          .unwrap()
          .result
          .result
          .unwrap()
      }
      Integrator::ClenshawCurtis { tolerance } => {
        let re = quadrature::clenshaw_curtis::integrate(
          |z| {
            quadrature::clenshaw_curtis::integrate(|w| func(z, w).into().re, c, d, *tolerance)
              .integral
          },
          a,
          b,
          *tolerance,
        )
        .integral;
        let im = quadrature::clenshaw_curtis::integrate(
          |z| {
            quadrature::clenshaw_curtis::integrate(|w| func(z, w).into().im, c, d, *tolerance)
              .integral
          },
          a,
          b,
          *tolerance,
        )
        .integral;
        Complex::new(re, im)
      }
      Integrator::GaussLegendre { degree } => {
        let quad = gauss_quad::GaussLegendre::new(*degree.max(&2)).unwrap();
        let re = quad.integrate(a, b, |z| quad.integrate(c, d, |w| func(z, w).into().re));
        let im = quad.integrate(a, b, |z| quad.integrate(c, d, |w| func(z, w).into().im));
        Complex::new(re, im)
      }
    }
  }
}

/// Get simpson weight for index
fn get_simpson_weight(n: usize, divs: usize) -> f64 {
  // n divs of...
  // 1, 4, 2, 4, 2, ..., 4, 1
  if n == 0 || n == divs {
    1.
  } else if n.is_odd() {
    4.
  } else {
    2.
  }
}

// /// Integrator that implements Simpson's rule
// pub struct SimpsonIntegration<F: Fn(f64) -> T, T> {
//   function: F,
// }

// impl<F: Fn(f64) -> T, T> SimpsonIntegration<F, T>
// where
//   T: Zero + std::ops::Mul<f64, Output = T> + std::ops::Add<T, Output = T>,
// {
//   /// Creates a new integrable function by using the supplied `Fn(f64) -> T` in
//   /// combination with numeric integration via simpson's rule to find the integral.
//   pub fn new(function: F) -> Self {
//     SimpsonIntegration { function }
//   }

//   /// Get simpson weight for index
//   pub fn get_weight(n: usize, divs: usize) -> f64 {
//     get_simpson_weight(n, divs)
//   }

//   /// Numerically integrate from `a` to `b`, in `divs` divisions
//   pub fn integrate(&self, a: f64, b: f64, divs: usize) -> T {
//     let divs = divs + divs % 2 - 2; // nearest even
//     assert!(divs >= 4, "Steps too low");

//     let dx = (b - a) / (divs as f64);

//     let result = (0..=divs)
//       .map(|n| Self::get_weight(n, divs))
//       .enumerate()
//       .fold(T::zero(), |acc, (i, a_n)| {
//         let x = a + (i as f64) * dx;

//         acc + (self.function)(x) * a_n
//       });

//     result * (dx / 3.)
//   }
// }

/// Simpson integration of function from x -> [a, b] using `divs` slices
pub fn simpson<F, Y>(func: F, a: f64, b: f64, divs: usize) -> Complex<f64>
where
  F: Send + Sync + Fn(f64) -> Y,
  Y: Into<Complex<f64>>,
{
  let divs = divs + divs % 2 - 2; // nearest even
  assert!(divs >= 4, "Steps too low");
  let dx = (b - a) / (divs as f64);

  let intg = |(i, a_n)| {
    let x = a + (i as f64) * dx;

    func(x).into() * a_n
  };

  // Parallelize if divs is large
  let result: Complex<f64> = if divs < 128 {
    (0..=divs)
      .map(|n| (n, get_simpson_weight(n, divs)))
      .map(intg)
      .sum()
  } else {
    (0..=divs)
      .into_par_iter()
      .map(|n| (n, get_simpson_weight(n, divs)))
      .map(intg)
      .sum()
  };

  result * (dx / 3.)
}

/// Two dimensional simpson integration in x -> [a, b] and y -> [c, d] with `divs` slices on each axis
pub fn simpson2d<F, Y>(func: F, ax: f64, bx: f64, ay: f64, by: f64, divs: usize) -> Complex<f64>
where
  F: Send + Sync + Fn(f64, f64) -> Y,
  Y: Into<Complex<f64>>,
{
  assert!(divs.is_even());
  assert!(divs >= 4);

  let steps = divs + 1;
  let dx = (bx - ax) / (divs as f64);
  let dy = (by - ay) / (divs as f64);

  let result: Complex<f64> = Steps(ay, by, steps)
    .into_par_iter()
    .enumerate()
    .map(|(ny, y)| {
      let sy: Complex<f64> = Steps(ax, bx, steps)
        .into_par_iter()
        .enumerate()
        .map(|(nx, x)| {
          let a_n = get_simpson_weight(nx, divs);
          func(x, y).into() * a_n
        })
        .sum();

      let a_n = get_simpson_weight(ny, divs);
      sy * a_n
    })
    .sum();

  result * (dx * dy / 9.)
}

// pub struct SimpsonIntegration2D<F: Fn(f64, f64, usize) -> T, T> {
//   function: F,
// }

// impl<F, T> SimpsonIntegration2D<F, T>
// where
//   T: Zero + std::ops::Mul<f64, Output = T> + std::ops::Add<T, Output = T> + Send + Sync + Sum<T>,
//   F: Send + Sync + Fn(f64, f64, usize) -> T,
// {
//   /// Creates a new integrable function by using the supplied `Fn(f64, f64) -> T` in
//   /// combination with numeric integration via simpson's rule to find the integral.
//   pub fn new(function: F) -> Self {
//     SimpsonIntegration2D { function }
//   }

//   /// Get simpson weight for index
//   pub fn get_weight(nx: usize, ny: usize, divs: usize) -> f64 {
//     get_simpson_weight(nx, divs) * get_simpson_weight(ny, divs)
//   }

//   /// Numerically integrate from `a` to `b`, in `divs` divisions
//   pub fn integrate(&self, x_range: (f64, f64), y_range: (f64, f64), divs: usize) -> T {
//     assert!(divs.is_even());
//     assert!(divs >= 4);

//     let steps = divs + 1;
//     let dx = (x_range.1 - x_range.0) / (divs as f64);
//     let dy = (y_range.1 - y_range.0) / (divs as f64);

//     let y_vals: Vec<_> = Steps(y_range.0, y_range.1, steps)
//       .into_iter()
//       .enumerate()
//       .collect();
//     let x_vals: Vec<_> = Steps(x_range.0, x_range.1, steps)
//       .into_iter()
//       .enumerate()
//       .collect();

//     let result: T = y_vals
//       .par_iter()
//       .fold(
//         || T::zero(),
//         |acc, (ny, y)| {
//           let sy: T = x_vals
//             .par_iter()
//             .fold(
//               || T::zero(),
//               |acc, (nx, x)| {
//                 let a_n = get_simpson_weight(*nx, divs);
//                 acc + (self.function)(*x, *y, get_1d_index(*ny, *nx, steps)) * a_n
//               },
//             )
//             .sum();

//           let a_n = get_simpson_weight(*ny, divs);
//           acc + sy * a_n
//         },
//       )
//       .sum();

//     result * (dx * dy / 9.)
//   }
// }

//
// Adaptive simpson's rule
//
/// Evaluates the Simpson's Rule, also returning m and f(m) to reuse
fn quad_simpsons_mem<F, Y>(
  f: F,
  a: f64,
  fa: Complex<f64>,
  b: f64,
  fb: Complex<f64>,
) -> (f64, Complex<f64>, Complex<f64>)
where
  F: Fn(f64) -> Y,
  Y: Into<Complex<f64>>,
{
  let m = (a + b) / 2.0;
  let fm = f(m).into();
  let simpson_value = (b - a).abs() / 6.0 * (fa + 4.0 * fm + fb);
  (m, fm, simpson_value)
}

#[allow(clippy::too_many_arguments)]
/// Efficient recursive implementation of adaptive Simpson's rule
/// Function values at the start, middle, end of the intervals are retained
fn quad_asr<F, Y>(
  f: &F,
  a: f64,
  fa: Complex<f64>,
  b: f64,
  fb: Complex<f64>,
  eps: f64,
  whole: Complex<f64>,
  m: f64,
  fm: Complex<f64>,
  max_depth: usize,
) -> Complex<f64>
where
  F: Fn(f64) -> Y,
  Y: Into<Complex<f64>>,
{
  if max_depth == 0 || (eps / 2.0) == eps || (b - a).abs() < f64::EPSILON {
    return whole;
  }

  let (lm, flm, left) = quad_simpsons_mem(f, a, fa, m, fm);
  let (rm, frm, right) = quad_simpsons_mem(f, m, fm, b, fb);
  let delta = left + right - whole;

  if delta.norm() <= 15.0 * eps {
    left + right + delta / 15.0
  } else {
    quad_asr(f, a, fa, m, fm, eps / 2.0, left, lm, flm, max_depth - 1)
      + quad_asr(f, m, fm, b, fb, eps / 2.0, right, rm, frm, max_depth - 1)
  }
}

/// Integrate f from a to b using Adaptive Simpson's Rule with max error of eps
pub fn simpson_adaptive<F, Y>(f: &F, a: f64, b: f64, eps: f64, max_depth: usize) -> Complex<f64>
where
  F: Fn(f64) -> Y,
  Y: Into<Complex<f64>>,
{
  let fa = f(a).into();
  let fb = f(b).into();
  // change eps to be relative to function values
  // let eps = eps * (fa.norm() + fb.norm()) / 2.0;
  let (m, fm, whole) = quad_simpsons_mem(f, a, fa, b, fb);
  quad_asr(f, a, fa, b, fb, eps, whole, m, fm, max_depth)
}

fn simpson_adaptive_2d<F, Y>(
  f: &F,
  ax: f64,
  bx: f64,
  ay: f64,
  by: f64,
  eps: f64,
  max_depth: usize,
) -> Complex<f64>
where
  F: Fn(f64, f64) -> Y,
  Y: Into<Complex<f64>>,
{
  let g = |x| simpson_adaptive(&|y| f(x, y), ay, by, eps, max_depth);
  simpson_adaptive(&g, ax, bx, eps, max_depth)
}

#[cfg(test)]
mod tests {
  use super::*;
  extern crate float_cmp;
  use float_cmp::*;
  use std::f64::consts::PI;

  #[test]
  fn integrator_test() {
    let actual = simpson(|x| x.sin(), 0., PI, 1000).re;

    let expected = 2.;

    // println!("divisions {}", calc_required_divisions_for_simpson_precision(|x| 10e6 * x.sin(), 0., PI, 1e-8));

    assert!(
      approx_eq!(f64, actual, expected, ulps = 2, epsilon = 1e-11),
      "actual: {}, expected: {}",
      actual,
      expected
    );
  }

  #[test]
  fn integrator_2d_test() {
    let divs = 1000;
    let actual = simpson2d(|x, y| x.sin() * y.powi(3), 0., PI, 0., 2., divs).re;
    // let actual = SimpsonIntegration2D::new(|x, y, _| {
    //   let zd = f64::cos(PI * f64::cos(x) / 2.0) * f64::cos(PI * (1.0 - f64::sin(x) * f64::cos(y)) / 4.0);
    //   (zd.powi(2)/f64::sin(x)).abs()
    // }).integrate((0.1e-8, PI), (0.1e-8, 2. * PI), divs);

    let expected = (2_f64).powi(4) * 2. / 4.;

    assert!(
      approx_eq!(f64, actual, expected, ulps = 2, epsilon = 1e-11),
      "actual: {}, expected: {}",
      actual,
      expected
    );
  }

  #[test]
  fn adaptive_simpsons_rule_test() {
    let f = |x: f64| x.sin() * x.cos() * x.powi(2);
    let actual = simpson_adaptive(&f, 0., PI, 1e-8, 1000).re;
    let expected = simpson(f, 0., PI, 1000).re;

    assert!(
      approx_eq!(f64, actual, expected, ulps = 2, epsilon = 1e-8),
      "actual: {}, expected: {}",
      actual,
      expected
    );
  }
}
