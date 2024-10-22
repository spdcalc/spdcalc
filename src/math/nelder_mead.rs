extern crate ndarray;
// extern crate optimize;
// use nelder_mead::{*, params::*, bounds::*};
// use argmin::prelude::*;
// use argmin::solver::neldermead::NelderMead;
// use serde::{Deserialize, Serialize};

use argmin::core::{CostFunction, Executor};
// use ndarray::{Array, ArrayView1};
use argmin::solver::neldermead::NelderMead;
// use optimize::*;

// #[derive(Clone, Default, Serialize, Deserialize)]
// struct NelderMead1d<F: FnMut(f64) -> f64> {
//   func: F,
// }
//
// impl<F : FnMut(f64) -> f64> ArgminOp for NelderMead1d<F>
// where
//   F: serde::Serialize + std::clone::Clone + std::marker::Sync +
// std::marker::Send {
//   /// Type of the parameter vector
//   type Param = Vec<f64>;
//   /// Type of the return value computed by the cost function
//   type Output = f64;
//   /// Type of the Hessian. Can be `()` if not needed.
//   type Hessian = ();
//   /// Type of the Jacobian. Can be `()` if not needed.
//   type Jacobian = ();
//
//   /// Apply the cost function to a parameter `p`
//   fn apply(&self, p: &Self::Param) -> Result<Self::Output, Error> {
//     Ok((self.func)(p[0]))
//   }
// }

struct Cost1d<F>
where
  F: Fn(f64) -> f64,
{
  func: F,
  min: f64,
  max: f64,
}
impl<F> CostFunction for Cost1d<F>
where
  F: Fn(f64) -> f64,
{
  type Param = f64;
  type Output = f64;
  fn cost(&self, x: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
    if x > &self.max || x < &self.min {
      Ok(f64::INFINITY)
    } else {
      let f = &self.func;
      Ok(f(*x))
    }
  }
}

/// nelder mead optimization. Returns x
pub fn nelder_mead_1d<F>(
  func: F,
  guess: (f64, f64),
  max_iter: u64,
  min: f64,
  max: f64,
  tolerance: f64,
) -> f64
where
  F: Fn(f64) -> f64,
{
  let cost = Cost1d { func, min, max };
  let solver = NelderMead::new(vec![guess.0, guess.1])
    .with_sd_tolerance(tolerance)
    .unwrap();

  let res = Executor::new(cost, solver)
    .configure(|state| state.max_iters(max_iter))
    .run()
    .unwrap();

  use argmin::core::State;
  *res.state().get_best_param().unwrap()

  // let minimizer = NelderMeadBuilder::default()
  //   .xtol(tolerance)
  //   .ftol(tolerance)
  //   .maxiter(max_iter as usize)
  //   .build()
  //   .unwrap();

  // let cost = |args : ArrayView1<f64>| {
  //   let x = args[0];
  //   if x > max || x < min {
  //     std::f64::INFINITY
  //   }
  //   // high cost if x outside bounds
  //   else {
  //     func(x)
  //   }
  // };

  // let ans = minimizer.minimize(&cost, Array::from_vec(vec![guess]).view());

  // ans[0]

  // let (x, _fx) = minimize(
  //   |args| func(args[0]),
  //   vec![guess],
  //   (max-min)/100.,
  //   Params::default(),
  //   Bounds { min: vec![min], max: vec![max] },
  //   max_iter
  // );
  //
  // x[0]

  // let cost = NelderMead1d { func };
  //
  // // Set up solver -- note that the proper choice of the vertices is very
  // important! let solver = NelderMead::new()
  //   .with_initial_params(vec![guess])
  //   .sd_tolerance(0.0001);
  //
  //   // Run solver
  // let res = Executor::new(cost, solver, vec![])
  //   .max_iters(100)
  //   .run()?;
  //
  // res
}
