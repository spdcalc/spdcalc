//! # Linear Interpolation Utility
//!
//! Provides a generic 1D linear interpolator that works with any output dimension.
//!
//! ## Example
//! ```
//! use spdcalc::math::Interpolator;
//!
//! // Create an interpolator with 2D output (e.g., for refractive indices)
//! let inputs = vec![400.0, 500.0, 600.0];
//! let outputs = vec![[1.66, 1.55], [1.65, 1.54], [1.64, 1.53]];
//! let interp = Interpolator::<2>::new(inputs, outputs).unwrap();
//!
//! // Interpolate at a specific input value
//! let result = interp.interpolate(450.0);
//! // result is approximately [1.655, 1.545] (midpoint between first two points)
//! ```

use super::lerp;
use crate::SPDCError;

/// A generic 1D linear interpolator
///
/// Supports any fixed-size output dimension via const generics.
/// Out-of-range inputs are clamped to boundary values (constant extrapolation).
#[derive(Debug, Clone, PartialEq)]
pub struct Interpolator<const N: usize> {
  inputs: Vec<f64>,
  outputs: Vec<[f64; N]>,
}

impl<const N: usize> Interpolator<N> {
  /// Create a new interpolator from input and output data
  ///
  /// ## Validation
  /// * Input and output arrays must have the same length
  /// * At least one data point is required
  /// * Input values must be unique and monotonically increasing
  ///
  /// ## Example
  /// ```
  /// use spdcalc::math::Interpolator;
  ///
  /// let inputs = vec![1.0, 2.0, 3.0];
  /// let outputs = vec![[10.0], [20.0], [30.0]];
  /// let interp = Interpolator::<1>::new(inputs, outputs).unwrap();
  /// ```
  pub fn new(inputs: Vec<f64>, outputs: Vec<[f64; N]>) -> Result<Self, SPDCError> {
    // Validate lengths match
    if inputs.len() != outputs.len() {
      return Err(SPDCError::new(format!(
        "Input and output arrays must have the same length: got {} inputs and {} outputs",
        inputs.len(),
        outputs.len()
      )));
    }

    // Validate non-empty
    if inputs.is_empty() {
      return Err(SPDCError::new("At least one data point is required"));
    }

    // Check monotonicity and uniqueness
    let is_monotonic = inputs.iter().skip(1).zip(inputs.iter()).all(|(a, b)| a > b);

    if !is_monotonic {
      return Err(SPDCError::new("Input values must be unique and monotonically increasing"));
    }

    Ok(Self {
      inputs,
      outputs,
    })
  }

  /// Interpolate output value(s) at a given input
  ///
  /// Uses linear interpolation between data points.
  /// Out-of-range inputs are clamped to boundary values.
  ///
  /// ## Arguments
  /// * `input` - Input value to interpolate at
  ///
  /// ## Returns
  /// Interpolated output value(s)
  ///
  /// ## Example
  /// ```
  /// use spdcalc::math::Interpolator;
  ///
  /// let inputs = vec![0.0, 10.0];
  /// let outputs = vec![[0.0, 100.0], [10.0, 200.0]];
  /// let interp = Interpolator::<2>::new(inputs, outputs).unwrap();
  ///
  /// let result = interp.interpolate(5.0);
  /// // result is [5.0, 150.0] (midpoint)
  /// ```
  pub fn interpolate(&self, input: f64) -> [f64; N] {
    // Handle single point case
    if self.inputs.len() == 1 {
      return self.outputs[0];
    }

    // Binary search for insertion point
    match self.inputs.binary_search_by(|probe| {
      probe
        .partial_cmp(&input)
        .unwrap_or(std::cmp::Ordering::Equal)
    }) {
      // Exact match
      Ok(idx) => self.outputs[idx],

      // Between points or out of range
      Err(insert_idx) => {
        // Below range - use first value (constant extrapolation)
        if insert_idx == 0 {
          self.outputs[0]
        }
        // Above range - use last value (constant extrapolation)
        else if insert_idx >= self.inputs.len() {
          self.outputs[self.inputs.len() - 1]
        }
        // Interpolate between insert_idx-1 and insert_idx
        else {
          let i0 = insert_idx - 1;
          let i1 = insert_idx;

          let input0 = self.inputs[i0];
          let input1 = self.inputs[i1];

          // Linear interpolation parameter
          let t = (input - input0) / (input1 - input0);

          // Interpolate each output dimension
          std::array::from_fn(|i| lerp(self.outputs[i0][i], self.outputs[i1][i], t))
        }
      }
    }
  }

  /// Get the input range covered by this interpolator
  ///
  /// Returns `(min_input, max_input)`
  pub fn input_bounds(&self) -> (f64, f64) {
    (self.inputs[0], self.inputs[self.inputs.len() - 1])
  }

  /// Get the number of data points
  pub fn len(&self) -> usize {
    self.inputs.len()
  }

  /// Check if the interpolator is empty
  pub fn is_empty(&self) -> bool {
    self.inputs.is_empty()
  }
}

impl Interpolator<1> {
  /// Convenience method to create a 1D output interpolator from f64 outputs
  pub fn new_1d(inputs: Vec<f64>, outputs: Vec<f64>) -> Result<Self, SPDCError> {
    let outputs_2d: Vec<[f64; 1]> = outputs.into_iter().map(|v| [v]).collect();
    Self::new(inputs, outputs_2d)
  }

  /// Convenience method to get a single f64 output instead of [f64; 1]
  pub fn interpolate_scalar(&self, input: f64) -> f64 {
    self.interpolate(input)[0]
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use float_cmp::approx_eq;

  #[test]
  fn test_new_valid_data() {
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![10.0, 20.0, 30.0];
    let result = Interpolator::new_1d(inputs, outputs);
    assert!(result.is_ok());
  }

  #[test]
  fn test_new_empty_data() {
    let inputs: Vec<f64> = vec![];
    let outputs: Vec<[f64; 1]> = vec![];
    let result = Interpolator::<1>::new(inputs, outputs);
    assert!(result.is_err());
    assert!(result.unwrap_err().0.contains("At least one data point"));
  }

  #[test]
  fn test_new_length_mismatch() {
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![[10.0], [20.0]]; // Only 2 outputs
    let result = Interpolator::<1>::new(inputs, outputs);
    assert!(result.is_err());
    assert!(result.unwrap_err().0.contains("same length"));
  }

  #[test]
  fn test_new_duplicate_inputs() {
    let inputs = vec![1.0, 2.0, 2.0, 3.0]; // Duplicate 2.0
    let outputs = vec![[10.0], [20.0], [20.0], [30.0]];
    let result = Interpolator::<1>::new(inputs, outputs);
    assert!(result.is_err());
  }

  #[test]
  fn test_new_unsorted() {
    let inputs = vec![3.0, 1.0, 2.0]; // Unsorted
    let outputs = vec![[30.0], [10.0], [20.0]];
    let result = Interpolator::<1>::new(inputs, outputs);
    // should fail
    assert!(result.is_err());
  }

  #[test]
  fn test_interpolate_exact_match() {
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![[10.0], [20.0], [30.0]];
    let interp = Interpolator::<1>::new(inputs, outputs).unwrap();

    let result = interp.interpolate(2.0);
    assert_eq!(result, [20.0]);
  }

  #[test]
  fn test_interpolate_midpoint() {
    let inputs = vec![0.0, 10.0];
    let outputs = vec![[0.0], [10.0]];
    let interp = Interpolator::<1>::new(inputs, outputs).unwrap();

    let result = interp.interpolate(5.0);
    assert!(approx_eq!(f64, result[0], 5.0, epsilon = 1e-10));
  }

  #[test]
  fn test_interpolate_quarter_point() {
    let inputs = vec![0.0, 100.0];
    let outputs = vec![[0.0], [100.0]];
    let interp = Interpolator::<1>::new(inputs, outputs).unwrap();

    let result = interp.interpolate(25.0);
    assert!(approx_eq!(f64, result[0], 25.0, epsilon = 1e-10));
  }

  #[test]
  fn test_interpolate_below_range() {
    let inputs = vec![10.0, 20.0];
    let outputs = vec![[100.0], [200.0]];
    let interp = Interpolator::<1>::new(inputs, outputs).unwrap();

    let result = interp.interpolate(5.0); // Below range
    assert_eq!(result, [100.0]); // Should clamp to first value
  }

  #[test]
  fn test_interpolate_above_range() {
    let inputs = vec![10.0, 20.0];
    let outputs = vec![[100.0], [200.0]];
    let interp = Interpolator::<1>::new(inputs, outputs).unwrap();

    let result = interp.interpolate(25.0); // Above range
    assert_eq!(result, [200.0]); // Should clamp to last value
  }

  #[test]
  fn test_interpolate_single_point() {
    let inputs = vec![5.0];
    let outputs = vec![42.0];
    let interp = Interpolator::new_1d(inputs, outputs).unwrap();

    assert_eq!(interp.interpolate_scalar(0.0), 42.0);
    assert_eq!(interp.interpolate_scalar(5.0), 42.0);
    assert_eq!(interp.interpolate_scalar(10.0), 42.0);
  }

  // Test with 2D output (our use case for refractive indices)
  #[test]
  fn test_interpolate_2d_output() {
    let inputs = vec![400.0, 500.0, 600.0];
    let outputs = vec![[1.66, 1.55], [1.65, 1.54], [1.64, 1.53]];
    let interp = Interpolator::<2>::new(inputs, outputs).unwrap();

    // Test exact match
    let result = interp.interpolate(500.0);
    assert_eq!(result, [1.65, 1.54]);

    // Test interpolation (midpoint between 400 and 500)
    let result = interp.interpolate(450.0);
    assert!(approx_eq!(f64, result[0], 1.655, epsilon = 1e-10));
    assert!(approx_eq!(f64, result[1], 1.545, epsilon = 1e-10));
  }

  // Test with 3D output (verify generics work)
  #[test]
  fn test_interpolate_3d_output() {
    let inputs = vec![0.0, 10.0];
    let outputs = vec![[0.0, 0.0, 0.0], [10.0, 20.0, 30.0]];
    let interp = Interpolator::<3>::new(inputs, outputs).unwrap();

    let result = interp.interpolate(5.0); // Midpoint
    assert!(approx_eq!(f64, result[0], 5.0, epsilon = 1e-10));
    assert!(approx_eq!(f64, result[1], 10.0, epsilon = 1e-10));
    assert!(approx_eq!(f64, result[2], 15.0, epsilon = 1e-10));
  }

  #[test]
  fn test_input_bounds() {
    let inputs = vec![100.0, 200.0, 300.0];
    let outputs = vec![[1.0], [2.0], [3.0]];
    let interp = Interpolator::<1>::new(inputs, outputs).unwrap();

    let (min, max) = interp.input_bounds();
    assert_eq!(min, 100.0);
    assert_eq!(max, 300.0);
  }

  #[test]
  fn test_len() {
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![[10.0], [20.0], [30.0]];
    let interp = Interpolator::<1>::new(inputs, outputs).unwrap();

    assert_eq!(interp.len(), 3);
    assert!(!interp.is_empty());
  }
}
