use dim::ucum;
use crate::math::{lerp};

/// Create a dimensioned vector3
pub fn dim_vector3<L, R>(unit_const : L, arr : &[R; 3]) -> na::Vector3<dim::typenum::Prod<L, R>>
where
  L : std::ops::Mul<R> + Copy,
  R : Copy,
  dim::typenum::Prod<L, R> : na::Scalar,
{
  na::Vector3::new(
    unit_const * arr[0],
    unit_const * arr[1],
    unit_const * arr[2],
  )
}

/// convert from celsius to kelvin
pub fn from_celsius_to_kelvin(c : f64) -> ucum::Kelvin<f64> {
  ucum::Kelvin::new(c + 273.15)
}

/// convert from kelvin to celsius
pub fn from_kelvin_to_celsius(k : ucum::Kelvin<f64>) -> f64 {
  *(k / ucum::K) - 273.15
}

#[derive(Debug, Copy, Clone)]
pub struct Steps<T>(pub T, pub T, pub u32);

impl<T> IntoIterator for Steps<T>
where T: std::ops::Mul<f64, Output=T> + std::ops::Add<T, Output=T> + Copy {
  type Item = T;
  type IntoIter = StepIterator<T>;

  fn into_iter(self) -> Self::IntoIter {
    StepIterator {
      steps: self,
      index: 0,
    }
  }
}

pub struct StepIterator<T> {
  steps: Steps<T>,
  index: u32,
}

impl<T> Iterator for StepIterator<T>
where T: std::ops::Mul<f64, Output=T> + std::ops::Add<T, Output=T> + Copy {
  type Item = T;

  fn next(&mut self) -> Option<Self::Item> {
    if self.index >= self.steps.2 {
      return None;
    }

    let progress = (self.index as f64) / ((self.steps.2 - 1) as f64);
    self.index += 1;

    Some(lerp(self.steps.0, self.steps.1, progress))
  }
}

/// An iterator that will iterate through rows and columns, giving you the
/// coordinates at every iteration. Like a 2d linspace.
pub struct Iterator2D<T> {
  x_steps : Steps<T>,
  y_steps : Steps<T>,
  index : u32,
  total : u32,
}

impl<T> Iterator2D<T> {
  /// Create a new 2d iterator
  pub fn new(
    x_steps : Steps<T>,
    y_steps : Steps<T>
  ) -> Self {
    let total = x_steps.2 * y_steps.2;
    Iterator2D {
      x_steps,
      y_steps,
      total,
      index: 0,
    }
  }

  // get the 2d indices (row, column) from the linear index
  pub fn get_2d_indices( index : u32, shape : (u32, u32) ) -> (u32, u32) {
    (
      (index % shape.0),
      (index / shape.1)
    )
  }
}

impl<T> Iterator for Iterator2D<T>
where T: std::ops::Mul<f64, Output=T> + std::ops::Add<T, Output=T> + Copy {
  type Item = (T, T); // x, y

  fn next(&mut self) -> Option<Self::Item> {
    if self.index >= self.total {
      return None;
    }

    let shape = (self.x_steps.2, self.y_steps.2);
    let (nx, ny) = Self::get_2d_indices(self.index, shape);
    let xt = (nx as f64) / ((shape.0 - 1) as f64);
    let yt = (ny as f64) / ((shape.1 - 1) as f64);
    let x = lerp(self.x_steps.0, self.x_steps.1, xt);
    let y = lerp(self.y_steps.0, self.y_steps.1, yt);

    self.index += 1;

    Some((x, y))
  }
}
