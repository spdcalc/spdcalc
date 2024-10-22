use dim::ucum::Meter2;

use crate::Wavelength;

/// Beam waist
///
/// Can be circular or elliptic.
///
/// ## Example
/// ```
/// use spdcalc::prelude::*;
/// use spdcalc::beam::BeamWaist;
/// let w = BeamWaist::new(100.0 * MICRO * M);
/// assert_eq!(w.x, 100.0 * MICRO * M);
/// assert_eq!(w.y, 100.0 * MICRO * M);
/// assert_eq!(w, (100.0 * MICRO * M).into());
/// ```
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct BeamWaist {
  /// wx at 1/e^2
  pub x: Wavelength,
  /// wy at 1/e^2
  pub y: Wavelength,
}

impl BeamWaist {
  /// Create a new circular beam waist
  pub fn new(wx: Wavelength) -> Self {
    Self { x: wx, y: wx }
  }

  /// Create a new elliptic beam waist
  ///
  /// Note elliptic beams are not supported in the calculations yet.
  #[cfg(feature = "elliptic")]
  pub fn new_elliptic(wx: Wavelength, wy: Wavelength) -> Self {
    Self { x: wx, y: wy }
  }

  #[cfg(feature = "elliptic")]
  /// Create a new elliptic beam waist
  pub fn ellipticity(&self) -> f64 {
    if self.x == self.y {
      1.
    } else if self.y < self.x {
      *(self.y / self.x)
    } else {
      *(self.x / self.y)
    }
  }

  /// Calculate the beam waist x by y
  pub fn x_by_y(&self) -> Meter2<f64> {
    self.x * self.y
  }
}

impl From<Wavelength> for BeamWaist {
  fn from(wx: Wavelength) -> Self {
    Self::new(wx)
  }
}

#[cfg(feature = "elliptic")]
impl From<(Wavelength, Wavelength)> for BeamWaist {
  fn from(wxwy: (Wavelength, Wavelength)) -> Self {
    Self::new_elliptic(wxwy.0, wxwy.1)
  }
}
