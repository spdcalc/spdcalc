//! Various iterators over the signal and idler in frequency/wavelength space
//!
use crate::{
  utils::{frequency_to_vacuum_wavelength, vacuum_wavelength_to_frequency, Steps2D},
  Frequency, Wavelength,
};
use rayon::prelude::*;

/// A range of signal and idler frequencies
#[derive(Debug, Clone, Copy)]
pub struct FrequencySpace(Steps2D<Frequency>);

impl FrequencySpace {
  /// Create a new frequency space
  pub fn new(xsteps: (Frequency, Frequency, usize), ysteps: (Frequency, Frequency, usize)) -> Self {
    Self(Steps2D(xsteps, ysteps))
  }

  /// Get the steps
  pub fn steps(&self) -> &Steps2D<Frequency> {
    &self.0
  }

  /// Convert to steps
  pub fn as_steps(self) -> Steps2D<Frequency> {
    self.0
  }

  /// Set the resolution
  pub fn set_resolution(&mut self, res: usize) -> &mut Self {
    self.0 .0 .2 = res;
    self.0 .1 .2 = res;
    self
  }

  /// Set the resolution by consuming self
  pub fn with_resolution(mut self, res: usize) -> Self {
    self.set_resolution(res);
    self
  }

  /// Iterate over the signal and idler frequencies
  pub fn from_wavelength_space(ws: WavelengthSpace) -> Self {
    let ws_min = vacuum_wavelength_to_frequency(ws.0 .0 .1);
    let ws_max = vacuum_wavelength_to_frequency(ws.0 .0 .0);
    let wi_min = vacuum_wavelength_to_frequency(ws.0 .1 .1);
    let wi_max = vacuum_wavelength_to_frequency(ws.0 .1 .0);
    Self::new((ws_min, ws_max, ws.0 .0 .2), (wi_min, wi_max, ws.0 .1 .2))
  }

  /// Convert to wavelength space
  pub fn as_wavelength_space(self) -> WavelengthSpace {
    let fs = self.as_steps();
    let ls_min = frequency_to_vacuum_wavelength(fs.0 .1);
    let ls_max = frequency_to_vacuum_wavelength(fs.0 .0);
    let li_min = frequency_to_vacuum_wavelength(fs.1 .1);
    let li_max = frequency_to_vacuum_wavelength(fs.1 .0);
    WavelengthSpace::new((ls_min, ls_max, fs.0 .2), (li_min, li_max, fs.1 .2))
  }

  /// Create a new frequency space from a sum-difference frequency space
  pub fn from_sum_diff_space(sdfs: SumDiffFrequencySpace) -> Self {
    sdfs.as_frequency_space()
  }

  /// Convert to sum-difference frequency space
  pub fn as_sum_diff_space(self) -> SumDiffFrequencySpace {
    SumDiffFrequencySpace::from_frequency_space(self)
  }
}

impl From<Steps2D<Frequency>> for FrequencySpace {
  fn from(steps: Steps2D<Frequency>) -> Self {
    Self(steps)
  }
}

impl From<WavelengthSpace> for FrequencySpace {
  fn from(ws: WavelengthSpace) -> Self {
    FrequencySpace::from_wavelength_space(ws)
  }
}

impl From<SumDiffFrequencySpace> for FrequencySpace {
  fn from(sd: SumDiffFrequencySpace) -> Self {
    sd.as_frequency_space()
  }
}

/// Provides an iterator over the signal and idler frequencies
pub trait IntoSignalIdlerIterator {
  /// Get an iterator over the signal and idler frequencies
  fn into_signal_idler_iterator(self) -> impl Iterator<Item = (Frequency, Frequency)>;
  /// Get a parallel iterator over the signal and idler frequencies
  fn into_signal_idler_par_iterator(self) -> impl ParallelIterator<Item = (Frequency, Frequency)>;
}

impl IntoSignalIdlerIterator for FrequencySpace {
  fn into_signal_idler_iterator(self) -> impl Iterator<Item = (Frequency, Frequency)> {
    self.0.into_iter()
  }

  fn into_signal_idler_par_iterator(self) -> impl ParallelIterator<Item = (Frequency, Frequency)> {
    self.0.into_par_iter()
  }
}

/// A 45 degree rotation in frequency space.
///
/// X-axis is half the sum of the signal and idler frequencies,
/// and the Y-axis is half the difference.
#[derive(Debug, Clone, Copy)]
pub struct SumDiffFrequencySpace(Steps2D<Frequency>);

impl SumDiffFrequencySpace {
  /// Create a new sum-difference frequency space
  pub fn new(xsteps: (Frequency, Frequency, usize), ysteps: (Frequency, Frequency, usize)) -> Self {
    Self(Steps2D(xsteps, ysteps))
  }

  /// Get the steps
  pub fn steps(&self) -> &Steps2D<Frequency> {
    &self.0
  }

  /// Convert to steps
  pub fn as_steps(self) -> Steps2D<Frequency> {
    self.0
  }

  /// Set the resolution
  pub fn set_resolution(&mut self, res: usize) -> &mut Self {
    self.0 .0 .2 = res;
    self.0 .1 .2 = res;
    self
  }

  /// Set the resolution by consuming self
  pub fn with_resolution(mut self, res: usize) -> Self {
    self.set_resolution(res);
    self
  }

  /// Create a new sum-difference frequency space from a frequency space
  pub fn from_frequency_space(frequencies: FrequencySpace) -> Self {
    let steps = frequencies.as_steps();
    let ws_min = steps.0 .0;
    let ws_max = steps.0 .1;
    let wi_min = steps.1 .0;
    let wi_max = steps.1 .1;
    //x: s = (wi + ws) / 2
    //y: d = (wi - ws) / 2
    let s_min = (wi_min + ws_min) / 2.;
    let s_max = (wi_max + ws_max) / 2.;
    let d_min = (wi_min - ws_max) / 2.;
    let d_max = (wi_max - ws_min) / 2.;

    Self(Steps2D(
      (s_min, s_max, steps.0 .2),
      (d_min, d_max, steps.1 .2),
    ))
  }

  /// Convert to frequency space
  pub fn as_frequency_space(self) -> FrequencySpace {
    let s_min = self.0 .0 .0;
    let s_max = self.0 .0 .1;
    let d_min = self.0 .1 .0;
    let d_max = self.0 .1 .1;
    let ws_min = 0.25 * (3. * s_min + s_max - d_min - 3. * d_max);
    let ws_max = 0.25 * (s_min + 3. * s_max - 3. * d_min - d_max);
    let wi_min = 0.25 * (3. * s_min + s_max + 3. * d_min + d_max);
    let wi_max = 0.25 * (s_min + 3. * s_max + d_min + 3. * d_max);
    FrequencySpace::new(
      (ws_min, ws_max, self.0 .0 .2),
      (wi_min, wi_max, self.0 .1 .2),
    )
  }

  /// Create a new sum-difference frequency space from a wavelength space
  pub fn from_wavelength_space(ws: WavelengthSpace) -> Self {
    Self::from_frequency_space(ws.as_frequency_space())
  }

  /// Convert to wavelength space
  pub fn as_wavelength_space(self) -> WavelengthSpace {
    WavelengthSpace::from_frequency_space(self.as_frequency_space())
  }
}

impl From<Steps2D<Frequency>> for SumDiffFrequencySpace {
  fn from(steps: Steps2D<Frequency>) -> Self {
    Self(steps)
  }
}

impl From<WavelengthSpace> for SumDiffFrequencySpace {
  fn from(ws: WavelengthSpace) -> Self {
    Self::from_wavelength_space(ws)
  }
}

impl From<FrequencySpace> for SumDiffFrequencySpace {
  fn from(fs: FrequencySpace) -> Self {
    Self::from_frequency_space(fs)
  }
}

impl IntoSignalIdlerIterator for SumDiffFrequencySpace {
  fn into_signal_idler_iterator(self) -> impl Iterator<Item = (Frequency, Frequency)> {
    self.0.into_iter().map(|(s, d)| {
      let w_s = s - d;
      let w_i = s + d;
      (w_s, w_i)
    })
  }

  fn into_signal_idler_par_iterator(self) -> impl ParallelIterator<Item = (Frequency, Frequency)> {
    self.0.into_par_iter().map(|(s, d)| {
      let w_s = s - d;
      let w_i = s + d;
      (w_s, w_i)
    })
  }
}

/// A range of signal and idler wavelengths
#[derive(Debug, Clone, Copy)]
pub struct WavelengthSpace(Steps2D<Wavelength>);

impl WavelengthSpace {
  /// Create a new wavelength space
  pub fn new(
    xsteps: (Wavelength, Wavelength, usize),
    ysteps: (Wavelength, Wavelength, usize),
  ) -> Self {
    Self(Steps2D(xsteps, ysteps))
  }

  /// Get the steps
  pub fn steps(&self) -> &Steps2D<Wavelength> {
    &self.0
  }

  /// Convert to steps
  pub fn as_steps(self) -> Steps2D<Wavelength> {
    self.0
  }

  /// Set the resolution
  pub fn set_resolution(&mut self, res: usize) -> &mut Self {
    self.0 .0 .2 = res;
    self.0 .1 .2 = res;
    self
  }

  /// Set the resolution by consuming self
  pub fn with_resolution(mut self, res: usize) -> Self {
    self.set_resolution(res);
    self
  }

  /// Create a new WavelengthSpace from a FrequencySpace
  pub fn from_frequency_space(fs: FrequencySpace) -> Self {
    fs.as_wavelength_space()
  }

  /// Convert to FrequencySpace
  pub fn as_frequency_space(self) -> FrequencySpace {
    FrequencySpace::from_wavelength_space(self)
  }

  /// Create a new WavelengthSpace from a SumDiffFrequencySpace
  pub fn from_sum_diff_space(sd: SumDiffFrequencySpace) -> Self {
    Self::from_frequency_space(sd.as_frequency_space())
  }

  /// Convert to SumDiffFrequencySpace
  pub fn as_sum_diff_space(self) -> SumDiffFrequencySpace {
    SumDiffFrequencySpace::from_wavelength_space(self)
  }
}

impl From<Steps2D<Wavelength>> for WavelengthSpace {
  fn from(steps: Steps2D<Wavelength>) -> Self {
    Self(steps)
  }
}

impl From<FrequencySpace> for WavelengthSpace {
  fn from(fs: FrequencySpace) -> Self {
    Self::from_frequency_space(fs)
  }
}

impl From<SumDiffFrequencySpace> for WavelengthSpace {
  fn from(sd: SumDiffFrequencySpace) -> Self {
    Self::from_sum_diff_space(sd)
  }
}

impl IntoSignalIdlerIterator for WavelengthSpace {
  fn into_signal_idler_iterator(self) -> impl Iterator<Item = (Frequency, Frequency)> {
    self.0.into_iter().map(|(ls, li)| {
      (
        vacuum_wavelength_to_frequency(ls),
        vacuum_wavelength_to_frequency(li),
      )
    })
  }

  fn into_signal_idler_par_iterator(self) -> impl ParallelIterator<Item = (Frequency, Frequency)> {
    self.0.into_par_iter().map(|(ls, li)| {
      (
        vacuum_wavelength_to_frequency(ls),
        vacuum_wavelength_to_frequency(li),
      )
    })
  }
}

/// A flat array holding signal and idler wavelengths
#[derive(Debug, Clone)]
pub struct SignalIdlerWavelengthArray(pub Vec<Wavelength>);

impl IntoSignalIdlerIterator for SignalIdlerWavelengthArray {
  fn into_signal_idler_iterator(self) -> impl Iterator<Item = (Frequency, Frequency)> {
    let chunked = self
      .0
      .chunks_exact(2)
      .map(|a| (a[0], a[1]))
      .collect::<Vec<_>>();
    chunked.into_iter().map(|(ls, li)| {
      (
        vacuum_wavelength_to_frequency(ls),
        vacuum_wavelength_to_frequency(li),
      )
    })
  }

  fn into_signal_idler_par_iterator(self) -> impl ParallelIterator<Item = (Frequency, Frequency)> {
    let chunked = self
      .0
      .chunks_exact(2)
      .map(|a| (a[0], a[1]))
      .collect::<Vec<_>>();
    chunked.into_par_iter().map(|(ls, li)| {
      (
        vacuum_wavelength_to_frequency(ls),
        vacuum_wavelength_to_frequency(li),
      )
    })
  }
}

/// A flat array holding signal and idler frequencies
#[derive(Debug, Clone)]
pub struct SignalIdlerFrequencyArray(pub Vec<Frequency>);

impl IntoSignalIdlerIterator for SignalIdlerFrequencyArray {
  fn into_signal_idler_iterator(self) -> impl Iterator<Item = (Frequency, Frequency)> {
    let chunked = self
      .0
      .chunks_exact(2)
      .map(|a| (a[0], a[1]))
      .collect::<Vec<_>>();
    chunked.into_iter()
  }

  fn into_signal_idler_par_iterator(self) -> impl ParallelIterator<Item = (Frequency, Frequency)> {
    let chunked = self
      .0
      .chunks_exact(2)
      .map(|a| (a[0], a[1]))
      .collect::<Vec<_>>();
    chunked.into_par_iter()
  }
}

#[cfg(test)]
mod test {
  use super::*;
  use crate::dim::ucum::*;
  use crate::math::Integrator;
  use crate::spdc::SPDC;

  #[test]
  fn test_si_arrays() {
    let spdc = SPDC::default();
    let spectrum = spdc.joint_spectrum(Integrator::default());
    let range = WavelengthSpace::new(
      (1400e-9 * M, 1600e-9 * M, 10),
      (1400e-9 * M, 1600e-9 * M, 10),
    );

    let jsi = spectrum.jsi_range(range);

    let values: Vec<Wavelength> = range
      .as_steps()
      .into_iter()
      .flat_map(|(s, i)| [s, i])
      .collect();
    let jsi2 = spectrum.jsi_range(SignalIdlerWavelengthArray(values));

    assert_eq!(jsi, jsi2);
  }
}
