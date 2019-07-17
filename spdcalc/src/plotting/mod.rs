use super::*;
use dim::ucum;
use math::lerp;
use na::Vector2;
use spd::*;

/// Holds configuration for drawing heatmaps
#[derive(Debug, Copy, Clone)]
pub struct HistogramConfig {
  /// the x axis range (min, max)
  pub x_range : (f64, f64),
  /// the y axis range (min, max)
  pub y_range : (f64, f64),
  /// the x axis number of bins
  pub x_count : usize,
  /// the y axis number of bins
  pub y_count : usize,
}

impl IntoIterator for HistogramConfig {
  type Item = Vector2<f64>; // x, y
  type IntoIter = HistogramConfigIterator;

  fn into_iter(self) -> Self::IntoIter {
    HistogramConfigIterator {
      cfg :   self,
      index : 0,
      total : self.x_count * self.y_count,
    }
  }
}

pub struct HistogramConfigIterator {
  cfg :   HistogramConfig,
  index : usize,
  total : usize,
}

impl Iterator for HistogramConfigIterator {
  type Item = Vector2<f64>; // x, y

  fn next(&mut self) -> Option<Self::Item> {
    if self.index >= self.total {
      return None;
    }

    let xt = ((self.index % self.cfg.x_count) as f64) / ((self.cfg.x_count - 1) as f64);
    let yt = ((self.index / self.cfg.y_count) as f64) / ((self.cfg.y_count - 1) as f64);
    let x = lerp(self.cfg.x_range.0, self.cfg.x_range.1, xt);
    let y = lerp(self.cfg.y_range.0, self.cfg.y_range.1, yt);

    self.index += 1;

    Some(Vector2::new(x, y))
  }
}

/// Create a JSI plot
pub fn plot_JSI(params : &SPD, cfg : &HistogramConfig) -> Vec<f64> {
  let norm_amp = phasematch_collinear(&params);
  // norm of intensity
  let norm = norm_amp.re.powi(2) + norm_amp.im.powi(2);

  cfg
    .into_iter()
    .map(|coords| {
      let l_s = coords.x;
      let l_i = coords.y;

      let mut signal = params.signal.clone();
      let mut idler = params.idler.clone();

      signal.set_wavelength(l_s * ucum::M);
      idler.set_wavelength(l_i * ucum::M);

      let spd = SPD {
        signal,
        idler,
        ..*params
      };

      let amplitude = phasematch(&spd);

      // intensity
      (amplitude.re.powi(2) + amplitude.im.powi(2)) / norm
    })
    .collect()
}