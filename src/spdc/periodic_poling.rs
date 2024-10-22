use crate::dim::ucum::{Meter, M, RAD};
use crate::math::lerp;
use crate::{
  beam::IdlerBeam,
  beam::PumpBeam,
  beam::SignalBeam,
  delta_k,
  math::{fwhm_to_sigma, nelder_mead_1d},
  CrystalSetup, SPDCError, Sign, TWO_PI,
};
use crate::{Distance, Wavenumber};
use core::f64::consts::PI;

const IMPOSSIBLE_POLING_PERIOD: &str = "Could not determine poling period from specified values";

/// Units for the poling period (meters)
pub type PolingPeriod = Meter<f64>;

/// Apodization for periodic poling
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
#[serde(
  try_from = "crate::ApodizationConfig",
  into = "crate::ApodizationConfig"
)]
pub enum Apodization {
  /// None
  #[default]
  Off,
  /// Gaussian
  Gaussian {
    /// Full-width half-max
    fwhm: Distance,
  },
  /// Bartlett apodization
  Bartlett(f64),
  /// Blackman apodization
  Blackman(f64),
  /// Connes apodization
  Connes(f64),
  /// Cosine apodization
  Cosine(f64),
  /// Hamming apodization
  Hamming(f64),
  /// Welch apodization
  Welch(f64),
  /// Custom apodization by specifying profile values directly
  Interpolate(Vec<f64>),
}

impl Apodization {
  /// Get the kind of apodization as a string
  pub fn kind(&self) -> &'static str {
    match self {
      Apodization::Off => "Off",
      Apodization::Gaussian { .. } => "Gaussian",
      Apodization::Bartlett(_) => "Bartlett",
      Apodization::Blackman(_) => "Blackman",
      Apodization::Connes(_) => "Connes",
      Apodization::Cosine(_) => "Cosine",
      Apodization::Hamming(_) => "Hamming",
      Apodization::Welch(_) => "Welch",
      Apodization::Interpolate(_) => "Interpolate",
    }
  }

  /// The integration constant for the apodization function at a given z position
  pub fn integration_constant(&self, z: f64, crystal_length: Distance) -> f64 {
    assert!((-1. ..=1.).contains(&z), "z must be between -1 and 1");
    match self {
      // https://mathworld.wolfram.com/ApodizationFunction.html
      Apodization::Off => 1.,
      &Apodization::Gaussian { fwhm } => {
        let bw = 2. * fwhm_to_sigma(fwhm) / crystal_length;
        f64::exp(-0.5 * (z / bw).powi(2))
      }
      Apodization::Bartlett(a) => 1. - z.abs() / a,
      Apodization::Blackman(a) => {
        21. / 50. + 0.5 * (PI * z / a).cos() + (2. / 25.) * (TWO_PI * z / a).cos()
      }
      Apodization::Connes(a) => (1. - (z / a).powi(2)).powi(2),
      Apodization::Cosine(a) => (0.5 * PI * z / a).cos(),
      Apodization::Hamming(a) => (27. + 23. * (PI * z / a).cos()) / 50.,
      Apodization::Welch(a) => 1. - (z / a).powi(2),
      Apodization::Interpolate(values) => {
        // todo!("Interpolation not implemented yet")
        let n = values.len();
        if n == 0 {
          return 1.;
        }
        let i = 0.5 * (z + 1.) * (n - 1) as f64;
        let above = i.ceil() as usize;
        let below = i.floor() as usize;
        let t = i - below as f64;
        lerp(values[below], values[above], t)
      }
    }
  }
}

/// Periodic Poling settings
#[derive(Debug, Clone, PartialEq, Default)]
pub enum PeriodicPoling {
  /// Off
  #[default]
  Off,
  /// On
  On {
    /// Period of the poling
    period: PolingPeriod,
    /// Sign of the poling
    sign: Sign,
    /// Apodization
    apodization: Apodization,
  },
}

impl AsRef<PeriodicPoling> for PeriodicPoling {
  fn as_ref(&self) -> &Self {
    self
  }
}

impl PeriodicPoling {
  /// Create a new instance of periodic poling. Sign of specified period will be used.
  pub fn new(period: PolingPeriod, apodization: Apodization) -> Self {
    Self::On {
      period: if period > 0. * M { period } else { -period },
      sign: if period > 0. * M {
        Sign::POSITIVE
      } else {
        Sign::NEGATIVE
      },
      apodization,
    }
  }

  /// Set the poling period and sign
  pub fn with_period(self, period: PolingPeriod) -> Self {
    match self {
      Self::Off => Self::new(period, Apodization::Off),
      Self::On { apodization, .. } => Self::new(period, apodization),
    }
  }

  /// Set the period
  pub fn assign_period(&mut self, period: PolingPeriod) -> &mut Self {
    if let Self::On {
      period: p, sign, ..
    } = self
    {
      use dim::Abs;
      *p = period.abs();
      *sign = if period > 0. * M {
        Sign::POSITIVE
      } else {
        Sign::NEGATIVE
      };
    };
    self
  }

  /// Get the number of domains for a given crystal length
  pub fn num_domains(&self, crystal_length: Distance) -> usize {
    if let Self::On { period, .. } = self {
      let num_domains = (crystal_length / *period).ceil();
      num_domains as usize
    } else {
      0
    }
  }

  /// Get the poling domains
  ///
  /// They are a list of fractions of poling period
  pub fn poling_domains(&self, crystal_length: Distance) -> Vec<(f64, f64)> {
    if let Self::On { apodization, .. } = self {
      let num_domains = self.num_domains(crystal_length);
      // for each domain
      (0..num_domains)
        .map(|i| {
          // 0.5 to get value in middle of domain
          let z = lerp(-1., 1., (i as f64 + 0.5) / num_domains as f64);
          let a = apodization.integration_constant(z, crystal_length);
          let x = (1. - 2. * a.powi(2)).acos() / TWO_PI;
          if z > 0. {
            (1. - x, x)
          } else {
            (x, 1. - x)
          }
        })
        .collect()
    } else {
      Vec::new()
    }
  }

  /// Get the poling domains as lengths
  pub fn poling_domain_lengths(&self, crystal_length: Distance) -> Vec<(Distance, Distance)> {
    let period = match self {
      Self::Off => 0. * M,
      Self::On { period, .. } => *period,
    };
    self
      .poling_domains(crystal_length)
      .iter()
      .map(|&(z1, z2)| (z1 * period, z2 * period))
      .collect()
  }

  /// Get the optimal periodic poling for the given signal and pump beams
  pub fn try_new_optimum(
    signal: &SignalBeam,
    pump: &PumpBeam,
    crystal_setup: &CrystalSetup,
    apodization: Apodization,
  ) -> Result<Self, SPDCError> {
    let period = optimum_poling_period(signal, pump, crystal_setup)?;
    Ok(Self::new(period, apodization))
  }

  /// Get the optimal periodic poling for the given signal and pump beams
  pub fn try_as_optimum(
    self,
    signal: &SignalBeam,
    pump: &PumpBeam,
    crystal_setup: &CrystalSetup,
  ) -> Result<Self, SPDCError> {
    match self {
      Self::Off => Self::try_new_optimum(signal, pump, crystal_setup, Apodization::Off),
      Self::On { apodization, .. } => {
        Self::try_new_optimum(signal, pump, crystal_setup, apodization)
      }
    }
  }

  /// calculate the sign of periodic poling
  pub fn compute_sign(signal: &SignalBeam, pump: &PumpBeam, crystal_setup: &CrystalSetup) -> Sign {
    let idler =
      IdlerBeam::try_new_optimum(signal, pump, crystal_setup, PeriodicPoling::Off).unwrap();
    let delkz = (delta_k(
      signal.frequency(),
      idler.frequency(),
      signal,
      &idler,
      pump,
      crystal_setup,
      PeriodicPoling::Off,
    ) * M
      / RAD)
      .z;

    // converts to sign
    delkz.into()
  }

  /// Set the apodization
  ///
  /// Does nothing if off
  pub fn set_apodization(&mut self, apodization: Apodization) {
    match self {
      Self::Off => {}
      Self::On { period, sign, .. } => *self = Self::new(*sign * *period, apodization),
    }
  }

  /// Set the apodization
  ///
  /// Does nothing if off
  pub fn with_apodization(self, apodization: Apodization) -> Self {
    match self {
      Self::Off => Self::Off,
      Self::On { period, sign, .. } => Self::new(sign * period, apodization),
    }
  }

  /// Get the apodization
  pub fn apodization(&self) -> &Apodization {
    match self {
      Self::Off => &Apodization::Off,
      Self::On { apodization, .. } => apodization,
    }
  }

  /// Get the signed period
  pub fn signed_period(&self) -> PolingPeriod {
    match self {
      Self::Off => f64::INFINITY * M,
      Self::On { period, sign, .. } => *sign * *period,
    }
  }

  /// Get the contribution of periodic poling to to delta_k
  pub fn k_eff(&self) -> Wavenumber {
    let m = 1.; // TODO: QPM order
    match self {
      Self::Off => 0. * RAD / M,
      &Self::On { period, sign, .. } => {
        assert!(
          period.value_unsafe > 0.,
          "Periodic Poling Period must be greater than zero"
        );
        TWO_PI * RAD * m / (sign * period)
      }
    }
  }

  /// Get the apodization integration constant for the given z position
  pub fn integration_constant(&self, z: f64, crystal_length: Distance) -> f64 {
    match &self {
      Self::Off => 1.,
      &Self::On { apodization, .. } => apodization.integration_constant(z, crystal_length),
    }
  }
}

/// Get the optimum poling period for the given signal and pump beams
pub fn optimum_poling_period(
  signal: &SignalBeam,
  pump: &PumpBeam,
  crystal_setup: &CrystalSetup,
) -> Result<PolingPeriod, SPDCError> {
  // z component of delta k, based on periodic poling
  let delta_kz = |pp| {
    let idler = IdlerBeam::try_new_optimum(signal, pump, crystal_setup, &pp).unwrap();
    let del_k = delta_k(
      signal.frequency(),
      idler.frequency(),
      signal,
      &idler,
      pump,
      crystal_setup,
      &pp,
    );

    let del_k_vec = *(del_k * M / RAD);

    del_k_vec.z
  };

  let z = delta_kz(PeriodicPoling::Off);

  if z == 0. {
    // z is already zero, that means there is already perfect phasematching
    // no poling period needed
    return Ok(f64::INFINITY * M);
  }

  // base our guess on the delta k calculation without periodic poling
  let guess = TWO_PI / z;
  // the sign of the z component of delta k gives the sign of pp
  let sign = z.into();

  // let spdc = Mutex::new(SPDC::new(
  //   crystal_setup.clone(),
  //   signal.clone(),
  //   IdlerBeam::try_new_optimum(signal, pump, crystal_setup, PeriodicPoling::Off).unwrap(),
  //   pump.clone(),
  //   5e-9 * M,
  //   1e-3 * W,
  //   1e-2,
  //   PeriodicPoling::Off,
  //   0. * M,
  //   0. * M,
  //   1e-12 * M / V,
  // ).with_optimal_waist_positions());

  // minimizable delta k function based on period (using predetermined sign)
  let pm = |period| {
    let pp = PeriodicPoling::On {
      period: period * M,
      sign,
      apodization: Apodization::Off,
    };

    delta_kz(pp).abs()

    // THIS was an attempt to use JSI at middle to optimize
    // the periodic poling. Since the min deltakz wasn't
    // lining up with maximum JSI.
    // let idler = IdlerBeam::try_new_optimum(signal, pump, crystal_setup, &pp).unwrap();
    // let wi = idler.frequency();

    // let mut spdc = spdc.lock().unwrap();
    // spdc.idler = idler;
    // spdc.pp = pp;

    // -(phasematch_fiber_coupling(
    //   signal.frequency(),
    //   wi,
    //   &spdc,
    //   None
    // ) / PerMeter4::new(1.)).norm_sqr()
  };

  // maximum period is the length of the crystal
  let max_period = *(crystal_setup.length / M);
  // minimum period... typical poling periods are on the order of microns
  let min_period = std::f64::MIN_POSITIVE;

  // minimize...
  let period = nelder_mead_1d(
    pm,
    (guess.abs(), guess.abs() + 1e-6),
    1000,
    min_period,
    max_period,
    1e-12,
  );

  if max_period < period || period < min_period {
    Err(SPDCError::new(IMPOSSIBLE_POLING_PERIOD.to_string()))
  } else {
    Ok(sign * period * M)
  }
}

#[cfg(feature = "pyo3")]
mod pyo3_impls {
  use super::*;
  use dim::f64prefixes::*;
  use pyo3::{
    exceptions::{PyKeyError, PyValueError},
    prelude::*,
    types::PyDict,
  };

  impl FromPyObject<'_> for Apodization {
    fn extract_bound(obj: &Bound<'_, PyAny>) -> PyResult<Self> {
      let dict = obj.downcast::<PyDict>()?;
      let kind = dict
        .get_item("kind")?
        .ok_or(PyKeyError::new_err(
          "Expecting 'kind' key in apodization config.",
        ))?
        .extract::<String>()?;

      let param = dict.get_item("parameter")?.ok_or(PyKeyError::new_err(
        "Expecting 'parameter' key in apodization config.",
      ))?;

      match kind.to_lowercase().as_str() {
        "off" => Ok(Apodization::Off),
        "gaussian" => {
          let fwhm = param
            .downcast::<PyDict>()?
            .get_item("fwhm_um")?
            .ok_or(PyKeyError::new_err(
              "Expecting dictionary with 'fwhm_um' key as gaussian apodization parameter value.",
            ))?
            .extract::<f64>()?;
          Ok(Apodization::Gaussian {
            fwhm: fwhm * MICRO * M,
          })
        }
        "interpolate" => {
          let points = param.extract::<Vec<f64>>()?;
          Ok(Apodization::Interpolate(points))
        }
        _ => serde_json::from_value(serde_json::json! {
          {
            "kind": kind,
            "parameter": param.extract::<f64>()?
          }
        })
        .map_err(|e| PyValueError::new_err(e.to_string())),
      }
    }
  }

  impl ToPyObject for Apodization {
    fn to_object(&self, py: Python<'_>) -> PyObject {
      let dict = PyDict::new_bound(py);

      match self {
        &Apodization::Off => {
          dict.set_item("kind", "off").unwrap();
        }
        &Apodization::Gaussian { fwhm } => {
          let fwhm_um = *(fwhm / MICRO / M);
          let param = PyDict::new_bound(py);
          param.set_item("fwhm_um", fwhm_um).unwrap();
          dict.set_item("kind", "gaussian").unwrap();
          dict.set_item("parameter", param).unwrap();
        }
        Apodization::Interpolate(points) => {
          dict.set_item("kind", "interpolate").unwrap();
          dict.set_item("parameter", points.clone()).unwrap();
        }
        &Apodization::Bartlett(p)
        | &Apodization::Blackman(p)
        | &Apodization::Connes(p)
        | &Apodization::Cosine(p)
        | &Apodization::Hamming(p)
        | &Apodization::Welch(p) => {
          dict.set_item("kind", self.kind()).unwrap();
          dict.set_item("parameter", p).unwrap();
        }
      }

      dict.into()
    }
  }

  impl IntoPy<PyObject> for Apodization {
    fn into_py(self, py: Python<'_>) -> PyObject {
      self.to_object(py).into_py(py)
    }
  }
}

#[cfg(test)]
mod test {
  use crate::beam::*;
  use crate::utils::testing::*;
  use crate::{dim::ucum::*, *};

  #[test]
  fn test_poling_period() {
    let signal = Beam::new(
      PolarizationType::Extraordinary,
      0. * RAD,
      0. * RAD,
      1550e-9 * M,
      30e-6 * M,
    )
    .into();

    let pump = Beam::new(
      PolarizationType::Extraordinary,
      0. * RAD,
      0. * RAD,
      775e-9 * M,
      100e-6 * M,
    )
    .into();

    let mut crystal_setup: CrystalSetup = CrystalConfig::default().into();
    crystal_setup.theta = 90. * DEG;
    crystal_setup.length = 20_000e-6 * M;

    let period = optimum_poling_period(&signal, &pump, &crystal_setup).unwrap();

    assert_nearly_equal!("poling period", *(period / M), -46.578592559e-6, 0.2);
  }

  #[test]
  fn test_poling_domains() {
    let crystal_length = 1000e-6 * M;
    let pp = PeriodicPoling::new(10e-6 * M, Apodization::Off);
    let domains = pp.poling_domains(crystal_length);
    assert_eq!(domains, vec![(0.5, 0.5); 100]);
  }

  #[test]
  fn test_poling_domains_step() {
    let crystal_length = 100e-6 * M;
    let step_function = vec![0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1.];
    let pp = PeriodicPoling::new(10e-6 * M, Apodization::Interpolate(step_function));
    let domains = pp.poling_domains(crystal_length);
    assert_eq!(
      domains,
      vec![
        (0., 1.),
        (0., 1.),
        (0., 1.),
        (0., 1.),
        (0., 1.),
        (0.5, 0.5),
        (0.5, 0.5),
        (0.5, 0.5),
        (0.5, 0.5),
        (0.5, 0.5)
      ]
    );
  }
}
