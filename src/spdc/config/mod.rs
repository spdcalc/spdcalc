use super::*;
use crate::{
  beam::Beam,
  beam::IdlerBeam,
  beam::PumpBeam,
  beam::SignalBeam,
  dim::{
    f64prefixes::{MICRO, NANO, PICO},
    ucum::{DEG, M, MILLIW, RAD, V},
  },
  math::sigfigs,
  utils::{self, from_kelvin_to_celsius},
  CrystalSetup, CrystalType, PMType, SPDCError,
};

mod periodic_poling_config;
pub use periodic_poling_config::PeriodicPolingConfig;
mod apodization;
pub use apodization::ApodizationConfig;
use serde_with::{serde_as, DisplayFromStr};

const SIG_FIGS_IN_CONFIG: u8 = 4;

/// A config parameter the could be signaled to be automatically calculated
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(untagged)]
pub enum AutoCalcParam<T>
where
  T: 'static,
{
  /// "auto"
  Auto(String),
  /// A parameter
  Param(T),
}

impl<T> AutoCalcParam<T> {
  /// Returns true if the parameter is set to be automatically calculated
  pub fn is_auto(&self) -> bool {
    matches!(self, Self::Auto(_))
  }
}

impl<T> Default for AutoCalcParam<T> {
  fn default() -> Self {
    Self::Auto("auto".into())
  }
}

impl<T> From<AutoCalcParam<T>> for Option<T> {
  fn from(param: AutoCalcParam<T>) -> Self {
    match param {
      AutoCalcParam::Param(x) => Some(x),
      AutoCalcParam::Auto(_) => None,
    }
  }
}

impl<T> From<Option<T>> for AutoCalcParam<T> {
  fn from(param: Option<T>) -> Self {
    match param {
      Some(x) => AutoCalcParam::Param(x),
      None => AutoCalcParam::Auto("auto".into()),
    }
  }
}

/// Flat configuration of crystal for ease of import/export
#[serde_as]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CrystalConfig {
  /// The type of crystal
  pub kind: CrystalType,
  /// The type of phase matching
  #[serde_as(as = "DisplayFromStr")]
  pub pm_type: PMType,
  /// The polar angle of the crystal in degrees
  #[serde(default)]
  pub phi_deg: f64,
  /// The azimuthal angle of the crystal in degrees
  #[serde(default)]
  pub theta_deg: AutoCalcParam<f64>,
  /// The length of the crystal in micrometers
  pub length_um: f64,
  /// The temperature of the crystal in degrees Celsius
  pub temperature_c: f64,
  /// Whether the downconversion is counter-propagating
  #[serde(default)]
  pub counter_propagation: bool,
}

impl Default for CrystalConfig {
  fn default() -> Self {
    Self {
      kind: CrystalType::KTP,
      pm_type: PMType::Type2_e_eo,
      phi_deg: 0.,
      theta_deg: AutoCalcParam::Auto("auto".into()),
      length_um: 2_000.,
      temperature_c: 20.,
      counter_propagation: false,
    }
  }
}

/// Converts without autocalculating theta
impl From<CrystalConfig> for CrystalSetup {
  fn from(cfg: CrystalConfig) -> Self {
    let theta = if let AutoCalcParam::Param(theta_deg) = cfg.theta_deg {
      theta_deg * DEG
    } else {
      0. * DEG
    };
    CrystalSetup {
      crystal: cfg.kind,
      pm_type: cfg.pm_type,
      phi: cfg.phi_deg * DEG,
      theta,
      length: cfg.length_um * MICRO * M,
      temperature: utils::from_celsius_to_kelvin(cfg.temperature_c),
      counter_propagation: cfg.counter_propagation,
    }
  }
}

impl From<CrystalSetup> for CrystalConfig {
  fn from(setup: CrystalSetup) -> Self {
    Self {
      kind: setup.crystal,
      pm_type: setup.pm_type,
      theta_deg: AutoCalcParam::Param(sigfigs(*(setup.theta / DEG), SIG_FIGS_IN_CONFIG)),
      phi_deg: sigfigs(*(setup.phi / DEG), SIG_FIGS_IN_CONFIG),
      length_um: sigfigs(*(setup.length / (MICRO * M)), SIG_FIGS_IN_CONFIG),
      temperature_c: sigfigs(
        from_kelvin_to_celsius(setup.temperature),
        SIG_FIGS_IN_CONFIG,
      ),
      counter_propagation: setup.counter_propagation,
    }
  }
}

/// Flat configuration of pump for ease of import/export
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PumpConfig {
  /// The wavelength of the pump in nanometers
  pub wavelength_nm: f64,
  /// The waist of the pump in micrometers
  pub waist_um: f64,
  /// The bandwidth of the pump in nanometers
  pub bandwidth_nm: f64,
  /// The average power of the pump in milliwatts
  pub average_power_mw: f64,
  /// If unset, defaults to 1e-2
  pub spectrum_threshold: Option<f64>,
}

impl Default for PumpConfig {
  fn default() -> Self {
    Self {
      wavelength_nm: 775.,
      waist_um: 100.,
      bandwidth_nm: 5.53,
      average_power_mw: 1.,
      spectrum_threshold: Some(1e-2),
    }
  }
}

impl PumpConfig {
  /// Converts to a [`PumpBeam`] with the given [`CrystalSetup`]
  pub fn as_beam(self, crystal_setup: &CrystalSetup) -> PumpBeam {
    Beam::new(
      crystal_setup.pm_type.pump_polarization(),
      0. * RAD,
      0. * RAD,
      self.wavelength_nm * NANO * M,
      self.waist_um * MICRO * M,
    )
    .into()
  }
}

/// Flat configuration of signal for ease of import/export
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SignalConfig {
  /// The wavelength of the signal in nanometers
  pub wavelength_nm: f64,
  /// The polar angle of the signal in degrees
  #[serde(default)]
  pub phi_deg: f64,
  // one of...
  /// The azimuthal angle of the signal in degrees (conflicts with `theta_external_deg`)
  pub theta_deg: Option<f64>,
  /// The azimuthal angle of the signal in degrees outside the crystal (conflicts with `theta_deg`)
  pub theta_external_deg: Option<f64>,
  /// The waist of the signal in micrometers
  pub waist_um: f64,
  /// The position of the waist of the signal in micrometers
  #[serde(default)]
  pub waist_position_um: AutoCalcParam<f64>,
}

impl Default for SignalConfig {
  fn default() -> Self {
    Self {
      wavelength_nm: 1550.,
      phi_deg: 0.,
      theta_deg: Some(0.),
      theta_external_deg: None,
      waist_um: 100.,
      waist_position_um: AutoCalcParam::default(),
    }
  }
}

impl SignalConfig {
  /// Converts to a [`SignalBeam`] with the given [`CrystalSetup`]
  pub fn try_as_beam(self, crystal_setup: &CrystalSetup) -> Result<SignalBeam, SPDCError> {
    let phi = self.phi_deg * DEG;
    let mut beam = Beam::new(
      crystal_setup.pm_type.signal_polarization(),
      phi,
      0. * RAD,
      self.wavelength_nm * NANO * M,
      self.waist_um * MICRO * M,
    );
    match (self.theta_deg, self.theta_external_deg) {
      (Some(theta), None) => beam.set_angles(phi, theta * DEG),
      (None, Some(theta_e)) => beam.set_theta_external(theta_e * DEG, crystal_setup),
      _ => {
        return Err(SPDCError(
          "Must specify one of theta_deg or theta_external_deg".into(),
        ))
      }
    };

    Ok(beam.into())
  }
}

/// Flat configuration of idler for ease of import/export
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct IdlerConfig {
  /// The wavelength of the idler in nanometers
  pub wavelength_nm: f64,
  /// The polar angle of the idler in degrees
  #[serde(default)]
  pub phi_deg: f64,
  // one of...
  /// The azimuthal angle of the idler in degrees (conflicts with `theta_external_deg`)
  pub theta_deg: Option<f64>,
  /// The azimuthal angle of the idler in degrees outside the crystal (conflicts with `theta_deg`)
  pub theta_external_deg: Option<f64>,
  /// The waist of the idler in micrometers
  pub waist_um: f64,
  /// The position of the waist of the idler in micrometers
  #[serde(default)]
  pub waist_position_um: AutoCalcParam<f64>,
}

impl IdlerConfig {
  /// Converts to a [`IdlerBeam`] with the given [`CrystalSetup`]
  pub fn try_as_beam(self, crystal_setup: &CrystalSetup) -> Result<IdlerBeam, SPDCError> {
    let phi = self.phi_deg * DEG;
    let mut beam = Beam::new(
      crystal_setup.pm_type.idler_polarization(),
      phi,
      0. * RAD,
      self.wavelength_nm * NANO * M,
      self.waist_um * MICRO * M,
    );
    match (self.theta_deg, self.theta_external_deg) {
      (Some(theta), None) => beam.set_angles(phi, theta * DEG),
      (None, Some(theta_e)) => beam.set_theta_external(theta_e * DEG, crystal_setup),
      _ => {
        return Err(SPDCError(
          "Must specify one of theta_deg or theta_external_deg".into(),
        ))
      }
    };

    Ok(beam.into())
  }
}

/// Flat configuration of SPDC for ease of import/export
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SPDCConfig {
  /// The crystal configuration
  pub crystal: CrystalConfig,
  /// The pump configuration
  pub pump: PumpConfig,
  /// The signal configuration
  pub signal: SignalConfig,
  /// The idler configuration
  #[serde(default)]
  pub idler: AutoCalcParam<IdlerConfig>,
  /// The configuration for the periodic poling
  #[serde(default)]
  pub periodic_poling: PeriodicPolingConfig,
  /// The deff of the crystal in pm/V
  pub deff_pm_per_volt: f64,
}

impl TryFrom<SPDCConfig> for SPDC {
  type Error = SPDCError;

  fn try_from(config: SPDCConfig) -> Result<Self, Self::Error> {
    config.try_as_spdc()
  }
}

impl SPDCConfig {
  /// Converts to a [`SPDC`] with the given [`CrystalSetup`]
  pub fn try_as_spdc(self) -> Result<SPDC, SPDCError> {
    let deff = self.deff_pm_per_volt * PICO * M / V;
    let pump_spectrum_threshold = self.pump.spectrum_threshold.unwrap_or(1e-2);
    let crystal_theta_autocalc = self.crystal.theta_deg.is_auto();
    let signal_waist_position_um = self.signal.waist_position_um.clone();
    let pump_bandwidth = self.pump.bandwidth_nm * NANO * M;
    let pump_average_power = self.pump.average_power_mw * MILLIW;
    let mut crystal_setup: CrystalSetup = self.crystal.into();
    let pump = self.pump.as_beam(&crystal_setup);
    let signal = self.signal.try_as_beam(&crystal_setup)?;
    let periodic_poling =
      self
        .periodic_poling
        .try_as_periodic_poling(&signal, &pump, &crystal_setup)?;

    if crystal_theta_autocalc {
      if periodic_poling == PeriodicPoling::Off {
        crystal_setup.assign_optimum_theta(&signal, &pump);
      } else {
        return Err(SPDCError("Can not autocalc theta when periodic poling is enabled. Provide an explicit value for crystal theta.".into()));
      }
    }
    let idler = match &self.idler {
      AutoCalcParam::Param(idler_cfg) => idler_cfg.clone().try_as_beam(&crystal_setup)?,
      AutoCalcParam::Auto(_) => {
        IdlerBeam::try_new_optimum(&signal, &pump, &crystal_setup, &periodic_poling)?
      }
    };
    let idler_waist_position = {
      let auto = AutoCalcParam::default();
      let autocalc_idler_waist = match self.idler {
        AutoCalcParam::Param(cfg) => cfg.waist_position_um,
        AutoCalcParam::Auto(_) => auto,
      };
      match autocalc_idler_waist {
        AutoCalcParam::Param(focus_um) => -focus_um.abs() * MICRO * M,
        AutoCalcParam::Auto(_) => {
          crystal_setup.optimal_waist_position(idler.vacuum_wavelength(), idler.polarization())
        }
      }
    };
    let signal_waist_position = match signal_waist_position_um {
      AutoCalcParam::Param(focus_um) => -focus_um.abs() * MICRO * M,
      AutoCalcParam::Auto(_) => {
        crystal_setup.optimal_waist_position(signal.vacuum_wavelength(), signal.polarization())
      }
    };

    Ok(SPDC::new(
      crystal_setup,
      signal,
      idler,
      pump,
      pump_bandwidth,
      pump_average_power,
      pump_spectrum_threshold,
      periodic_poling,
      signal_waist_position,
      idler_waist_position,
      deff,
    ))
  }
}

impl Default for SPDCConfig {
  fn default() -> Self {
    Self {
      crystal: CrystalConfig::default(),
      pump: PumpConfig::default(),
      signal: SignalConfig::default(),
      idler: AutoCalcParam::default(),
      periodic_poling: PeriodicPolingConfig::default(),
      deff_pm_per_volt: 1.,
    }
  }
}

impl From<SPDC> for PumpConfig {
  fn from(spdc: SPDC) -> Self {
    Self {
      wavelength_nm: sigfigs(
        *(spdc.pump.vacuum_wavelength() / (NANO * M)),
        SIG_FIGS_IN_CONFIG,
      ),
      bandwidth_nm: sigfigs(*(spdc.pump_bandwidth / (NANO * M)), SIG_FIGS_IN_CONFIG),
      waist_um: sigfigs(*(spdc.pump.waist().x / (MICRO * M)), SIG_FIGS_IN_CONFIG),
      average_power_mw: sigfigs(*(spdc.pump_average_power / MILLIW), SIG_FIGS_IN_CONFIG),
      spectrum_threshold: Some(spdc.pump_spectrum_threshold),
    }
  }
}

impl From<SPDC> for SignalConfig {
  fn from(spdc: SPDC) -> Self {
    Self {
      wavelength_nm: sigfigs(
        *(spdc.signal.vacuum_wavelength() / (NANO * M)),
        SIG_FIGS_IN_CONFIG,
      ),
      theta_deg: Some(sigfigs(
        *(spdc.signal.theta_internal() / DEG),
        SIG_FIGS_IN_CONFIG,
      )),
      theta_external_deg: None,
      phi_deg: sigfigs(*(spdc.signal.phi() / DEG), SIG_FIGS_IN_CONFIG),
      waist_um: sigfigs(*(spdc.signal.waist().x / (MICRO * M)), SIG_FIGS_IN_CONFIG),
      waist_position_um: AutoCalcParam::Param(sigfigs(
        *(spdc.signal_waist_position / (MICRO * M)),
        SIG_FIGS_IN_CONFIG,
      )),
    }
  }
}

impl From<SPDC> for IdlerConfig {
  fn from(spdc: SPDC) -> Self {
    Self {
      wavelength_nm: sigfigs(
        *(spdc.idler.vacuum_wavelength() / (NANO * M)),
        SIG_FIGS_IN_CONFIG,
      ),
      theta_deg: Some(sigfigs(
        *(spdc.idler.theta_internal() / DEG),
        SIG_FIGS_IN_CONFIG,
      )),
      theta_external_deg: None,
      phi_deg: sigfigs(*(spdc.idler.phi() / DEG), SIG_FIGS_IN_CONFIG),
      waist_um: sigfigs(*(spdc.idler.waist().x / (MICRO * M)), SIG_FIGS_IN_CONFIG),
      waist_position_um: AutoCalcParam::Param(*(spdc.idler_waist_position / (MICRO * M))),
    }
  }
}

impl From<SPDC> for SPDCConfig {
  fn from(spdc: SPDC) -> Self {
    let crystal = spdc.crystal_setup.into();
    let pump = PumpConfig {
      wavelength_nm: sigfigs(
        *(spdc.pump.vacuum_wavelength() / (NANO * M)),
        SIG_FIGS_IN_CONFIG,
      ),
      bandwidth_nm: sigfigs(*(spdc.pump_bandwidth / (NANO * M)), SIG_FIGS_IN_CONFIG),
      waist_um: sigfigs(*(spdc.pump.waist().x / (MICRO * M)), SIG_FIGS_IN_CONFIG),
      average_power_mw: sigfigs(*(spdc.pump_average_power / MILLIW), SIG_FIGS_IN_CONFIG),
      spectrum_threshold: Some(spdc.pump_spectrum_threshold),
    };
    let signal = SignalConfig {
      wavelength_nm: sigfigs(
        *(spdc.signal.vacuum_wavelength() / (NANO * M)),
        SIG_FIGS_IN_CONFIG,
      ),
      theta_deg: Some(sigfigs(
        *(spdc.signal.theta_internal() / DEG),
        SIG_FIGS_IN_CONFIG,
      )),
      theta_external_deg: None,
      phi_deg: sigfigs(*(spdc.signal.phi() / DEG), SIG_FIGS_IN_CONFIG),
      waist_um: sigfigs(*(spdc.signal.waist().x / (MICRO * M)), SIG_FIGS_IN_CONFIG),
      waist_position_um: AutoCalcParam::Param(sigfigs(
        *(spdc.signal_waist_position / (MICRO * M)),
        SIG_FIGS_IN_CONFIG,
      )),
    };

    let idler = AutoCalcParam::Param(IdlerConfig {
      wavelength_nm: sigfigs(
        *(spdc.idler.vacuum_wavelength() / (NANO * M)),
        SIG_FIGS_IN_CONFIG,
      ),
      theta_deg: Some(sigfigs(
        *(spdc.idler.theta_internal() / DEG),
        SIG_FIGS_IN_CONFIG,
      )),
      theta_external_deg: None,
      phi_deg: sigfigs(*(spdc.idler.phi() / DEG), SIG_FIGS_IN_CONFIG),
      waist_um: sigfigs(*(spdc.idler.waist().x / (MICRO * M)), SIG_FIGS_IN_CONFIG),
      waist_position_um: AutoCalcParam::Param(*(spdc.idler_waist_position / (MICRO * M))),
    });

    let periodic_poling = spdc.pp.into();
    let deff_pm_per_volt = sigfigs(*(spdc.deff / (PICO * M / V)), SIG_FIGS_IN_CONFIG);

    Self {
      crystal,
      pump,
      signal,
      idler,
      periodic_poling,
      deff_pm_per_volt,
    }
  }
}

#[cfg(test)]
mod test {
  use super::*;
  use crate::utils::testing::*;
  use crate::PolarizationType;
  use dim::ucum;
  use float_cmp::approx_eq;
  use serde_json::json;

  #[test]
  fn auto_idler_auto_focus_test() {
    let json = json!({
      "crystal": {
        "kind": "BBO_1",
        "pm_type": "e->eo",
        "phi_deg": 0,
        "theta_deg": 0,
        "length_um": 2000,
        "temperature_c": 20
      },
      "pump": {
        "wavelength_nm": 775,
        "waist_um": 100,
        "bandwidth_nm": 5.35,
        "average_power_mw": 1
      },
      "signal": {
        "wavelength_nm": 1550,
        "phi_deg": 0,
        "theta_external_deg": 0,
        "waist_um": 100,
        "waist_position_um": "auto"
      },
      "idler": "auto",
      "deff_pm_per_volt": 1
    });

    let config: SPDCConfig = serde_json::from_value(json).expect("Could not unwrap json");
    let actual = config
      .try_as_spdc()
      .expect("Could not convert to SPDC instance");

    let expected = {
      let crystal_setup = CrystalSetup {
        crystal: CrystalType::BBO_1,
        pm_type: PMType::Type2_e_eo,
        phi: 0. * RAD,
        theta: 0. * RAD,
        length: 2000. * MICRO * M,
        temperature: 293.15 * ucum::K,
        counter_propagation: false,
      };
      let signal = Beam::new(
        PolarizationType::Extraordinary,
        0. * DEG,
        0. * RAD,
        1550. * NANO * M,
        100. * MICRO * M,
      )
      .into();
      let idler = Beam::new(
        PolarizationType::Ordinary,
        180. * DEG,
        0. * DEG,
        1550. * NANO * M,
        100. * MICRO * M,
      )
      .into();
      let pump = Beam::new(
        PolarizationType::Extraordinary,
        0. * DEG,
        0. * DEG,
        775. * NANO * M,
        100. * MICRO * M,
      )
      .into();
      let pump_average_power = 1. * MILLIW;
      let pump_bandwidth = 5.35 * NANO * M;
      let periodic_poling = PeriodicPoling::Off;
      let signal_waist_position = -0.0006073170564963904 * M;
      let idler_waist_position = -0.0006073170564963904 * M;
      let deff = 1. * PICO * M / V;
      let pump_spectrum_threshold = 1e-2;
      SPDC::new(
        crystal_setup,
        signal,
        idler,
        pump,
        pump_bandwidth,
        pump_average_power,
        pump_spectrum_threshold,
        periodic_poling,
        signal_waist_position,
        idler_waist_position,
        deff,
      )
    };
    dbg!(&actual);
    assert_eq!(actual, expected);
  }

  #[test]
  fn autocalc_crystal_theta_test() {
    let json = json!({
      "crystal": {
        "kind": "KTP",
        "pm_type": "e->eo",
        "phi_deg": 0,
        "theta_deg": "auto",
        "length_um": 2000,
        "temperature_c": 20
      },
      "pump": {
        "wavelength_nm": 775,
        "bandwidth_nm": 5.35,
        "waist_um": 100,
        "average_power_mw": 1
      },
      "signal": {
        "wavelength_nm": 1550,
        "phi_deg": 0,
        "theta_external_deg": 0,
        "waist_um": 100,
        "waist_position_um": "auto"
      },
      "idler": "auto",
      "deff_pm_per_volt": 7.6,
    });

    let config: SPDCConfig = serde_json::from_value(json).expect("Could not unwrap json");
    let spdc = config
      .try_as_spdc()
      .expect("Could not convert to SPDC instance");

    let actual = spdc.crystal_setup.theta;
    let expected = spdc.crystal_setup.optimum_theta(&spdc.signal, &spdc.pump);

    assert_nearly_equal!("crystal theta", actual.value_unsafe, expected.value_unsafe);
  }

  #[test]
  fn auto_pp_test() {
    let json = json!({
      "crystal": {
        "kind": "KTP",
        "pm_type": "e->eo",
        "phi_deg": 0,
        "theta_deg": 90,
        "length_um": 2000,
        "temperature_c": 20
      },
      "pump": {
        "wavelength_nm": 775,
        "bandwidth_nm": 5.35,
        "waist_um": 100,
        "average_power_mw": 1
      },
      "signal": {
        "wavelength_nm": 1550,
        "phi_deg": 0,
        "theta_external_deg": 0,
        "waist_um": 100,
        "waist_position_um": "auto"
      },
      "idler": "auto",
      "periodic_poling": {
        "poling_period_um": "auto",
      },
      "deff_pm_per_volt": 7.6,
    });

    let config: SPDCConfig = serde_json::from_value(json).expect("Could not unwrap json");
    let spdc = config
      .try_as_spdc()
      .expect("Could not convert to SPDC instance");

    let actual = match spdc.pp {
      PeriodicPoling::Off => panic!("Periodic poling should be on"),
      PeriodicPoling::On { period, .. } => period,
    };
    let expected = 46.52987937008885 * MICRO * M;

    assert!(approx_eq!(
      f64,
      actual.value_unsafe,
      expected.value_unsafe,
      epsilon = 1e-6
    ));
  }

  #[test]
  fn spdc_config_test_interpolated_crystal() {
    let json_cfg = json!({
      "crystal": {
        "kind": {
          "name": "InterpolatedUniaxial",
          "wavelengths_nm": [400.0, 500.0, 600.0],
          "no": [1.66, 1.65, 1.64],
          "ne": [1.55, 1.54, 1.53]
        },
        "pm_type": "Type2_e_eo",
        "theta_deg": 90.0,
        "length_um": 2000.0,
        "temperature_c": 20.0,
        "counter_propagation": false
      },
      "pump": {
        "wavelength_nm": 775.0,
        "waist_um": 100.0,
        "bandwidth_nm": 5.35,
        "average_power_mw": 1.0,
        "spectrum_threshold": 0.01
      },
      "signal": {
        "wavelength_nm": 1550.0,
        "phi_deg": 0.0,
        "theta_deg": 0.0,
        "waist_um": 100.0
      },
      "idler": "auto",
      "periodic_poling": {
        "poling_period_um": "auto"
      },
      "deff_pm_per_volt": 1.0
    });
    let config: SPDCConfig = serde_json::from_value(json_cfg).unwrap();
    let result = config.try_as_spdc();
    assert!(result.is_ok());
  }
}
