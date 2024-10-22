# SPDCalc Library

![build](https://github.com/spdcalc/spdcalc/actions/workflows/ci.yaml/badge.svg?branch=master)
![Crates.io Version](https://img.shields.io/crates/v/spdcalc)

SPDCalc is a fast design tool for spontaneous parametric downconversion sources.

This is the base rust library used by the [web app](https://app.spdcalc.org/)
and the python library [spdcalc-py](https://pypi.org/project/spdcalc-py/).

# Introduction

The process of spontaneous parametric downconversion is an important source
of single-photons and quantum states of light. Much recent theoretical and
experimental work focuses on engineering the properties of the photons
emitted from downconversion. For multi-photon experiments generating photon
pairs that are spectrally pure, meaning the signal and idler photons are
uncorrelated, has been an active area of investigation. For fundamental
tests of quantum mechanics, such as loop-hole-free Bell tests, optimizing
the coupling of photons to single-mode fiber is an important consideration.
//!
While the principles behind the nonlinear process of downconversion are well
understood, calculating the properties of the photons can be difficult and
time-consuming when considering photons emitted at an angle with respect
to the pump (non collinear) when the signal, idler, and pump can experience
birefringent walk off. Here we present a library designed to simplify these
calculations that builds on previous work. The library allows
a user to easily determine the phasematching properties as well as the emission
characteristics and joint spectrum of the signal and idler photons. The library
also allows a user to quickly find phasematching parameters that generate
spectrally pure states. The Hong-Ou-Mandel dip between the signal and idler
is also calculated, and the joint spectrum in the dip can also be easily
examined. Finally, the spatial mode structure and single-mode fiber coupling
efficiency can be computed. The library is designed to be flexible, and can
handle periodically-poled crystals that have been apodized with a Gaussian
function.

**There is also a web-based application that uses this library available at
[app.spdcalc.org](https://app.spdcalc.org).**

# Usage

The easiest way to use this library is to use the [SPDCConfig] struct.
Either create a default setup with `SPDCConfig::default()` or create a
custom setup by specifying the parameters you need. This is made to be
used with [serde](https://serde.rs/), so you can deserialize a config from a JSON file or YAML
or similar.

```rs
use serde_json::json;
use spdcalc::prelude::*;

let json = serde_json::json!({
  "crystal": {
    "kind": "KTP",
    "pm_type": "e->eo",
    "phi_deg": 0,
    "theta_deg": 90,
    "length_um": 14000,
    "temperature_c": 20
  },
  "pump": {
    "wavelength_nm": 775,
    "waist_um": 200,
    "bandwidth_nm": 0.5,
    "average_power_mw": 300
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
    "poling_period_um": "auto"
  },
  "deff_pm_per_volt": 7.6
});

let config : SPDCConfig = serde_json::from_value(json).expect("Could not unwrap json");
let spdc = config.try_as_spdc().expect("Could not convert to SPDC instance");

// Now you can use the spdc instance to calculate things like the jsa

let range = WavelengthSpace::new(
  (1400e-9 * M, 1600e-9 * M, 100),
  (1400e-9 * M, 1600e-9 * M, 100)
);
let spectrum = spdc.joint_spectrum(Integrator::default());
let results = spectrum.jsa_range(range);
```

## Development

```sh
cargo build
```
