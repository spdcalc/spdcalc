use crate::{SPDC, jsa::FrequencySpace};
use dim::ucum::{S, Hertz};

pub struct Efficiencies {
  pub symmetric: f64,
  pub signal: f64,
  pub idler: f64
}

pub fn efficiencies(spdc : &SPDC, ranges: FrequencySpace, integration_steps: Option<usize>) -> Efficiencies {
  let coincidences_rate = spdc.counts_coincidences(ranges, integration_steps);
  let signal_singles_rate = spdc.counts_singles_signal(ranges, integration_steps);
  let idler_singles_rate = spdc.counts_singles_idler(ranges, integration_steps);
  let signal_efficiency = if idler_singles_rate == Hertz::new(0.) { 0. } else { *(coincidences_rate / idler_singles_rate) };
  let idler_efficiency = if signal_singles_rate == Hertz::new(0.) { 0. } else { *(coincidences_rate / signal_singles_rate) };
  let symmetric_efficiency = if signal_singles_rate == Hertz::new(0.) || idler_singles_rate == Hertz::new(0.) {
    0.
  } else {
    let denom : f64 = *(signal_singles_rate * idler_singles_rate * S * S);
    *(coincidences_rate * S / denom.sqrt())
  };
  Efficiencies {
    symmetric: symmetric_efficiency,
    signal: signal_efficiency,
    idler: idler_efficiency
  }
}