use dimensioned::f64prefixes::{MICRO, NANO};
use spdcalc::{
  beam::Beam,
  beam::IdlerBeam,
  beam::SignalBeam,
  dim::ucum::{DEG, K, M},
  Apodization, CrystalSetup, CrystalType, PMType, PeriodicPoling, PolarizationType,
};

fn main() {
  let signal: SignalBeam = Beam::new(
    PolarizationType::Extraordinary,
    0. * DEG,
    180. * DEG,
    1550e-9 * M,
    100e-6 * M,
  )
  .into();

  let pump = Beam::new(
    PolarizationType::Extraordinary,
    0. * DEG,
    0. * DEG,
    775e-9 * M,
    100e-6 * M,
  )
  .into();

  let mut crystal_setup = CrystalSetup {
    crystal: CrystalType::KTP,
    pm_type: PMType::Type2_e_eo,
    phi: 0. * DEG,
    theta: 90. * DEG,
    length: 9000. * MICRO * M,
    temperature: 293.15 * K,
    counter_propagation: true,
  };

  let data: Vec<(_, _)> = (1..20)
    .map(|count| {
      let crystal_length = (1000 * count) as f64 * MICRO * M;
      crystal_setup.length = crystal_length;
      let pp =
        PeriodicPoling::try_new_optimum(&signal, &pump, &crystal_setup, Apodization::Off).unwrap();
      if let PeriodicPoling::On { period, sign, .. } = pp {
        (
          *(crystal_length / (MICRO * M)),
          *(sign * period / (MICRO * M)),
        )
      } else {
        (0., 0.)
      }
    })
    .collect();

  // print as csv
  data.iter()
  // also get difference
  .scan(0., |prev, (length, period)| {
    let diff = period - *prev;
    *prev = *period;
    Some((*length, *period, diff))
  })
  .for_each(|(length, period, diff)| {
    println!("{},{},{}", length, period, diff);
  });

  let pp =
    PeriodicPoling::try_new_optimum(&signal, &pump, &crystal_setup, Apodization::Off).unwrap();
  if let PeriodicPoling::On { period, sign, .. } = pp {
    let optimum_idler = IdlerBeam::try_new_optimum(&signal, &pump, &crystal_setup, pp).unwrap();
    println!("Crystal Setup: {:#?}", crystal_setup);
    println!("Signal: {:#?}", signal);
    println!("Optimum idler: {:#?}", optimum_idler);
    println!("Poling Period: {}nm", sign * period / (NANO * M));
  };
}
