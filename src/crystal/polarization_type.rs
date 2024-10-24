use crate::SPDCError;
use std::fmt;
use std::str::FromStr;

/// The polarization type (ordinary or extraordinary)
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum PolarizationType {
  /// Ordinary polarization
  Ordinary,
  /// Extraordinary polarization
  Extraordinary,
}

impl fmt::Display for PolarizationType {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{:?}", self)
  }
}

impl FromStr for PolarizationType {
  type Err = SPDCError;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    let mut s = String::from(s);
    s.make_ascii_lowercase();
    match s.as_str() {
      "o" | "ordinary" => Ok(Self::Ordinary),
      "e" | "extraordinary" => Ok(Self::Extraordinary),
      _ => Err(SPDCError::new(
        "Can not parse polarization type".to_string(),
      )),
    }
  }
}
