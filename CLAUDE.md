# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SPDCalc is a Rust library for designing and analyzing spontaneous parametric downconversion (SPDC) sources used in quantum optics. SPDC generates pairs of entangled photons and is crucial for quantum information science. The library calculates phasematching properties, joint spectrum analysis (JSA/JSI), spectral purity, Hong-Ou-Mandel interference, fiber coupling efficiency, and photon count rates.

This is the core Rust library (published to crates.io) that powers both a web application at [app.spdcalc.org](https://app.spdcalc.org) and a Python wrapper [spdcalc-py](https://pypi.org/project/spdcalc-py/).

## Commands

### Build and Development
```bash
cargo build                      # Build the library
cargo build --release            # Build with optimizations (opt-level=3)
cargo test                       # Run all unit tests
cargo test --verbose             # Run tests with detailed output
cargo test <test_name>           # Run a specific test
cargo bench                      # Run performance benchmarks
cargo doc --open                 # Generate and view documentation
```

### Code Quality
```bash
cargo clippy --all --all-targets -- -D warnings   # Lint (must pass CI)
cargo fmt --all -- --check                        # Check formatting (must pass CI)
cargo fmt --all                                   # Auto-format code
```

### Feature Flags
```bash
cargo build --features wasm-bindgen   # Build with WebAssembly support
cargo build --features pyo3           # Build with Python bindings
cargo build --features elliptic       # Enable elliptic integral support
```

### Examples and Benchmarks
```bash
cargo run --example hom                  # Run Hong-Ou-Mandel example
cargo run --example counter_propagation # Run counter-propagation example
cargo bench --bench crystals            # Benchmark crystal calculations
cargo bench --bench heralding           # Benchmark heralding efficiency
cargo bench --bench hom                 # Benchmark HOM calculations
```

### CI Requirements
All PRs must pass:
- Clippy linting (Rust 1.78.0, no warnings)
- rustfmt formatting check
- Tests on stable, beta, and nightly toolchains

## Architecture

### Core Design Pattern

The library uses a **configuration-driven API** centered around the `SPDCConfig` struct:

1. **Configuration Phase**: Create or deserialize `SPDCConfig` from JSON/YAML using serde
2. **Validation Phase**: Call `config.try_as_spdc()` to validate and construct the `SPDC` object
3. **Calculation Phase**: Use the `SPDC` object to compute physical properties

The `SPDCConfig` supports "auto" values that are calculated during conversion. For example, `"idler": "auto"` automatically computes the idler beam properties based on energy and momentum conservation, and `"waist_position_um": "auto"` optimizes the beam waist position.

### Type Safety Through Dimensioned Units

All physical quantities use type-safe units via the `dimensioned` crate. Never use raw floats for physical values:

```rust
use spdcalc::prelude::*;

// Correct:
let wavelength = 1550e-9 * M;
let angle = 45.0 * DEG;

// Incorrect:
let wavelength = 1550e-9;  // Unitless, will not compile
```

Common unit types defined in [src/types.rs](src/types.rs): `Angle`, `Wavelength`, `Frequency`, `Length`, `Power`, `Time`, `Temperature`, `Intensity`, `ElectricField`, `RefractiveIndex`.

### Main Component Hierarchy

```
SPDC (Main Object)
  ├── CrystalSetup (crystal properties & configuration)
  │   ├── CrystalType enum (15 implementations: BBO_1, KTP, LiNbO3_1, etc.)
  │   └── Sellmeier equations (refractive indices)
  ├── Beams (Pump, Signal, Idler)
  │   ├── Gaussian beam propagation
  │   └── Wavevector calculations
  ├── PeriodicPoling (quasi-phasematching)
  │   └── Apodization support
  └── JointSpectrum (calculation engine)
      ├── Phasematching functions
      └── JSA/JSI computation
```

### Key Module Organization

**[src/spdc/spdc_obj.rs](src/spdc/spdc_obj.rs)** (16KB): The main `SPDC` struct providing the public API. All calculations start here.

**[src/crystal/](src/crystal/)**: Crystal properties and physics
- Each crystal is implemented as a module (e.g., [bbo_1.rs](src/crystal/bbo_1.rs), [ktp.rs](src/crystal/ktp.rs)) containing:
  - A `META` constant with crystal metadata
  - A `get_indices(wavelength, temperature)` function returning refractive indices
- `CrystalType` enum in [crystal_type.rs](src/crystal/crystal_type.rs) dispatches to the correct implementation
- `CrystalSetup` in [crystal_setup.rs](src/crystal/crystal_setup.rs) ties together crystal type, angles (θ, φ), length, and temperature
- `PMType` enum defines phasematching type (Type 0/1/2 with polarization combinations like "e->eo", "e->oo")
- Custom crystal expressions supported via JSON/HJSON with variables `l` (wavelength in μm) and `T` (temperature offset from 20°C)

**[src/beam/](src/beam/)**: Beam optics and propagation
- Three beam types (Pump, Signal, Idler) using a wrapper pattern for type safety
- Refractive indices, group/phase velocities, wavevectors
- External-to-internal angle conversion via Snell's law in [beam_waist.rs](src/beam/beam_waist.rs)
- Spatial walk-off calculations

**[src/phasematch/](src/phasematch/)**: Phasematching function calculations
- [delta_k.rs](src/phasematch/delta_k.rs): Wavevector mismatch Δk(ωs, ωi)
- [coincidences.rs](src/phasematch/coincidences.rs), [singles.rs](src/phasematch/singles.rs): Different phasematching scenarios
- [normalization.rs](src/phasematch/normalization.rs): JSI normalization constants

**[src/jsa/](src/jsa/)**: Joint Spectrum Analysis
- `JointSpectrum` struct in [joint_spectrum.rs](src/jsa/joint_spectrum.rs) computes JSA/JSI over frequency/wavelength ranges
- Provides iterators for signal-idler parameter spaces in [si_iterator.rs](src/jsa/si_iterator.rs)
- Integrates phasematching functions with pump spectrum

**[src/math/](src/math/)**: Numerical methods
- [integration.rs](src/math/integration.rs): Multiple integration algorithms (Simpson, Gauss-Legendre, Gauss-Konrod, Clenshaw-Curtis)
- [differentiation.rs](src/math/differentiation.rs): Numerical derivatives
- [nelder_mead.rs](src/math/nelder_mead.rs): Optimization algorithms
- [schmidt.rs](src/math/schmidt.rs): Schmidt decomposition for entanglement quantification

**[src/spdc/config/](src/spdc/config/)**: Configuration and deserialization
- Flexible JSON/YAML input with "auto" keyword support
- Automatic calculation of dependent parameters
- Validation during `try_as_spdc()` conversion

### Parallel Computation

The library uses `rayon` for parallel computation in hot paths. When working on performance-sensitive code, consider whether parallelization makes sense (typically for loops over independent calculations).

### Integration Methods

The library provides multiple numerical integration strategies via the `Integrator` enum:
- `Simpson`: Fast, good for smooth functions (default)
- `Romberg`: Adaptive Richardson extrapolation
- `GaussLegendre`: High accuracy for smooth functions
- `GaussKonrod`: Adaptive with error estimation
- `ClenshawCurtis`: Efficient for periodic functions

Choose the integrator based on the smoothness and behavior of the integrand.

## Testing Strategy

**Unit Tests**: Inline with `#[cfg(test)]` modules. Found in 17+ source files. Use `float-cmp` crate for floating-point assertions:

```rust
use float_cmp::approx_eq;
assert!(approx_eq!(f64, result, expected, epsilon = 1e-10));
```

**Benchmarks**: Located in [benches/](benches/). Use Criterion framework for statistical analysis.

**Examples as Integration Tests**: The [examples/](examples/) directory contains runnable examples that serve as both documentation and integration tests.

**Legacy Tests**: [test/tests.js](test/tests.js) contains JavaScript tests from the older implementation, used for cross-validation. Not actively maintained.

## Working with Crystals

Adding a new crystal requires:

1. Create a new module file in [src/crystal/](src/crystal/) (e.g., `my_crystal.rs`)
2. Define a `META` constant of type `CrystalMeta` with crystal metadata:
   - `id`: String identifier (e.g., "MyXtal_1")
   - `name`: Human-readable name
   - `reference_url`: Source for Sellmeier equations
   - `axis_type`: `OpticAxisType` (Uniaxial/Biaxial, Positive/Negative)
   - `point_group`: Symmetry group
   - `transmission_range`: Optional wavelength range
   - `temperature_dependence_known`: Boolean flag
3. Implement a `get_indices(wavelength: Wavelength, temperature: Kelvin<f64>) -> Indices` function
   - Use wavelength in micrometers for Sellmeier equations: `wavelength / (MICRO * M)`
   - Return `Indices::new(na::Vector3::new(nx, ny, nz))` or `(no, no, ne)` for uniaxial
   - Apply temperature corrections relative to 20°C: `(temperature - from_celsius_to_kelvin(20.0)) / K`
4. Register in `CrystalType` enum in [src/crystal/crystal_type.rs](src/crystal/crystal_type.rs):
   - Add variant to enum
   - Add case to `from_string()` match
   - Add to `get_all_meta()` vector
   - Add cases to `get_indices()` and `get_meta()` match statements
5. Add to module exports in [src/crystal/mod.rs](src/crystal/mod.rs)
6. Add tests comparing against published refractive index data

Reference existing crystals like [src/crystal/bbo_1.rs](src/crystal/bbo_1.rs) or [src/crystal/ktp.rs](src/crystal/ktp.rs) for the pattern.

## Configuration File Format

Typical JSON configuration structure:

```json
{
  "crystal": {
    "kind": "KTP",           // Crystal type (or custom expression)
    "pm_type": "e->eo",      // Phasematching: pump_pol->signal_pol idler_pol
    "phi_deg": 0,            // Crystal azimuthal angle
    "theta_deg": 90,         // Crystal polar angle
    "length_um": 14000,      // Crystal length
    "temperature_c": 20      // Operating temperature
  },
  "pump": {
    "wavelength_nm": 775,
    "waist_um": 200,
    "bandwidth_nm": 0.5,
    "average_power_mw": 300
  },
  "signal": {
    "wavelength_nm": 1550,
    "theta_external_deg": 0,
    "phi_deg": 0,
    "waist_um": 100,
    "waist_position_um": "auto"  // Auto-calculate optimal position
  },
  "idler": "auto",  // Auto-calculate from conservation laws
  "periodic_poling": {
    "poling_period_um": "auto"  // Auto-calculate from phasematching
  },
  "deff_pm_per_volt": 7.6  // Effective nonlinear coefficient
}
```

The `"auto"` keyword triggers automatic calculation during `config.try_as_spdc()`.

Custom crystal expressions can be used in place of named crystals:
```json
{
  "crystal": {
    "kind": {
      "no": "sqrt(2.7359+0.01878/(l^2-0.01822)-0.01354*l^2) - 9.3e-6 * T",
      "ne": "sqrt(2.3753+0.01224/(l^2-0.01667)-0.01516*l^2) - 16.6e-6 * T"
    }
  }
}
```
Variables: `l` = wavelength in μm, `T` = temperature offset from 20°C in Kelvin.

## Code Style

- Follow standard Rust conventions (enforced by rustfmt)
- Use descriptive variable names, especially for physics quantities
- Include doc comments (`///`) for public APIs
- Add inline comments for non-obvious physics or math
- Keep functions focused and modular
- Prefer type safety (dimensioned units) over raw numeric types
- Use 2 spaces for indentation size

## Key Dependencies

- **nalgebra**: Linear algebra (vectors, matrices, cross products)
- **ndarray**: N-dimensional arrays for JSA/JSI data
- **serde**: Serialization/deserialization of configs
- **dimensioned**: Type-safe dimensional analysis
- **rayon**: Parallel iterators
- **argmin**: Optimization algorithms
- **quad-rs**, **gauss-quad**, **quadrature**: Numerical integration
- **meval**: Mathematical expression parsing for custom crystals

## Excluded from Build

The [_scratchwork/](_scratchwork/) and [test/](test/) directories are excluded from the crate (see `Cargo.toml` exclude list). Development notes and legacy tests live here but are not shipped.

## Larger Ecosystem Context

This library is the backbone of two other spdcalc libraries:

- [spdcalc-ui](https://github.com/spdcalc/spdcalc-ui), the web-based SPA providing nice graphical tools for spdcalc including plots, calculations, importing and exporting data. It uses spdcalc as a wasm module.
- [spdcalc-py](https://github.com/spdcalc/spdcalc-py), the python bindings for spdcalc.

**Important:** Any modifications should take into account that those libraries depend on spdcalc and would break if the API changes, so breaking API changes should be called out.

## Other notes

