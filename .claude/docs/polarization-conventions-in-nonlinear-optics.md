# Polarization Conventions in Nonlinear Optics

Reference document covering the standard conventions for labeling polarization
eigenmodes in nonlinear optical crystals, with emphasis on SPDC source design.
This document captures research into textbook definitions, community standards,
and how they relate to the spdcalc codebase.

## Uniaxial Crystals: Ordinary and Extraordinary

Uniaxial crystals (e.g., BBO, LiNbO3) have a single optic axis (conventionally
the Z-axis) and two principal refractive indices:

- **Ordinary (o)**: polarization perpendicular to the optic axis; refractive
  index `n_o` is independent of propagation direction.
- **Extraordinary (e)**: polarization in the plane containing the optic axis and
  the propagation direction; refractive index `n_e(theta)` depends on the angle
  theta between the propagation direction and the optic axis.

Two sub-types based on index ordering:

| Type | Index relation | Example crystals |
|---|---|---|
| **Negative uniaxial** | n_o > n_e | BBO, LiNbO3, KDP |
| **Positive uniaxial** | n_o < n_e | Quartz |

For negative uniaxial: ordinary = slow (high n), extraordinary = fast (low n).
For positive uniaxial: ordinary = fast (low n), extraordinary = slow (high n).

## Biaxial Crystals: Fast and Slow

Biaxial crystals (e.g., KTP, KTA, LBO, BiBO) have three distinct principal
refractive indices along three orthogonal crystal axes:

```
n_x < n_y < n_z    (standard convention, IEEE/ANSI Std. 176)
```

For biaxial crystals, the terms "ordinary" and "extraordinary" **do not have
rigorous physical meaning**. Instead, the two polarization eigenmodes for a
given propagation direction are labeled:

- **Fast (f)**: the eigenmode with the lower refractive index (higher phase
  velocity)
- **Slow (s)**: the eigenmode with the higher refractive index (lower phase
  velocity)

This convention is used by:
- Dmitriev, Gurzadyan, Nikogosyan: *Handbook of Nonlinear Optical Crystals*
  (Springer)
- Smith, A.V.: *Crystal Nonlinear Optics with SNLO Examples* (uses "lo"/"hi")
- Boeuf, Migdall et al., NIST phasematching code (Optical Engineering, 2000)
- Roberts, D.A.: "Simplified characterization of uniaxial and biaxial nonlinear
  optical crystals" (IEEE J. Quantum Electron. 28, 2057, 1992) — the closest
  document to a formal standard

Many SPDC papers still informally use "o" and "e" for biaxial crystals by
analogy, but this is a **convention of convenience**, not a physically precise
description.

## Phase-Matching Type Classification

### General Definition

| Type | Definition (for downconversion: pump -> signal + idler) |
|---|---|
| **Type 0** | All three waves have the same polarization |
| **Type I** | Signal and idler have the same polarization, different from pump |
| **Type II** | Signal and idler have orthogonal polarizations |

### Notation for Uniaxial Crystals

Uses `o` and `e`:

- Type 0: `o->oo` or `e->ee`
- Type I: `e->oo` (negative uniaxial) or `o->ee` (positive uniaxial)
- Type II: `e->eo`, `e->oe` (negative uniaxial) or `o->oe`, `o->eo` (positive
  uniaxial)

### Notation for Biaxial Crystals

Uses `s` (slow) and `f` (fast):

- Type 0: `s->ss` (all slow) or `f->ff` (all fast)
- Type I: `f->ss` (pump fast, daughters slow) or `s->ff` (pump slow, daughters
  fast)
- Type II: `f->fs`, `f->sf` (pump fast) or `s->sf`, `s->fs` (pump slow)

### Type 0 and Quasi-Phase-Matching

Type 0 phase matching (all same polarization) **cannot be achieved with
birefringent phase matching** — it requires quasi-phase-matching (QPM) via
periodic poling. The advantage is access to the largest nonlinear coefficient
(e.g., d_33 in both LiNbO3 and KTP).

## QPM in Periodically Poled Crystals

For quasi-phase-matching with periodic poling, propagation is along a principal
crystal axis. This simplifies the eigenmodes to pure principal-axis
polarizations.

### ppKTP (propagation along crystal X-axis)

For KTP with theta=90deg, phi=0deg, propagation is along the crystal X-axis.
The two eigenmodes are:

| Eigenmode | Polarization direction | Index | Wave type |
|---|---|---|---|
| Y-polarized | Crystal Y-axis | n_y (middle) | **Fast** |
| Z-polarized | Crystal Z-axis | n_z (highest) | **Slow** |

The physically allowed interactions (determined by KTP's mm2 point group
symmetry and the d-tensor) for collinear propagation along X are:

| d coefficient | Value | Interaction | PM Type |
|---|---|---|---|
| d_33 | ~10.7 pm/V | Z->ZZ (all slow) | Type 0: `s->ss` |
| d_24 | ~3.64 pm/V | Y->YZ (fast->fast+slow) | Type II: `f->fs` |
| d_32 | ~2.4 pm/V | Z->YY (slow->fast+fast) | Type I: `s->ff` |

The d_24 interaction (Type II) is the most commonly used for polarization-
entangled photon pair generation in ppKTP.

### ppLN (propagation perpendicular to crystal Z-axis)

For LiNbO3 (negative uniaxial, n_o > n_e), the standard QPM configuration
propagates perpendicular to the optic axis (Z). All beams polarized along Z
access d_33 (~27 pm/V), giving Type 0: `e->ee`.

## KTP Crystal Specifics

KTP (KTiOPO4) is a **positive biaxial** crystal with orthorhombic mm2 point
symmetry. Crystallographic axes a, b, c correspond to optical axes X, Y, Z:

```
n_x(1550nm) ~ 1.740    (a-axis)
n_y(1550nm) ~ 1.747    (b-axis)
n_z(1550nm) ~ 1.830    (c-axis)
```

For Type II SPDC in ppKTP (the standard entangled photon source configuration):
- Propagation along X (crystal a-axis)
- Pump: Y-polarized (fast, n_y) — labeled H in lab frame
- Signal: Y-polarized (fast, n_y) — labeled H
- Idler: Z-polarized (slow, n_z) — labeled V
- Nonlinear coefficient: d_24
- The birefringence (n_z - n_y) creates temporal walk-off between H and V
  photons

## Crystal Theta Convention in spdcalc

Setting `theta_deg: 90` with `phi_deg: 0` in spdcalc applies a rotation
`R_y(90deg)` to map between lab and crystal frames:

```
Lab Z (propagation) -> Crystal X
Lab Y               -> Crystal Y  (unchanged)
Lab X               -> Crystal -Z
```

This is the standard configuration for collinear QPM in ppKTP and ppLN: the
beam propagates along the crystal X-axis, and the crystal Z-axis (the poling
direction) lies perpendicular to propagation in the lab horizontal plane.

## Formal Standards

There is no single unified IEEE/OSA/SPIE standard for polarization labeling in
nonlinear optics. The closest references are:

1. **IEEE/ANSI Std. 176-1987** (Standard on Piezoelectricity): Defines
   crystallographic axis conventions. Widely adopted for crystal orientation.

2. **Roberts (1992)**, IEEE J. Quantum Electron. 28, 2057: "A plea for
   standardization of nomenclature and conventions." Proposes the Positive
   Nonlinear Optics Frame, recommends IEEE Std. 176 with a modification for
   mm2 crystals.

3. **Hobden (1967)**, J. Appl. Phys. 38, 4365: Original classification of
   phase-matched SHG directions in biaxial crystals (up to 72 classes).

The de facto community convention:
- Uniaxial crystals: `o`/`e`
- Biaxial crystals: `f`/`s` (or equivalently `lo`/`hi`)
- Index ordering: `n_x < n_y < n_z` for biaxial
- Type I/II defined by whether the two lower-frequency waves share the same or
  different polarization

## Key References

- Boyd, R.W. *Nonlinear Optics* (4th ed., Academic Press) — standard textbook;
  primarily covers uniaxial crystals
- Dmitriev, V.G., Gurzadyan, G.G., Nikogosyan, D.N. *Handbook of Nonlinear
  Optical Crystals* (Springer) — uses f/s for biaxial; comprehensive crystal
  data
- Smith, A.V. *Crystal Nonlinear Optics with SNLO Examples* — uses hi/lo index
  surfaces
- Roberts, D.A. IEEE J. Quantum Electron. 28, 2057 (1992) — standardization
  proposal
- Boeuf et al. Opt. Eng. 39, 1016 (2000) — NIST phasematching code; uses f/s
  for biaxial
- Couteau, C. "Spontaneous parametric down-conversion" (arXiv:1809.00127) —
  SPDC review
