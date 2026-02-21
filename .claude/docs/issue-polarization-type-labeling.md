# Issue: PolarizationType and PMType Labeling

## Problem Summary

The `PolarizationType` enum uses "Ordinary" and "Extraordinary" labels, but the
Fresnel equation solver in `CrystalSetup::index_along()` actually implements a
**fast/slow** mapping that is only correct for negative uniaxial crystals. This
creates conceptual confusion, prevents correct handling of positive uniaxial
and custom expression crystals, and is missing PM types needed for certain
biaxial crystal interactions.

See [polarization-conventions-in-nonlinear-optics.md](polarization-conventions-in-nonlinear-optics.md)
for the physics background.

## Current Behavior

### The Fresnel Solver (`crystal_setup.rs:73-85`)

The `index_along()` method solves a quadratic (Eq. 11 from the NIST
phasematching paper) for the two refractive index eigenvalues along a given
propagation direction. The two roots always map as:

```
PolarizationType::Ordinary      -> high-n root (SLOW wave)
PolarizationType::Extraordinary -> low-n root  (FAST wave)
```

**The inline comments (lines 78-81) have fast/slow reversed.** They say
"fast" for Ordinary and "slow" for Extraordinary, but the math gives the
opposite. Verified numerically with BBO and KTP.

### Why 77 tests pass

All existing tests use either:
- Negative uniaxial crystals (BBO, LiNbO3) where ordinary = slow and
  extraordinary = fast — matching the code's convention
- Biaxial crystals (KTP) where there is no physical "ordinary"/"extraordinary"
  and the fast/slow mapping is internally consistent

The convention works because the PM type labels are interpreted consistently
with the solver: `e` always means "fast" and `o` always means "slow",
regardless of crystal type.

### Why the expression crystal test fails

The test `index_along_expr_crystal_test` creates a positive uniaxial crystal
with `no=1, ne=2` (indices `[1, 1, 2]`). It expects:

```
Extraordinary -> 2.0 (the ne value)
Ordinary      -> 1.0 (the no value)
```

But the code returns the opposite because it maps Extraordinary -> fast(low-n)
= 1.0 and Ordinary -> slow(high-n) = 2.0. The test was added in commit
78d297b as part of an attempted fix that was reverted in 85a2350 because
changing the solver broke all other tests.

### Missing PM Types

The enum only supports pump-fast configurations:

```rust
Type0_o_oo,   // pump slow, daughters slow
Type0_e_ee,   // pump fast, daughters fast
Type1_e_oo,   // pump fast, daughters slow
Type2_e_eo,   // pump fast, signal fast, idler slow
Type2_e_oe,   // pump fast, signal slow, idler fast
```

Missing pump-slow-with-mixed-daughters types:

```
Type1_o_ee    // pump slow, daughters fast    (needed for d_32 in ppKTP)
Type2_o_oe    // pump slow, signal slow, idler fast
Type2_o_eo    // pump slow, signal fast, idler slow
```

For ppKTP (propagation along X), the d_32 coefficient (~2.4 pm/V) enables
a Type I interaction: pump(Z/slow) -> signal(Y/fast) + idler(Y/fast). This
cannot be expressed with the current PMType enum.

## Mapping Between Conventions

For reference, here is how the current code labels translate to standard
conventions for common configurations:

### BBO (negative uniaxial, n_o > n_e)

| Code label | Physical meaning | Standard convention |
|---|---|---|
| `Ordinary` | n_o (slow, high-n) | Ordinary (correct) |
| `Extraordinary` | n_e (fast, low-n) | Extraordinary (correct) |
| `e->oo` (Type I) | pump(e/fast) -> daughters(o/slow) | Type I: `e->oo` (correct) |
| `e->eo` (Type II) | pump(e/fast) -> sig(e/fast) + idl(o/slow) | Type II: `e->eo` (correct) |

### KTP at theta=90 (biaxial, propagation along crystal X)

| Code label | Physical meaning | Standard biaxial convention |
|---|---|---|
| `Ordinary` | n_z = 1.830 (slow, Z-polarized) | Slow |
| `Extraordinary` | n_y = 1.747 (fast, Y-polarized) | Fast |
| `o->oo` (Type 0) | all Z/slow, uses d_33 | `s->ss` |
| `e->eo` (Type II) | pump(Y/fast) -> sig(Y/fast) + idl(Z/slow), uses d_24 | `f->fs` |
| N/A (missing) | pump(Z/slow) -> sig(Y/fast) + idl(Y/fast), uses d_32 | `s->ff` |

### Positive uniaxial crystal (e.g., custom no=1, ne=2)

| Code label | Physical meaning | Standard convention |
|---|---|---|
| `Ordinary` | ne = 2 (slow, high-n) | **Extraordinary** (WRONG) |
| `Extraordinary` | no = 1 (fast, low-n) | **Ordinary** (WRONG) |

This is the source of the failing test and conceptual confusion.

## Recommendations

### Phase 1: Immediate fixes (no behavioral changes, no API breakage)

1. **Fix the comments** in `index_along` (lines 78-81):
   ```rust
   // slow (higher refractive index)
   PolarizationType::Ordinary => -x2,
   // fast (lower refractive index)
   PolarizationType::Extraordinary => -x1,
   ```

2. **Fix the failing test** to match the actual code convention. The test
   expectations should be swapped, or (better) the test should be rewritten
   to demonstrate the actual mapping with a clear comment explaining why
   `Ordinary` gives the `ne` value for this positive uniaxial crystal.

3. **Add documentation** to `PolarizationType` and `PMType` explaining the
   actual semantics:
   ```rust
   /// The polarization eigenmode type.
   ///
   /// **Important**: In this codebase, `Ordinary` always selects the eigenmode
   /// with the *higher* refractive index (the "slow" wave), and `Extraordinary`
   /// always selects the eigenmode with the *lower* refractive index (the "fast"
   /// wave). This matches the standard physics convention for negative uniaxial
   /// crystals (like BBO), but is reversed for positive uniaxial crystals and
   /// does not correspond to the physical ordinary/extraordinary distinction
   /// for biaxial crystals.
   ///
   /// For biaxial crystals, think of `Ordinary` as "slow" and `Extraordinary`
   /// as "fast".
   ```

### Phase 2: Add fast/slow aliases (backward-compatible)

1. **Add string aliases** to `PolarizationType::from_str`:
   - `"s"` / `"slow"` -> `Ordinary` (slow)
   - `"f"` / `"fast"` -> `Extraordinary` (fast)

2. **Add string aliases** to `PMType::from_str`:
   - `"s->ss"` -> `Type0_o_oo`
   - `"f->ff"` -> `Type0_e_ee`
   - `"f->ss"` -> `Type1_e_oo`
   - `"f->fs"` -> `Type2_e_eo`
   - `"f->sf"` -> `Type2_e_oe`

3. **Add the missing PM types** to the enum:
   ```rust
   Type1_o_ee,   // s->ff: pump slow, daughters fast
   Type2_o_oe,   // s->sf: pump slow, signal slow, idler fast
   Type2_o_eo,   // s->fs: pump slow, signal fast, idler slow
   ```
   With corresponding string parsing for both `o`/`e` and `s`/`f` notations.

### Phase 3: Rename to fast/slow (breaking API change, major version)

1. **Rename** `PolarizationType` variants:
   ```rust
   pub enum PolarizationType {
     /// Eigenmode with lower refractive index (faster phase velocity).
     Fast,
     /// Eigenmode with higher refractive index (slower phase velocity).
     Slow,
   }
   ```

2. **Rename** `PMType` variants:
   ```rust
   pub enum PMType {
     Type0_s_ss,   // All slow
     Type0_f_ff,   // All fast
     Type1_f_ss,   // Pump fast, daughters slow
     Type1_s_ff,   // Pump slow, daughters fast
     Type2_f_fs,   // Pump fast, signal fast, idler slow
     Type2_f_sf,   // Pump fast, signal slow, idler fast
     Type2_s_sf,   // Pump slow, signal slow, idler fast
     Type2_s_fs,   // Pump slow, signal fast, idler slow
   }
   ```

3. **Keep backward compatibility** in parsing: `"e->eo"` still parses to
   `Type2_f_fs`, `"o->oo"` still parses to `Type0_s_ss`.

4. **Coordinate with downstream**: spdcalc-ui and spdcalc-py depend on the
   public API. The serialized form should remain backward-compatible (accept
   old format on input, emit new format on output).

### Optional Enhancements

- **Principal-axis notation for QPM**: Accept `"Y->YZ"` or `"Z->ZZ"` as PM
  type strings when theta=90 (propagation along a principal axis). This is
  natural for the ppKTP/ppLN community.

- **PM type validation**: Warn or error when a PM type has zero effective
  nonlinear coefficient for the given crystal symmetry and orientation (e.g.,
  `Type0_f_ff` / all-Y in KTP is forbidden because d_222 = 0 for mm2).

- **Display methods**: Add a method that returns the PM type in the
  crystal-axis notation (e.g., `"Y->YZ"`) given the crystal setup, for use in
  UI display and logging.

## Impact Assessment

### Phase 1 (immediate fixes)
- No behavioral changes
- No API breakage
- Fixes 1 failing test
- Improves developer understanding

### Phase 2 (aliases + missing types)
- Backward-compatible addition
- Enables new physics (d_32 interactions, positive uniaxial crystals with
  correct labeling via f/s notation)
- Minor additions to spdcalc-py to expose new types

### Phase 3 (rename)
- Breaking API change (major version bump)
- Requires updates to spdcalc-ui and spdcalc-py
- Significant improvement to conceptual clarity
- Aligns with established nonlinear optics conventions
