# Validation log

Each stage records what was checked, what was predicted *before* measuring, and
what the model returned. Values from our simulation are **numerically
determined**; nothing here is experimentally demonstrated.

---

## Stage 3 — Baseline A0 and null test (2026-10-03)

Model: `zemax/baseline/A0_conventional.zmx`. Settings and raw numbers:
`zemax/baseline/run_config.json`, `zemax/baseline/a0_baseline_popd.csv`.
OpticStudio 2026 R1.00, build 260127, Interactive Extension (Plugin mode).

### Model as built

| Surface | Role | Thickness | Material |
|---|---|---|---|
| 0 OBJ | unused, POP starts at surface 1 | infinity | - |
| 1 STO | Gaussian source plane (grating-coupler equivalent mode) | 0.020 mm (20 um gap) | air |
| 2 IMA | SMF facet, receiver mode defined here | - | - |

No refracting surface, so Resample After Refraction and AutoResample are False
on every surface. Lens units are millimetres (the OpticStudio enum has no micron
option), so 4.6 um is entered as 0.0046 mm.

POP: start surface 1, end surface 2, wavelength 1.31 um, Gaussian Waist beam,
waist X = Y = 0.0046 mm, grid 512 x 512, window 0.08 mm x 0.08 mm (point spacing
156 nm). Receiver: Gaussian fiber, waist X = Y = 0.0046 mm (= MFD 9.2 um / 2).

### Null test (the factor-of-two check) — PASSED

Prediction: identical modes at zero distance give eta = 1, S = T = 1.

| Run | Prediction | Result |
|---|---|---|
| Identical 4.6 um modes, zero distance | 1.000 | **1.000000** (S 1.000000, T 1.000000) |
| Control 1: receiver entered as diameter, 9.2 um | 4 w1^2 w2^2 / (w1^2+w2^2)^2 = 0.640 | **0.640000** |
| Control 2: receiver decentred by d = w = 4.6 um | exp(-1) = 0.367879 | **0.367879** |

A zero-distance null test alone is blind to a *common-mode* factor of two (if
POP read both entries as diameters, both would be wrong equally and eta would
still be 1). Control 2 closes that gap: it returns exp(-d^2/w^2), which only
holds if the entered number is the 1/e^2 intensity **radius**. Independent
confirmation: POP's pilot-beam Rayleigh range is 50.745 um, which equals
pi w^2 / lambda with w = 4.6 um as a radius (it would be 12.7 um if read as a
diameter). Convention confirmed: **Waist entry = MFD / 2**.

### A0 nominal result

Predicted before building (Gaussian overlap, z = 20 um, z_R = 50.745 um):
eta = 1/(1 + (z/2 z_R)^2) = 0.9626, loss 0.165 dB, beam radius at fiber 4.944 um.

| Quantity | Predicted | POP (POPD) |
|---|---|---|
| eta_total (POPD 0) | 0.9626 | **0.962617754** |
| S, system (POPD 1) | 1 | **1.000000000** |
| T, receiver (POPD 2) | 0.9626 | **0.962617754** |
| Loss | 0.165 dB | **0.1655 dB** |
| Beam radius at fiber face | 4.944 um | **4.9444 um** |

S = 1 and T = eta: the loss is entirely modal (the beam has diverged slightly
from the receiver mode), none of it geometric. Expected, since A0 has nothing
that can truncate the beam.

### Convergence — PASSED, with a caveat

| Grid | Window | eta | Change vs 512 |
|---|---|---|---|
| 512 | 80 um | 0.962618 | - |
| 1024 | 80 um | 0.962618 | 0.0000% |
| 2048 | 80 um | 0.962618 | 0.0000% |
| 1024 | 160 um | 0.962618 | 0.0000% |

Requirement was < 0.5%. **Caveat:** the result is identical to six digits because
A0 is an ideal Gaussian through free space with no aperture, so sampling cannot
distort it. This test is weakly discriminating for A0. It becomes meaningful for
A1 and B, where the lens surface and truncation can alias, and it must be rerun
there rather than assumed.

### Prop Report

The POP results object reports 0 messages and no warnings or aliasing text in the
report body. The ZOS-API does not expose the Prop Report tab itself, so the GUI
tab was not read programmatically. **To be confirmed by eye in the GUI.**

### Prediction for Stage 5 (recorded before any sweep)

Lateral 1-dB tolerance of A0 at the 20 um nominal gap, measured as the extra 1 dB
from the on-axis peak:

- `d_1dB = 0.339 * sqrt(w1^2 + w2^2)` with w1 = 4.944 um (beam at fiber face),
  w2 = 4.6 um: **2.29 um**.
- Including the arriving beam's wavefront curvature (R = 149 um), exact Gaussian
  overlap gives **2.25 um** (eta falls as exp(-d^2 / 21.98 um^2)).
- Expectation: a Zemax lateral sweep lands within about 2.2 to 2.3 um. This is a
  prediction to test, not a value to tune toward. A0 depends strongly on the gap
  (Rayleigh range only 51 um), so Stage 5 reports tolerance against gap, not a
  single number.

### ZOS-API facts learned this stage (Stage 3)

- Lens units cannot be micrometres (mm, cm, in, m only).
- POP settings are set with `ModifySettings(cfgFile, token, value)`. It returns
  True for **any** token, valid or not, so every setting must be read back from
  the report. Tokens (found in `ZemaxCore.dll`): `POP_START`, `POP_END`,
  `POP_BEAMTYPE` (0 = Gaussian Waist), `POP_PARAM1..8` (waist X/Y, decenter X/Y,
  aperture X/Y, order X/Y), `POP_SAMPX/Y` (index: 5 = 512, 6 = 1024, 7 = 2048),
  `POP_WIDEX/Y`, `POP_COMPUTE`, `POP_FIBERTYPE`, `POP_FPARAM1..8` (same layout
  for the receiver; FPARAM3 = X decenter).
- Apply a setting file with `LoadFrom(cfg)`; run with `ApplyAndWaitForCompletion`;
  read the report with `GetResults().GetTextFile(path)`.
- **POPD reads the saved default POP settings**, not the open analysis. It returns
  0 until `LoadFrom(cfg)` then `Save()` is called.
- POPD columns: Param1 = surface, Param2 = wavelength, Param3 = field,
  Param4 = Data (0 total, 1 system, 2 receiver).
- The Interactive Extension arming is consumed by the first connection that uses
  it. Keep one process connected for the whole session.

---

## Stage 4 — Validation of A0 against the literature (PARTIAL)

**Stage 4 validation of A0 is partial.** The primary paper (Mangal et al. 2021)
reports no conventional top-side measurement, only the backside expanded-beam
device, so there is no paper number to compare A0 against. A0 is therefore
checked against the analytic Gaussian-overlap relations (`methodology.md`
section 5) and the general literature figure of roughly +/-2 um for conventional
top-side coupling only. The comparison against the paper is deferred to Stage 8,
where A1 and B can be set against the paper's measured +/-7 um, +/-0.6 deg and
700 um.

---

## Stage 5 — A0 lateral (X) tolerance (2026-10-03)

Script: `python/alignment_sweep.py`. Outputs: `results/coupling_curves/`
(per-gap sweep CSVs, `a0_tolerance_vs_gap.csv`, `a0_tolerance_vs_gap.png`,
`run_config_a0_lateral.json`). All values **numerically determined**.

### Method

- Receiver X decenter (`POP_FPARAM3`) swept from -6 to +6 um: 0.05 um steps at
  the 20 um nominal gap, 0.1 um elsewhere. POPD 0, 1 and 2 logged at every point.
- 1-dB tolerance = the decenter at which loss rises 1 dB above its on-axis value,
  found by cubic-spline interpolation of loss(d) plus a root find, on each side
  separately. Not the nearest sample.
- The extractor is self-tested on analytic data first (recovers 2.249773 um).
- Thinning the nominal sweep from 0.05 to 0.2 um changes the tolerance by
  6.5e-7 um, so the step size is far finer than three significant figures need.

### Review fixes made before the sweep

1. **Read-back is now a permanent assertion.** The POP .CFG is binary, and
   `ModifySettings` returns True for any token. The script locates each
   setting's byte offset in the CFG (by writing two marker values and
   diffing), then after every write reads the value back and raises
   `SettingsMismatch` on any difference. The POP text report is also checked
   against the configuration: grid size, window, wavelength, end surface,
   source waist, pilot position = gap, and beam radius = analytic. The on-axis
   eta must match the analytic eta0 to 1e-4, which checks the receiver waist.
   926 read-backs passed in the full run. A deliberate mismatch was confirmed
   to raise.
2. **The CFG is generated from `run_config.json`** (`pop_tokens`,
   `surface_settings`). This matters because POPD does not read the analysis
   window. It reads the *saved default* POP settings, a binary file outside
   git. Without generating it from text, the committed `.zmx` would silently
   give POPD = 0, or results computed with whatever settings were saved last.
   The generated CFG reproduces the committed Stage 3 eta (0.962617754)
   exactly, and the script checks this before every sweep.
3. Stage 4 partial-validation note added above.

### Result at the nominal 20 um gap

| | Tolerance | Measured vs it |
|---|---|---|
| **Zemax POP, interpolated** | **2.2498 um** | — |
| Prediction, 0.339 sqrt(w1^2 + w2^2) | 2.2894 um | -1.73% |
| Prediction, exact overlap incl. wavefront curvature | 2.2498 um | 0.000% |

Both sides agree (symmetry eta(+d) = eta(-d) within 1e-6, enforced). Within the
10% threshold, so no investigation was triggered. The 0.339 relation is slightly
optimistic because it ignores the curvature of the arriving wavefront.

### Tolerance against gap

| Gap (um) | eta0 | L0 (dB) | d_1dB measured (um) | Exact | 0.339 relation | vs 0.339 | Absolute 1 dB (um) |
|---|---|---|---|---|---|---|---|
| 0 | 1.0000 | 0.000 | 2.2073 | 2.2073 | 2.2053 | +0.09% | 2.207 |
| 5 | 0.9976 | 0.011 | 2.2100 | 2.2100 | 2.2107 | -0.03% | 2.198 |
| 10 | 0.9904 | 0.042 | 2.2180 | 2.2180 | 2.2266 | -0.39% | 2.171 |
| 20 | 0.9626 | 0.165 | 2.2498 | 2.2498 | 2.2894 | -1.73% | 2.055 |
| 50 | 0.8047 | 0.944 | 2.4607 | 2.4607 | 2.6878 | -8.45% | 0.584 |
| 100 | 0.5074 | 2.947 | 3.0988 | 3.0988 | 3.7824 | **-18.07%** | none (L0 > 1 dB) |

"Absolute 1 dB" is the decenter at which total loss reaches 1.000 dB.

**The 100 um disagreement with the 0.339 relation (-18%) was investigated, not
adjusted.** It is explained in full by wavefront curvature. The 0.339 relation
assumes both modes have flat phase. By 100 um the beam is past its Rayleigh
range (51 um) and strongly curved (R = 126 um). Including the curvature term in
the overlap reproduces Zemax to better than 1e-4 um at every gap. The simple
relation is valid only for gaps well inside z_R.

**Interpretation.** Measured relative to its own peak, the A0 lateral tolerance
depends only *weakly* on gap: 2.21 to 3.10 um over 0 to 100 um. The strong gap
dependence is in the **loss**: eta0 falls from 1.00 to 0.51. A larger gap
"buys" tolerance only by first spending loss. On the absolute convention, the
1-dB window collapses from 2.2 um to 0.58 um at 50 um and does not exist at
100 um. Reporting only the relative number would make a large gap look
better than it is.

### Agreement level, and what it does and does not show

POP agrees with the exact Gaussian overlap to about 1e-12 in eta at every gap,
and grid doubling changes nothing (0.0000% at d = 0 and at d = tolerance, every
gap). This is consistent with a correct grid calculation: the trapezoid rule on
a well-sampled Gaussian converges exponentially. It confirms the *setup*
(waist convention, decenter, gap, wavelength) is what we think it is. It does
**not** exercise POP's numerical propagation in any demanding way, because A0
contains no aperture, lens or aberration. That test starts with A1 and B.

### POP propagator: a silent grid change found and fixed

At gaps beyond the Rayleigh range, the default POP propagator **rescales its own
grid**: the window went from the configured 0.08 mm to 0.67 mm at 80 um and
0.84 mm at 100 um (1.64 um pixels). The report assertion caught this and
stopped the first sweep. eta was unaffected (identical to 1e-12 with the
alternative), but the configured window and the convergence test no longer
described what POP was doing. **Fix:** surface 1 uses
`UseAngularSpectrumPropagator = True` (in `run_config.json` →
`surface_settings`, applied and read back by the script). The window then stays
at 0.08 mm (>= 4x the largest beam radius, 10.2 um) at every gap. This is not
tuning: the results are identical in both modes, and the change keeps the stated
sampling true. The Stage 3 run used the default propagator; at 20 um the two
give the same value.

### Prop Report

The API reports 0 POP messages at every gap. **By-eye check of the GUI Prop
Report tab: pending.** Claude has no view of the OpticStudio screen, so this
line is to be completed by the user.
