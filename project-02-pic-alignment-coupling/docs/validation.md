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

---

## Stage 6 — A0 X-Y tolerance map (machinery check) (2026-10-03)

Script: `python/tolerance_map.py`. Outputs: `results/tolerance_maps/`
(`A0_xy_loss_map.csv`, `A0_radial_cuts.csv`, `A0_xy_tolerance_map.png`,
`run_config_A0_map.json`). A0 at the 20 um gap.

Purpose: A0's source and receiver are both circular, so the 1-dB contour must
be a circle with the Stage 5 radius. Any ellipticity would mean the X and Y
decenters, the grid or the code treat the two axes differently. No new physics.

Token check first: `POP_FPARAM4` is the receiver Y decenter (a 4.6 um Y offset
at zero gap gives 0.367879441, which is e^-1). `POP_TILTX`/`POP_TILTY` are
receiver tilts in degrees. 2 deg gives 0.862084, against 0.862188 from the
small-angle formula; POP matches exactly if the tilt enters as tan(theta),
a difference of order theta^2/2, which is negligible below 2 deg.

| Check | Result |
|---|---|
| 1-dB radius along 24 directions (spline-interpolated cuts) | 2.249773 um in every direction |
| Stage 5 1-D tolerance | 2.249773 um |
| Circularity, r_max / r_min - 1 | 1.9e-13 (limit set at 5e-3) |
| A_1dB from radial cuts (integral of r^2/2) | 15.9011 um^2 |
| pi r^2 with r = 2.249773 um | 15.9011 um^2 (difference 2e-5 %) |
| A_1dB from the grid contour (0.25 um grid) | 15.839 um^2 (-0.39 %) |

**PASSED: the contour is a circle.** The grid-contour area is 0.4% low because
straight-line contour segments between grid points cut inside a curved
contour. The radial-cut value is the accurate one.

---

## Stage 7 — A1 (flat exit) and B (micro-lens) (2026-10-03)

Models: `zemax/microlens/A1/A1_flat_exit.zmx`, `zemax/microlens/B/B_microlens.zmx`,
built in code from their `run_config.json` (`coupling_analysis.build_model`).
Driver: `python/alignment_sweep.py` (`run_stage7`). Outputs:
`results/coupling_curves/` (A1_*, B_*, `stage7_A1_vs_B.png`,
`run_config_stage7.json`). All values **numerically determined**.

### Configuration

| Surface | A1 | B |
|---|---|---|
| 1 | Gaussian source, w0 = 4.6 um, 630 um of SILICON_1310 | same |
| 2 | Si exit face, **R = infinity**, 20 um air | Si exit face, **R = -480 um**, 20 um air |
| 3 | TEC fiber, MFD 34 um (waist entry 17 um) | same |

The code refuses to run unless the two configs differ in exactly one cell
(`surfaces[2].radius_mm`). Checked: that is the only difference.

**Why 630 um and 480 um are a matched pair.** After 630 um of silicon the
Gaussian wavefront has radius R = 680 um. A spherical Si-air surface collimates
when (n - 1)/|R_lens| = n/R_wavefront, i.e. R_lens = 486 um. 480 um is 1.2% from
that, so the beam leaves essentially flat (ABCD: wavefront radius -44 mm at the
fiber, waist 31 um past the exit face).

### Silicon index: an extrapolation, found and recorded

OpticStudio's SILICON (INFRARED.AGF, Salzberg & Villa 1957) is defined only for
1.36 to 11 um. At 1.31 um it is out of range, and OpticStudio **silently used
n = 1.0** (POP then failed). Fix: `zemax/glasscat/PROJECT02.AGF` holds a copy of
the same entry, renamed SILICON_1310, with only the lower wavelength limit moved
to 1.30 um. n(1.31 um) = 3.50391, matching the CSV's 3.504. This is a 50 nm
extrapolation of a smooth Sellmeier fit in silicon's transparent region, now
labelled as such in `extracted_parameters.csv`. `build_model` installs the file
into OpticStudio's glass directory.

### Lens sign — verified empirically

| Exit radius | Beam radius at fiber (POP pilot) | ABCD | Wavefront R | eta |
|---|---|---|---|---|
| **-480 um** | **16.920 um** | 16.920 um | -44 mm (flat) | **0.99992** |
| +480 um | 20.449 um | 20.453 um | 116 um (diverging) | 0.071 (11.5 dB) |

Negative radius collimates. In OpticStudio's convention a surface bulging toward
+z (toward the fiber) has its centre of curvature behind it, so R < 0, and its
power (n2 - n1)/R = (1 - 3.504)/(-480 um) is positive. The wrong sign makes the
beam diverge faster, as expected.

### Sampling, window and propagator

- **Window 400 um**, sized for the largest beam: B at a 2000 um gap has a
  ~52 um radius (window = 7.7 x radius, rule >= 4x). The decentred receiver stays
  inside too (24 um offset + 3 x 17 um = 75 um, against a 200 um half-width).
- **Grid 1024, pixel 0.39 um**, sized for the smallest beam: the 4.6 um source
  waist gets about 12 pixels per radius.
- **Propagators agree.** Angular-spectrum vs default: eta identical to 1e-13
  (A1) and 2.4e-7 (B). The default again rescales its grid (A1: 0.40 -> 0.67 mm;
  B: 0.40 -> 0.60 mm). Angular-spectrum keeps the configured window (B: 0.39999
  mm, a 0.0025% adjustment after the curved surface; the report check now allows
  0.1%). The 630 um of silicon does not change the Stage 5 choice.
- **Resample After Refraction stays OFF.** With RAR on and Auto off, POP
  re-grids onto the surface's default 32 x 32 grid over 1 mm and the result is
  meaningless (A1: eta 0.9996 instead of 0.2369). With Auto on it agrees to
  0.03%. The methodology's plan to set RAR on the lens would only be safe with
  Auto on. Nothing here needs it, so it is off and recorded.

### Convergence — passes, and the test can fail

| Point | A1, grid x2 / grid+window x2 | B, grid x2 / grid+window x2 |
|---|---|---|
| Nominal | 0.0000% / 0.0000% | 0.0000% / 0.0000% |
| Lateral 8 um | 0.0000% / 0.0000% | +0.0004% / -0.0010% |
| Tilt 0.6 deg | 0.0000% / 0.0000% | -0.0003% / +0.0008% |
| Largest beam (A1 300 um gap, B 2000 um gap) | 0.0000% / 0.0000% | +0.0021% / -0.0061% |

All far below 0.5%. To show this is not a test that always passes, B at nominal
was run with deliberately bad settings: a 64-point grid gives eta 0.939 (-6%);
a 40 um window gives 0.963 (-3.7%); a 60 um window gives -0.08%; a 128 grid
gives -0.001%. The chosen settings sit well inside the converged region.

POP messages: 0 for both models. 801 (A1) and 996 (B) settings read-backs passed.
POP pilot radius matches the ABCD prediction to < 0.003% in both models.

### Results against the predictions in CLAUDE.md

"rel" = 1 dB above the value at the reference point; "abs" = total loss reaches
1.000 dB. Both are reported because A1 and B have very different nominal losses.

| Quantity | Predicted (CLAUDE.md) | Zemax | vs prediction | Paper | vs paper |
|---|---|---|---|---|---|
| A1 nominal loss | ~4.0 dB | **6.255 dB** | **+56%, investigated** | - | - |
| B nominal loss | ~0 dB | **0.0004 dB** (eta 0.999918) | agrees | - | - |
| B lateral 1-dB | 8.13 um | **8.140 um** (abs 8.138) | +0.1% | +/-7 um | **+16%** |
| B angular 1-dB | 0.677 deg | **0.6758 deg** (abs 0.6757) | -0.2% | +/-0.6 deg | **+13%** |
| B longitudinal 1-dB, from 20 um nominal | 700 um | **713.3 um** (abs 713.2) | +1.9% | 700 um | +1.9% |
| B longitudinal 1-dB, from re-optimised gap (30.9 um) | 700 um | **702.3 um** | +0.3% | 700 um | +0.3% |

| A1 | rel | abs |
|---|---|---|
| Lateral 1-dB | 6.188 um | **none**: nominal is already 6.26 dB |
| Angular 1-dB | 1.318 deg | none |
| Longitudinal 1-dB (gap increasing) | 240 um | none |
| Best gap | 0 um (6.20 dB) | |

Every Zemax number agrees with the independent ABCD + numerical-overlap
calculation to within 0.2%.

**A1 loss: 6.26 dB, not 4 dB — investigated, not adjusted.** The ~4 dB in
CLAUDE.md is A1's loss with the *best possible* receiver for its beam, a 9.3 um
mode radius (18.6 um MFD). Zemax gives exactly that: **4.036 dB** (analytic
4.033). With the 34 um TEC receiver specified for both A1 and B it is 6.26 dB,
and with standard SMF 6.88 dB. The prediction was right for a different fiber.
Physics: A1's wavefront arrives curved (R = 213 um) while every fiber mode is
flat. A smaller mode spans less of the curved wavefront, so the best fiber is
*smaller* than the beam. The loss is entirely in T (receiver/modal); S = 1.000
in both models, so none of it is geometric.

**B vs the paper: lateral +16% and angular +13% — investigated, not adjusted.**
The model matches its own analytic prediction to 0.2%, so the gap is between
the idealised model and the experiment, not a setup error. Candidate
explanations, none tested yet:

1. **Receiver MFD is ASSUMED** (34 um, matched). The paper does not report it.
   A smaller real TEC mode would narrow the lateral tolerance. For example, a
   ~24 um MFD gives about 7.0 um lateral at a 0.55 dB mismatch penalty. This is
   a sensitivity to check (Stage 9), not a value to adopt.
2. **The real grating beam is elliptical** (`limitations.md` section 1), so the
   measured tolerance is axis-dependent. A circular model cannot reproduce that.
   Stage 8b tests it directly.
3. **Angular measurement pivot.** If the experiment rotated the fiber about a
   point behind its facet, every tilt also adds a lateral offset, which lowers
   the measured angular tolerance. The model tilts about the facet centre.
4. **Real lens imperfections** (sag error, roughness), absent from the model.

Longitudinal agrees with the paper to 0.3-1.9%.

**Invariant check.** For B, lateral x angular = 8.140 um x 0.011795 rad =
0.0960 um rad. The predicted invariant 0.0733 lambda/n = 0.0960 um rad (air).
Confirmed numerically. For A1 the product is 0.142: the invariant holds only
for flat-phase, matched modes, which A1 is not.

### The revised Stage 2 hypothesis, tested

1. *A1 loss >> A0 and >> B:* **supported.** 6.26 dB vs 0.17 dB vs 0.0004 dB.
2. *A1 and B lateral similar (~8 um):* **not supported.** Relative: A1 6.19 um vs
   B 8.14 um. Absolute: A1 has no 1-dB window at all. The curvature that costs
   A1 its loss also narrows its lateral tolerance.
3. *B's advantage grows with gap:* **supported.** A1 goes from 6.20 dB at 0 um to
   7.46 dB at 300 um. B stays within 1 dB out to 713 um.

Overall: the substrate supplies the beam *width* (A1 and B reach 18.7 and 16.9
um); the lens makes it *usable* by flattening the phase. On these numbers that
is worth 6.25 dB of loss and turns a configuration with no 1-dB window into one
with +/-8.1 um, +/-0.68 deg and 700 um. The relative-only view would have hidden
most of that: it rates A1's angular tolerance (1.32 deg) *better* than B's
(0.68 deg), because a coupling that is already losing 76% of the light changes
little, proportionally, when tilted.

**Mismatched control (B with SMF, 9.2 um MFD):** 5.92 dB nominal, lateral 5.95 um
relative and no absolute window. Expanding the beam without expanding the
receiver buys little and costs 6 dB, consistent with the Stage 2 decision.

---

## Stage 8 — data gaps filled (2026-10-03)

Driver: `python/alignment_sweep.py` (`run_stage8`); maps: `python/tolerance_map.py`.
Outputs: `results/coupling_curves/A0_angular.csv`, `A0_longitudinal.csv`,
`B592_*.csv`, `run_config_stage8.json`; `results/tolerance_maps/A1_*`, `B_*`.
Synthesis, comparison table and headline claim: **`docs/results.md`**.

| Quantity | Predicted (prompt) | Predicted (ABCD + overlap, before running) | Zemax |
|---|---|---|---|
| A0 angular 1-dB, 20 um gap | 2.32 deg | 2.448 deg (abs 2.236) | **2.446 deg** (abs 2.234) |
| A0 longitudinal 1-dB, from best gap (0 um) | 52 um | 51.6 um | **51.6 um** |
| A0 longitudinal 1-dB, from 20 um nominal | - | rel 36.3 / abs 31.6 um | **36.3 / 31.6 um** |
| A0 lateral x angular | 0.0911 um rad | 0.0961 | **0.0960** |
| B area A_1dB | ~208 um^2 | 208.2 | **208.2 um^2** (abs 208.1) |
| A1 area A_1dB | - | 120.3 rel, 0 abs | **120.3 um^2** rel, **0** abs |
| B592 lateral 1-dB | ~7.7 um | 7.677 um | **7.678 um** (abs 7.674) |
| B592 angular 1-dB | - | 0.717 deg | **0.7165 deg** |

The B and A1 maps are circular (circularity 1.5e-9 and 2.6e-13 over 16
directions), and X and Y tolerances are identical. Convergence: A0 grid
doubling at 2.4 deg tilt and at a 150 um gap changes eta by 3e-11%. B592 meets
the < 0.5% rule at nominal, at 8 um lateral and at 0.6 deg. B592's POP pilot
radius (15.999 um) matches ABCD to < 0.01%.

**A0 invariant prediction (0.0911, 5% below ideal) did not hold.** A0 sits on the
ideal 0.0960 um rad, as B does. The ideal is a minimum for matched flat modes;
mismatch raises the product (A1: 0.142), so curvature cannot pull A0 below it.
Explained in `results.md` section 2.

**B592 (2021-paper pair: 592 um Si, R = -461 um, 32 um MFD receiver)** is a
separate config, not a one-cell variant of B. Lateral lands at 7.68 um (+9.7%
against the paper's +/-7 um, against B's +16.3%), so 6.6 percentage points of
the lateral discrepancy come from the literature pair. Angular moves the other
way (+12.6% to +19.4%).
