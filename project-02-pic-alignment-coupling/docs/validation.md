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

### ZOS-API facts learned this stage

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
