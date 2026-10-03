# Project 2 — Alignment-Tolerant Fiber-to-PIC Coupling

Master's project, Optical Technologies. Zemax OpticStudio + Lumerical FDTD via
Python. This file is the standing context for every session in this repo.

---

## Who I am

I am a Master's student. I am **new to photonics** and a
**beginner/intermediate Zemax user**. Teach, do not just deliver.

For every new Zemax or photonics operation:
1. What we are doing
2. Why it is necessary
3. The underlying physics, in plain language before any maths
4. Exactly where to click or what to run
5. What the important parameters mean
6. What result I should expect
7. How to verify it

**Stop and wait for my confirmation before moving to the next major stage.**
Do not run ahead. Do not dump a whole stage's work in one response.

---

## The research question

Can a micro-optical fiber-to-PIC coupling interface increase alignment tolerance
while keeping coupling loss acceptably low?

**Do not assume the micro-lens wins.** A negative or nuanced result is valid.

The expected answer, already derived (see `docs/methodology.md`): beam expansion
does not create tolerance, it *moves* it. Lateral tolerance improves linearly
with expansion factor, longitudinal improves quadratically, angular degrades
linearly, and the lateral-angular product is invariant at `0.0733*lambda/n`.
The project's job is to demonstrate this numerically and quantify the exchange
rate — not to rediscover it from scratch, and not to assume it without proof.

---

## Architecture (settled)

Grating coupler, vertical / backside configuration. Not edge coupling.

```
PIC grating coupler -> beam expands through Si substrate -> micro-lens
  -> free space -> fiber
```

Modelled in the **PIC -> fiber ("OUT") direction**, with reciprocity invoked for
the fiber -> PIC case. Reason: Ansys documents that POP cannot accurately cost
the fiber -> grating direction, because the receiver mode would have to be the
grating's own field.

---

## Literature

**Primary:** Mangal, Snyder, Van Campenhout, Van Steenberge & Missinne,
"Monolithic integration of microlenses on the backside of a silicon photonics
chip for expanded beam coupling," *Opt. Express* 29(5), 7601–7615 (2021).
DOI 10.1364/OE.412353.

Chosen because every parameter the model needs is published or derivable, and
because its own design table was generated in OpticStudio.

**Supporting:** Mangal et al. *IEEE JSTQE* 26(2) (companion numbers);
Gradkowski & O'Brien *Appl. Opt.* 63(32), 8407 (2024) (analytic tolerance
equations); Gradkowski & O'Brien *Appl. Opt.* 64(14), 4014 (2025) (lens design
equations); Scarcella et al. *IEEE PTL* 29(22), 1943 (2017) (benchmark:
1.7 dB, +/-30 um).

**Methodology template:** Ansys Optics application gallery, "Integrated
microlens and grating coupler for photonic integrated circuits."

Full screening in `literature/papers.md`. Parameters with provenance in
`literature/extracted_parameters.csv`.

---

## Environment

| | |
|---|---|
| OpticStudio | 2026 R1.00, build 260127, **EnterpriseEdition** |
| ZOS-API | valid, standalone (Server) mode confirmed working |
| Lumerical | licensed, FDTD available |
| Machine | Windows, 16 CPUs |

### ZOS-API facts verified on THIS install — do not re-derive

- `ZOSAPI_NetHelper.dll` lives at the **install root** in 2026 R1, not under
  `ZOS-API\Libraries\`. Probe both.
- **There is no `app.Edition` property.** `app.LicenseStatus` returns the
  `LicenseStatusType` enum and *is* the edition
  (`StandardEdition`/`ProfessionalEdition`/`PremiumEdition`/`EnterpriseEdition`/...).
- `ConnectAsExtension(0)` returns a stub with `LicenseStatus = Unknown` and
  `conn.IsAlive = False` unless OpticStudio is open with
  Programming > Interactive Extension armed.
- `CreateNewApplication()` works headless with no GUI. Mode reports `Server`.
- **Never run two standalone instances concurrently** — they collide over the
  licence seat and produce `FRU__delta_init()` errors that look like bugs.
  One instance, loop inside it. POP already uses multiple cores.

Use **Interactive Extension** for model-building stages so I can watch the GUI.
Use **standalone** for long sweeps.

---

## How coupling efficiency is computed

POP fiber-coupling integral. Not a detector-power metric.

```
eta_total = S (system efficiency) x T (receiver efficiency / mode overlap)
L = -10 * log10(eta)          # 10, not 20 — eta is a power ratio
```

`POPD` operands: `Data=0` total, `Data=1` system, `Data=2` receiver,
`10` waist, `23` radius, `26` M-squared.

**Log all of 0, 1 and 2 on every run.** The split between S and T tells us
whether a tolerance loss is geometric (beam walking off an aperture) or modal
(mode mismatch). Without it we have an observation, not a result.

### The factor-of-two trap

POP **Fiber Data** takes `Waist X`/`Waist Y` as the **1/e^2 intensity radius**.
Datasheets quote mode field **diameter**.

```
Waist entry = MFD / 2
```

Getting this wrong produces plausible, wrong numbers silently. The null test
(fiber to identical fiber, zero distance, `eta ~ 1`) exists to catch it.

---

## Analytic ground truth

Derived in `docs/methodology.md`. For modes of 1/e^2 intensity radius `w`:

```
lateral:       eta(d) = exp(-d^2 / w^2)          ->  d_1dB    = 0.48 * w
angular:       eta(th) = exp(-(k*th*w)^2 / 4)    ->  th_1dB   = 0.1528 * lambda / (n*w)
longitudinal:  eta(z) = 1/(1 + (z/(2*z_R))^2)    ->  z_1dB    = 1.018 * z_R
               with z_R = pi * w^2 * n / lambda
```

These reproduce all three of the primary paper's measured tolerances to ~20%.
**Every Zemax sweep gets compared against these before it is believed.**

---

## Model boundary — state this, never blur it

Zemax POP models: Gaussian launch, scalar diffraction, refraction and
aberration at the micro-lens, aperture truncation, rigid-body misalignment via
Coordinate Breaks, mode overlap.

Zemax does **not** model: grating period, etch depth, duty cycle, directionality,
Bloch modes, silicon waveguide modes, polarisation at sub-wavelength scale.

The scalar approximation is **justified, not apologised for**: the 32 um
expanded beam at 1310 nm diverges at roughly 1.5 degrees, far inside scalar
validity.

In all writing, distinguish:
- "experimentally demonstrated" — only when the literature directly supports it
- "numerically determined" — our simulation results
- interpretation vs. limitation, kept separate

---

## Stage plan

Phase A is Zemax only. Do not introduce FDTD before Stage 8.

| Stage | Content | Gate to pass |
|---|---|---|
| 1 | Understand the primary paper | I can explain the architecture back |
| 2 | Define the physical model, coordinate system | Optical path agreed before building |
| 3 | Baseline, no lens + **null test** (`eta ~ 1`) | Null test passes |
| 3a | Analytic model in Python (~20 lines) | Predicts paper's 3 tolerances |
| 4 | Validate baseline vs paper and vs analytics | Disagreements explained, not tuned away |
| 5 | X / Y / Z / angular sweeps, 1-dB tolerances | Agrees with Stage 3a |
| 6 | 2D X-Y tolerance map, area `A_1dB` | — |
| 7 | Micro-lens architecture | — |
| 8 | Compare A vs B, full table | — |
| 8b | **Lumerical FDTD source swap** (Phase B) | Only after Stage 8 is validated |
| 9 | Parameter study, one variable at a time | Physical expectation stated *first* |
| 10 | Robust optimisation, Pareto front | — |
| 11 | Python automation / ZOS-API consolidation | — |

### Non-negotiable rules

- **Never silently invent a parameter.** Every value is labelled REPORTED,
  DERIVED, ASSUMED, or ZEMAX_DATA in `extracted_parameters.csv`.
- **Never tune the model to force agreement** with the paper. Investigate the
  disagreement instead.
- **POP convergence test on every new configuration**: double the grid, require
  < 0.5% change in `eta`. Check the Prop Report tab for warnings.
- Nominal z is re-optimised before measuring 1-dB tolerance (the peak moves once
  the lens is inserted). Report both conventions.
- Keep the software simple. Sophistication comes from the physics and the
  validation, not from architecture.

---

## Repo layout

```
literature/   papers.md, extracted_parameters.csv
zemax/        baseline/ conventional/ microlens/ optimized/
python/       coupling_analysis.py alignment_sweep.py tolerance_map.py
              parameter_study.py optimization.py
data/         literature/ raw/ processed/
results/      coupling_curves/ tolerance_maps/ parameter_studies/ optimization/
docs/         methodology.md validation.md limitations.md
report/       project_report.pdf
```

Do not create files that are not necessary.

## Reproducibility

Every run writes `run_config.json` beside its results CSV, recording: wavelength,
OpticStudio version and build, POP grid sampling X/Y, analysis window width,
beam type and waist, receiver waist, Resample-After-Refraction state per surface,
all radii/thicknesses/materials, coordinate-break values, sweep range and step,
optimisation operands and bounds, and raw POPD 0/1/2.
---

## Stage 2 decisions (settled 2026-10-03)

### Three configurations, not two

The improvement is decomposed rather than lumped, so that the substrate's
contribution is separated from the lens's.

| | Configuration | Isolates |
|---|---|---|
| **A0** | Conventional top-side: fiber above the grating, short air gap, no substrate traversal | The industry baseline (the ~±2 um case everyone quotes) |
| **A1** | Backside through the Si substrate, flat back face, no lens | What free beam expansion through existing material contributes |
| **B** | Backside with the monolithic etched micro-lens | What the lens adds on top of A1 |

Headline comparison is A0 -> B. The decomposition A0 -> A1 -> B is a result the
source papers do not report.

Working hypothesis to be tested, not assumed: **most of the lateral tolerance
gain comes from the substrate (A0 -> A1), not the lens (A1 -> B)**, because the
substrate takes the beam from 4.6 um to ~16 um radius while the lens only
collimates. If true, the lens's real contribution is to the longitudinal and
angular axes.

### Air / adhesive gap: 20 um nominal, swept

Chip back face (or top face, for A0) to fiber. Tagged ASSUMED.

Sourced from published optical-adhesive bond-line data: 3-50 um is the
documented range for fiber terminations and high-performance lens assemblies.
20 um sits inside that range at neither extreme.

The gap is swept 0-100 um regardless, which *is* the Stage 5 Z-tolerance
analysis, so it costs nothing extra and demonstrates the choice did not
determine the result.

Quantitative justification that the gap is not load-bearing for A1/B: the beam
leaving the back face has a radius of ~16 um, giving a Rayleigh range in air of
~614 um. A 20 um gap is ~3% of that, changing the beam radius by under 0.1%.

For A0 the gap *is* load-bearing (Rayleigh range only ~51 um for a 4.6 um
beam, so a 50 um gap grows it ~40%). A0 baseline tolerance must therefore be
reported as a curve against gap, not as a single number.

### Receiver mode is part of the configuration

Coupling efficiency is maximised when the receiving fiber mode matches the
arriving beam, so each configuration is paired with the fiber a packaging
engineer would actually choose:

| Config | Arriving beam radius | Matched receiver |
|---|---|---|
| A0 | 4.6 um | standard SMF, 9.2 um MFD |
| A1 | ~16 um, diverging | TEC fiber, ~32 um MFD (ASSUMED) |
| B | ~16 um, collimated | TEC fiber, ~32 um MFD (ASSUMED) |

Consequence, from `d_1dB = 0.339*sqrt(w1^2 + w2^2)` and
`eta_0 = 4*w1^2*w2^2/(w1^2+w2^2)^2`:

- A1/B with matched TEC: d_1dB = 7.7 um, eta_0 = 1. Reproduces the paper's +/-7 um.
- A1/B with SMF retained: d_1dB = 5.6 um but eta_0 = 0.28, i.e. a 5.5 dB penalty.

So the published tolerance **requires** the matched large-mode fiber. Expanding
the beam without also expanding the receiver buys some tolerance at a severe
efficiency cost.

**Both receiver cases are run** — matched as the headline, mismatched as an
instructive control, because the mismatched case is the clearest single
demonstration of the loss-versus-tolerance trade.

---

## Stage 7 predictions, recorded BEFORE building (2026-10-03)

### Substrate thickness and lens ROC are a matched pair — use 630 um

Collimation requires the surface to cancel the Gaussian wavefront curvature
arriving at it, `R(z) = z(1 + (z_R/z)^2)`, with `z_R = 177.8 um` in silicon.
Condition: `(n_Si - 1)/|R_surf| = n_Si/R(z)`.

| Substrate | Beam diameter | Wavefront R | ROC needed for collimation |
|---|---|---|---|
| 592 um | 32.0 um | 645 um | 461 um |
| **630 um** | **33.9 um** | **680 um** | **486 um** |

The paper's 480 um ROC therefore pairs with the companion paper's **630 um**
substrate (1.2% from ideal), **not** with the 592 um that the 2021 paper's 32 um
beam implies (4.1% off). **Use 630 um substrate + 480 um ROC + 34 um TEC
receiver.** Mixing 592 um with 480 um ROC leaves the beam imperfectly
collimated and would look like a modelling fault.

### Predicted results

| | A0 | A1 (flat exit) | B (lens) |
|---|---|---|---|
| Beam radius at fiber | 4.94 um | 18.7 um | 16.93 um |
| Wavefront at fiber | nearly flat, R = 149 um | strongly curved, R = 214 um | flat |
| Nominal loss | 0.17 dB (measured) | **~4.0 dB** | **~0 dB** |
| Lateral 1-dB | 2.25 um (measured) | ~8 um | **8.13 um** |
| Angular 1-dB | - | - | **0.677 deg** |
| Longitudinal 1-dB | - | - | **700 um** |

Paper measured, for comparison with B: +/-7 um lateral, +/-0.6 deg angular,
700 um longitudinal.

Coupling formula used (validated: reproduces POP's A0 eta = 0.962618 to 6 dp):

```
eta = 4 / [ (w/wf + wf/w)^2 + (pi*w*wf/(lambda*R))^2 ]
```

with `w` the arriving beam radius, `R` its wavefront radius, `wf` the fiber mode
radius. The second term is the **phase** mismatch, and it is what kills A1.

### The Stage 2 hypothesis is REVISED

Stage 2 predicted that most of the lateral tolerance gain would come from the
substrate (A0 -> A1) rather than the lens (A1 -> B). The calculation above says
that framing is wrong.

A1 has nearly the same beam size as B and a similar lateral tolerance, but
**~4 dB of loss** even with an optimally chosen fiber, because the beam arrives
with a wavefront radius of ~214 um while a fiber mode has flat phase. Expanding
the beam without flattening its phase produces a wide beam no fiber can accept.

**Revised prediction, to be tested not assumed:**

> The substrate supplies the beam *width*; the lens makes that width *usable* by
> flattening the wavefront. Neither alone is a working interface. The lens's
> contribution is primarily **phase**, not size — which is why its benefit shows
> up in nominal loss and in longitudinal/working-distance tolerance rather than
> in lateral tolerance at close range.

Falsifiable consequences to check at Stage 7/8:
1. A1 nominal loss >> A0 and >> B (predicted ~4 dB).
2. A1 and B lateral 1-dB tolerances similar (~8 um) despite that loss gap.
3. B's advantage over A1 grows with gap, since A1's beam keeps diverging while
   B's does not.

---

## Stage 9 findings (2026-10-03) — read before Stage 10

### Sensitivity ranking, measured

| Rank | Parameter | Effect | Predicted? |
|---|---|---|---|
| 1 | Substrate thickness (ROC re-matched) | Sets beam width: lateral ∝ w, longitudinal ∝ w². 400→630 um gives lateral 5.44→8.13 um, working distance 294→681 um | yes |
| 2 | Receiver MFD | 20 um: 1.15 dB / 6.67 um. 34 um: 0.0004 dB / 8.14 um. 48 um: 0.52 dB / 9.97 um | yes, to 2 dp |
| 3 | Lens ROC alone | 1 dB at about **-21% / +42%** — **asymmetric** | magnitude yes, asymmetry NO |
| 4 | Incidence angle | 0.033 dB at 6 deg off normal. Negligible | yes |
| 5 | Gap | 0.011 dB at 100 um, 0.065 at 200 um | yes |
| 6 | Aperture | Flat above ~80 um diameter (**2.4x** beam radius, not the 3x rule of thumb); 0.11 dB at 50 um | threshold yes, value too conservative |
| 7 | Conic constant | **Nothing.** -10 to +10 moves loss 0.0002→0.0019 dB, lateral 8.142→8.135 um | yes |
| 8 | Lens material | Not a free parameter (monolithic Si) | - |

### Three findings that matter for Stage 10

**1. A deliberately over-strong lens buys working distance.** Going *below* the
collimating ROC trades loss for working distance; going *above* it is strictly
worse on both counts. One-sided trade.

| ROC | Loss | Working distance |
|---|---|---|
| 380 um | 0.91 dB | 966 um |
| **420 um** | **0.30 dB** | **965 um** |
| 450 um | 0.071 dB | 862 um |
| 480 um (collimating) | 0.0004 dB | 713 um |
| 510 um | 0.042 dB | 578 um |
| 600 um | 0.51 dB | 366 um |

**This is the first real Pareto point and it must appear on the Stage 10 front.**
If the optimiser does not find the ROC < 480 branch, it has not searched properly.

**2. The optimisation is two-dimensional, not eight.** From the ranking above,
conic, aperture (above the cliff), gap, incidence angle and material are all
flat. The only real freedoms are:

- **Beam size** (substrate thickness + matched receiver MFD) — chooses *where on
  the lateral/angular invariant curve* you sit. Cannot change the product.
- **ROC relative to collimation** — chooses how much *loss* you pay for *working
  distance*.

Stage 10 should optimise over those two axes and state explicitly that the others
were shown flat. A high-dimensional search would be dishonest about Stage 9.

**3. The invariant product is a mismatch meter, confirmed across all 8 sweeps.**
Flat at 0.09601 whenever the configuration stays matched (conic: all 7 points;
gap: 0 to 200 um). Rises whenever it does not (receiver 20 um: 0.1096; ROC
300 um: 0.1172). **Never falls below 0.0960.** Use it as a diagnostic in Stage 10:
any Pareto point whose product exceeds 0.0960 is mismatched, and the amount of
excess quantifies by how much.

### The paper discrepancy is now closed as far as this model can close it

A 24 um receiver MFD gives lateral 7.04 um — matching the paper's +/-7 um
exactly — but its angular tolerance rises to 0.83 deg, now **+38%** off the
paper's 0.6 deg (worse than B's +12.6%).

So **two independent attempts** to close the lateral gap (B592's smaller beam,
and a smaller receiver) both fail the same way: lateral improves, angular
worsens. That is the invariant. Combined with the paper's own pair sitting 24%
below the theoretical floor, the remaining discrepancy is attributable to the
experiment's definitions or measurement geometry, **not** to anything this model
omits. Do not keep hunting for a parameter that closes both.
