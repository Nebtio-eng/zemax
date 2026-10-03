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