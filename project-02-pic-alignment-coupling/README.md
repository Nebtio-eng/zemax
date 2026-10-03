# Alignment-Tolerant Fiber-to-PIC Coupling Using Micro-Optics

Master's project, Optical Technologies. System-level optical modelling of a
fiber-to-photonic-integrated-circuit coupling interface in Ansys Zemax
OpticStudio, with Lumerical FDTD used for the grating coupler.

**Status:** Zemax phase complete through Stage 11 (baseline, three
architectures, parameter study, Pareto front, consolidated code). Stage 8b
(Lumerical FDTD source) is next. Results: [`docs/results.md`](docs/results.md);
stage-by-stage record: [`docs/validation.md`](docs/validation.md).

---

## Research question

Can a micro-optical fiber-to-PIC coupling interface increase alignment tolerance
while maintaining acceptably low coupling loss?

The project does not assume the micro-lens approach succeeds. A result showing
no useful improvement would be a valid outcome.

## Motivation

Light leaves a photonic chip's grating coupler as a beam roughly 9 um across. A
fiber misaligned by about 2 um already loses a significant fraction of that
light. Hitting that tolerance requires active alignment — a machine that powers
the chip, searches for peak coupling and then bonds — which dominates the cost of
photonic packaging.

Placing a micro-lens in the path expands the beam, which should relax the
positional tolerance. The question is what that costs.

## Literature basis

**Primary reference:** N. Mangal, B. Snyder, J. Van Campenhout,
G. Van Steenberge and J. Missinne, "Monolithic integration of microlenses on the
backside of a silicon photonics chip for expanded beam coupling," *Optics
Express* **29**(5), 7601–7615 (2021). [doi:10.1364/OE.412353](https://doi.org/10.1364/OE.412353)

Selected as the only screened candidate reporting enough quantitative detail to
build the model without inventing parameters, and because its own design values
were produced in OpticStudio.

Full screening, including rejected candidates and the reasons, is in
[`literature/papers.md`](literature/papers.md). Every model parameter is listed
with its provenance (reported / derived / assumed / Zemax material data) in
[`literature/extracted_parameters.csv`](literature/extracted_parameters.csv).

## Physical model

```
PIC grating coupler
   |  beam expands through the Si substrate
micro-lens  (collimation)
   |  free-space propagation
optical fiber
```

Modelled in the PIC-to-fiber direction, with reciprocity invoked for the reverse
case. Coupling efficiency is computed as the Physical Optics Propagation
fiber-coupling integral, `eta = S x T` (system efficiency x mode overlap), with
loss `L = -10 log10(eta)`.

Method, model boundary and the resolved specification gaps are documented in
[`docs/methodology.md`](docs/methodology.md).

## Analytic ground truth

Closed-form Gaussian overlap results used to validate every simulated sweep, for
modes of 1/e^2 intensity radius `w`:

| Axis | 1-dB tolerance |
|---|---|
| Lateral | `d = 0.48 w` |
| Angular | `theta = 0.1528 lambda / (n w)` |
| Longitudinal | `z = 1.018 z_R`, `z_R = pi w^2 n / lambda` |

These reproduce the primary paper's three measured tolerances to within roughly
20%. The product `d x theta = 0.0733 lambda / n` is independent of beam size —
beam expansion trades positional tolerance against angular tolerance at a fixed
exchange rate.

## Validation standard

Results are validated on three independent legs: closed-form theory, the Zemax
POP model, and the primary paper's measurements. The report distinguishes
throughout between what was *experimentally demonstrated* in the literature and
what was *numerically determined* here.

## Current results (numerically determined, 1310 nm, 4.6 um source waist)

| | A0 top-side | A1 backside, flat | B backside + lens | B paper (measured) |
|---|---|---|---|---|
| Nominal loss | 0.165 dB | 6.26 dB | 0.0004 dB | - |
| Lateral 1-dB (rel / abs) | 2.25 / 2.06 um | 6.19 um / none | 8.14 / 8.14 um | +/-7 um |
| Angular 1-dB | 2.45 deg | 1.32 deg / none | 0.676 deg | +/-0.6 deg |
| Longitudinal 1-dB | 52 um | 249 um / none | ~700 um | 700 um |
| 2-D area A_1dB (abs) | 13.3 um^2 | 0 | 208 um^2 | - |

- The lens converts positional into angular tolerance at a fixed exchange rate:
  lateral x angular = 0.0960 um rad for every matched design. Longitudinal is
  the only net gain.
- The substrate supplies the beam width; the lens makes it usable by
  flattening the phase. Without it (A1) there is no 1-dB window at all.
- A spherical lens is optically sufficient; the conic constant has no effect.
  The lens ROC is the critical manufacturing tolerance (1 dB at -21.7% / +40.2%).
- Passive alignment (+/-10 um lateral) is reached at >= 786 um of silicon with
  the ROC and fibre scaled to it, paid for in angular tolerance (0.55 deg),
  not in loss.
- Excluded from every number: grating directionality and Fresnel reflection
  (an uncoated Si-air face would cost 1.6 dB). See
  [`docs/limitations.md`](docs/limitations.md).

## Repository layout

```
literature/   paper screening and extracted parameters (with provenance)
zemax/        OpticStudio models + run_config.json per model
              baseline/ (A0)  microlens/ (A1, B, B592, B_incidence)  glasscat/ (SILICON_1310)
python/       analytic.py          independent Gaussian model (never imports the Zemax path)
              coupling_analysis.py OpticStudio session, verified POP settings, convergence test
              alignment_sweep.py   Stages 3/5/7/8 sweeps     tolerance_map.py  2-D maps
              parameter_study.py   Stage 9                   optimization.py   Stage 10
              run.py               single entry point        test_validated.py regression tests
results/      coupling_curves/ tolerance_maps/ parameter_studies/ optimization/
docs/         methodology, validation, results, limitations, brief-questions,
              zosapi_gotchas (every API pitfall met in this project)
```

## Environment

Ansys Zemax OpticStudio 2026 R1.00 (Enterprise), driven through ZOS-API with
`pythonnet`; Ansys Lumerical FDTD via `lumapi`. Windows.

## Reproducing

Requirements: Windows, OpticStudio with ZOS-API (Professional or above), 64-bit
Python 3 (pythonnet will not load the ZOS-API assemblies from 32-bit).

```
python -m venv .venv
.venv\Scripts\pip install -r python\requirements.txt
cd python
..\.venv\Scripts\python test_validated.py                 # analytic checks, no OpticStudio
set ZEMAX_MODE=standalone
..\.venv\Scripts\python test_validated.py                 # + null test, A0 0.962617754, B 8.140 um, invariant
..\.venv\Scripts\python run.py 7                          # any stage: 3 5 6 7 8 9 10, or all
```

Every stage reads its settings from a `run_config.json` (models, POP tokens,
surface settings, sweep ranges) and writes one beside its results. The binary
POP settings file is generated from that JSON and every value is read back and
asserted. Default mode is a headless standalone instance (never run two at
once); `--mode extension` drives an open OpticStudio with the Interactive
Extension armed. Runs append to `results/progress.log`.

## Limitations

Zemax POP uses scalar diffraction theory and does not model grating period, etch
depth, directionality, Bloch modes or silicon waveguide modes. The grating is
represented by its equivalent Gaussian output; electromagnetic accuracy at the
grating is addressed separately in Lumerical FDTD. No claim of PIC-level
electromagnetic accuracy is made from the Zemax model alone.

## References

1. Mangal et al., *Opt. Express* **29**(5), 7601 (2021). doi:10.1364/OE.412353
2. Mangal et al., *IEEE JSTQE* **26**(2), 1–7.
3. Gradkowski & O'Brien, *Appl. Opt.* **63**(32), 8407 (2024). doi:10.1364/AO.540682
4. Gradkowski & O'Brien, *Appl. Opt.* **64**(14), 4014 (2025). doi:10.1364/AO.557834
5. Scarcella et al., *IEEE Photon. Technol. Lett.* **29**(22), 1943 (2017). doi:10.1109/LPT.2017.2757082
6. Ansys Optics, "Integrated microlens and grating coupler for photonic integrated circuits," application gallery.