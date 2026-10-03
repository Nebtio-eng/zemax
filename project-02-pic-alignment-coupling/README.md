# Alignment-Tolerant Fiber-to-PIC Coupling Using Micro-Optics

Master's project, Optical Technologies. System-level optical modelling of a
fiber-to-photonic-integrated-circuit coupling interface in Ansys Zemax
OpticStudio, with Lumerical FDTD used for the grating coupler.

**Status:** Stage 0 complete (literature selection and methodology). Zemax
implementation not yet started.

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

## Repository layout

```
literature/   paper screening and extracted parameters
zemax/        OpticStudio models (baseline, conventional, microlens, optimized)
python/       analysis, sweeps, tolerance maps, parameter studies, optimisation
data/         literature values, raw output, processed results
results/      coupling curves, tolerance maps, parameter studies, optimisation
docs/         methodology, validation, limitations
report/       final report
```

## Environment

Ansys Zemax OpticStudio 2026 R1.00 (Enterprise), driven through ZOS-API with
`pythonnet`; Ansys Lumerical FDTD via `lumapi`. Windows.

## Reproducing

Not yet reproducible — implementation begins at Stage 1. Each run will write a
`run_config.json` recording wavelength, software version, POP sampling and window
width, beam and receiver definitions, geometry, sweep ranges and raw POPD values.

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