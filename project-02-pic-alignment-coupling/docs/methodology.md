# Methodology — Zemax OpticStudio Model

Stage 0 deliverable. Defines what the simulation actually computes, where the
model boundary lies, and the gaps in the original project specification that had
to be resolved before any modelling can begin.

---

## 1. How coupling efficiency is computed

This project does **not** use a detector-power metric. It uses the Physical
Optics Propagation (POP) fiber-coupling integral, which is a mode-overlap
calculation.

OpticStudio reports:

```
eta_total  =  S  x  T
```

| Term | Name in OpticStudio | Meaning |
|---|---|---|
| `S` | System efficiency | Fraction of launched energy that survives the system (vignetting, aperture truncation, Fresnel losses if coatings are defined) |
| `T` | Receiver efficiency | Normalised overlap integral between the arriving complex field and the receiving mode |
| `eta_total` | Coupling efficiency | The product |

The Ansys application-gallery example shows this arithmetic explicitly:
`0.593864 x 0.66287 = 0.39365`, i.e. ~40%.

Coupling loss follows the project definition:

```
L = -10 * log10(eta)
```

`10 * log10` — not 20 — because `eta` is a **power** ratio, not a field ratio.

### POPD operands (for merit functions and Python automation)

| `Data` | Returns |
|---|---|
| 0 | Total coupling efficiency `eta_total` |
| 1 | System efficiency `S` |
| 2 | Receiver efficiency `T` |
| 10 | Beam waist size |
| 23 | Beam radius |
| 26 | M² |

**All three of 0, 1 and 2 must be logged on every run.** If tolerance degrades,
`S` vs `T` tells us immediately whether the cause is geometric (the beam is
walking off an aperture) or modal (the beam no longer matches the receiver). That
distinction is the difference between a result and an observation.

### Critical convention: the receiver mode is a radius, not a diameter

The POP **Fiber Data** tab takes `Waist X`/`Waist Y` as the **1/e² intensity
radius**. Fiber datasheets quote **mode field diameter**. Therefore:

```
Waist entry  =  MFD / 2
```

Getting this wrong by a factor of two is the single most common error in POP
fiber-coupling models and silently produces plausible-looking but wrong
efficiencies. It will be checked in Stage 3 against a trivial case: fiber to
identical fiber at zero distance must give `eta ≈ 1`.

---

## 2. Model boundary — what Zemax does and does not simulate

This follows §5 of the project brief and is confirmed by Ansys' own description
of the workflow: Lumerical handles "microscopic light interactions with the
grating coupler"; OpticStudio handles "macroscopic propagation and tolerancing".

### Zemax POP **does** model

- Propagation of a scalar complex field using **scalar diffraction theory**
- Gaussian mode launch and free-space / in-substrate diffraction
- Refraction at the micro-lens surface, including its aberrations
- Aperture truncation and vignetting
- Rigid-body misalignment of any element via Coordinate Breaks
- The overlap integral against a defined receiving mode

### Zemax POP **does not** model

- The grating itself: period, etch depth, duty cycle, Bloch modes, directionality
- Silicon waveguide modes or the taper/spot-size converter
- Vector/polarisation effects at sub-wavelength features; POP is **scalar**
- Back-reflection, resonance or wavelength-dependent grating directionality
- Anything where feature size approaches λ

### Consequence for the model, and how it is handled

The grating coupler is replaced by its **validated equivalent Gaussian source**:
a 4.6 µm waist radius at 1310 nm, taken directly from the primary paper. This is
exactly what the primary paper itself did to design the lens, so we are
reproducing their method, not improvising one.

**The scalar approximation is safe here, and we can say why.** The expanded beam
is 32 µm across at 1310 nm, giving a divergence half-angle of roughly
λ/(πw) ≈ 0.026 rad ≈ 1.5°. Scalar diffraction theory is accurate well beyond
this. The approximation would be questionable for a tightly focused high-NA
system; it is not questionable for this one. This is stated as a *justified*
approximation rather than an apologised-for one.

### Propagation direction: the model runs "OUT" (PIC → lens → fiber)

The Ansys documentation warns explicitly that the "IN" direction (fiber →
grating) cannot be accurately costed by POP alone, because the receiver mode
would have to be the grating coupler's own field profile, which POP does not
know. The vendor example resorts to handing the field back to FDTD for that
direction.

We therefore model the **OUT** direction and invoke **optical reciprocity** for
the fiber-to-PIC case. The alignment tolerance curve is reciprocal; the absolute
grating efficiency is not, and is excluded from our loss budget (see §4).

---

## 3. Gaps in the original specification, and their resolution

The project brief is strong but leaves eight things undefined. Each is resolved
here so that no parameter is silently invented later.

| # | Gap in the brief | Resolution |
|---|---|---|
| G1 | "Produce a measurable coupling quantity" — unspecified | POP fiber-coupling integral, `eta = S x T`, logged as POPD 0/1/2 |
| G2 | Does not choose edge coupler vs grating coupler | Grating / vertical backside architecture, to match the primary paper and the Ansys workflow |
| G3 | Does not say how the PIC interface is represented | Equivalent Gaussian mode, 4.6 µm waist, from the primary paper |
| G4 | Does not fix the propagation direction | OUT direction + reciprocity argument (§2) |
| G5 | "1-dB tolerance = L(0) + 1 dB" assumes the peak is at zero | After inserting the lens the peak moves. Nominal z is re-optimised first, *then* displacement is measured from the re-optimised peak. Both conventions will be reported in Stage 5 so the difference is visible, not hidden |
| G6 | POP sampling and window width not mentioned | Convergence test mandated: double the grid, require <0.5% change in `eta`. "Resample After Refraction" set on the micro-lens surface per Ansys guidance. Prop Report tab checked for warnings on every configuration |
| G7 | Assumes Python/ZOS-API is available | ZOS-API availability depends on licence tier. Fallback path defined: ZPL macro + Universal Plot CSV export. Note that the POP beam waist is *not* exposed to the Universal Plot, so any waist sweep needs ZPL or the API regardless |
| G8 | Angular tolerance treated as optional | It is mandatory. It is the axis that decides the research question (§5) |

---

## 4. Loss budget — what our number does and does not include

Our `eta` is the **packaging-interface coupling efficiency**: mode mismatch plus
propagation plus geometric misalignment. It excludes grating directionality,
grating scattering, and waveguide/taper loss.

This means **our absolute loss will be lower than the paper's measured insertion
loss**, and that is correct, not an error to be tuned away. Stage 4 will compare
the two with the grating contribution explicitly separated, and §19 of the brief
is respected: our figure is "numerically determined", the paper's is
"experimentally demonstrated".

---

## 5. Analytical ground truth (the project's third validation leg)

For two Gaussian modes of 1/e² intensity radius `w`, the normalised power overlap
under misalignment has closed-form solutions. Deriving these ourselves gives an
independent check that does not depend on Zemax being set up correctly.

**Lateral offset `d`:**

```
eta(d) = exp(-d^2 / w^2)
=>  d_1dB = 0.48 * w
```

(The coefficient 0.48 is confirmed independently by Gradkowski & O'Brien 2024.)

**Angular tilt `theta`:**

```
eta(theta) = exp(-(k * theta * w)^2 / 4),   k = 2*pi*n/lambda
=>  theta_1dB = 0.1528 * lambda / (n * w)
```

**Longitudinal offset `z`:**

```
eta(z) = 1 / (1 + (z / (2*z_R))^2),   z_R = pi * w^2 * n / lambda
=>  z_1dB = 1.018 * z_R
```

### Sanity check against the primary paper

| Quantity | Analytic prediction | Paper (measured) |
|---|---|---|
| Lateral 1-dB, 32 µm expanded beam (w = 16 µm) | 0.48 × 16 = **7.7 µm** | **±7 µm** |
| Angular 1-dB, same beam | 0.1528 × 1.31/16 = 0.0125 rad = **0.72°** | **±0.6°** |
| Longitudinal 1-dB, same beam (z_R = 614 µm) | **≈625 µm** | **700 µm** (JSTQE companion paper; the 2021 paper reports 0.2 dB over 300 µm with one fiber retracted, 1 dB at 400 µm with both) |

The analytic framework already reproduces all three published tolerances to
within ~20% before a single Zemax surface has been entered. This is the correct
place to start: it means that when the Zemax model disagrees, we will know the
disagreement is informative rather than a setup error.

---

## 6. The actual research question, sharpened

Combining the lateral and angular results above:

```
d_1dB * theta_1dB  =  (0.48 * w) * (0.1528 * lambda / (n * w))  =  0.0733 * lambda / n
```

**The beam radius `w` cancels.** The product of lateral and angular tolerance is
a constant set only by wavelength and index — it cannot be improved by beam
expansion. Expanding the beam buys lateral tolerance and spends angular tolerance
at exactly the same rate.

This is not new. The trade between offset tolerance, angular tolerance and
coupling efficiency for Gaussian beams is due to Joyce & DeLoach, "Alignment of
Gaussian beams", *Appl. Opt.* 23(23), 4187 (1984); the 1-dB tolerance
equations are in Gradkowski & O'Brien, *Appl. Opt.* 63(32), 8407 (2024). This
project **confirms it numerically** in POP. For unequal flat Gaussian modes the
product becomes `0.0733 lambda / (n sqrt(eta_0))`, with `eta_0` the
mode-mismatch efficiency: a larger product is paid for in efficiency.

So the honest answer to "does the micro-lens improve alignment tolerance?" is
**not** a yes. It is:

- **Lateral tolerance** improves **linearly** with the beam expansion factor
- **Longitudinal tolerance** improves **quadratically** (via `z_R ∝ w²`) — the
  largest and most under-reported gain
- **Angular tolerance** degrades **linearly** — the hidden cost
- **Nominal loss** gets slightly worse (lens aberration, truncation, an extra
  surface)

The project therefore has a real engineering answer: the micro-lens converts a
*positional* alignment problem into an *angular* one. Whether that is a win
depends entirely on which tolerance the packaging process can actually hold —
and mechanical fixtures hold angles far better than they hold microns, which is
why the approach works in industry despite the invariant above.

Stage 10's Pareto plot becomes meaningful rather than decorative: it is a genuine
multi-objective trade, not a search for one best number.

---

## 7. Reproducibility record (per §18 of the brief)

Every simulation run records: wavelength, OpticStudio version and build, POP grid
sampling in X and Y, analysis window width, beam definition type and waist,
receiver mode waist, "Resample After Refraction" state per surface, all surface
radii/thicknesses/materials, coordinate-break values, sweep ranges and step
sizes, optimisation operands and bounds, and the raw POPD 0/1/2 values. A
`run_config.json` is written alongside every results CSV.
