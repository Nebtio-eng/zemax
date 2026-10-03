# Limitations

What this model cannot claim, why the simplifications are nevertheless
defensible, and how each one is tested rather than asserted.

This document answers question 10 of the project brief: *what limitations
prevent us from claiming complete PIC-level electromagnetic accuracy?*

---

## 1. The central simplification

The grating coupler is not simulated. It is replaced by its equivalent
Gaussian output: a single beam waist radius of 4.6 um at 1310 nm, taken from
the primary paper.

Nothing in the Zemax model represents the grating period, etch depth, duty
cycle, fill factor, Bloch modes, or the directionality that determines how much
power leaves upward versus downward.

### Why this is legitimate

**Separation of scales is physical, not a convenience.** Once a field exists at
the grating plane, its propagation through the substrate is determined by the
wave equation acting on that field and the distance travelled. It does not
depend on the structure that produced it. The etch geometry determined *what
field emerged*; having emerged, it is no longer in the problem.

**The measured quantity is the one least sensitive to the omission.** Alignment
tolerance is a mode-overlap problem, set by beam width, wavefront curvature and
receiver mode size. The grating's internal physics sets the *absolute*
efficiency, not the *rate at which efficiency degrades with displacement*.
Absolute grating efficiency is explicitly excluded from this project's loss
budget (`methodology.md` section 4), so the omitted physics and the reported
quantity do not overlap.

**The primary authors used the same reduction.** Mangal et al. designed their
microlens by propagating an equivalent Gaussian in OpticStudio, and the
fabricated lens collimated the beam as predicted. The reduction is validated by
the hardware working.

**The closed-form version already reproduces three independent measurements.**
Gaussian overlap theory, containing no grating whatsoever, predicts 7.7 um
lateral against a measured +/-7 um, 0.72 deg angular against +/-0.6 deg, and
about 625 um longitudinal against 700 um. A simplification that still predicts
three independent axes to within roughly 20% is capturing the governing physics,
not evading it.

### Where it genuinely fails

A real grating beam is **not a circular Gaussian**. Light leaks out gradually
along the length of the grating, so the emitted beam is **elliptical** —
elongated along the in-plane propagation direction — with non-Gaussian tails and
a small wavefront curvature, emitted at an angle off normal rather than
straight down.

Consequence: **this model predicts equal tolerance in X and Y. The real device
is anisotropic.** Any X/Y asymmetry in the published or measured behaviour
cannot be reproduced by a circularly symmetric source, and must not be
explained away by tuning other parameters.

---

## 2. Accuracy, output by output

Not all results from this model deserve equal confidence. Stating which is
which is part of the result.

| Output | Confidence | Reason |
|---|---|---|
| Alignment tolerance (lateral, longitudinal, angular) | **Good** | Governed by beam geometry and mode overlap, which the model represents properly |
| Relative comparison between architectures A and B | **Good** | Both share the same approximations, so systematic error largely cancels |
| Absolute coupling loss | **Optimistic by construction** | Excludes grating directionality, grating scattering and waveguide/taper loss. Will be lower than the paper's measured insertion loss. This is expected and is not to be corrected by tuning. |
| X versus Y asymmetry | **Wrong** | Source is circularly symmetric; the real beam is elliptical. Addressed at Stage 8b. |
| Any claim about grating design | **Out of scope** | No grating physics in the model at all |

---

## 3. Other limitations

**Scalar diffraction.** POP propagates a scalar field and does not track
polarisation. Justified here: the 32 um expanded beam at 1310 nm diverges at
roughly 1.5 degrees, far inside the regime where scalar theory is accurate. This
would not hold for a tightly focused high-NA interface.

**Monochromatic.** Each POP run is at a single wavelength. Wavelength-dependent
grating directionality — a known weakness of grating couplers and the main
reason edge couplers are preferred for broadband WDM — is invisible to this
model. No bandwidth claim can be made from it.

**Direction and reciprocity.** The model runs PIC-to-fiber. The fiber-to-PIC
case is inferred by reciprocity, which holds for the tolerance curves but not
for absolute grating efficiency. Ansys documents that POP cannot cost the
fiber-to-grating direction directly, because the receiver mode would have to be
the grating's own field.

**No back-reflection or resonance.** POP as configured propagates forward only.
Fabry-Perot effects between chip facets, and reflection at the silicon-air
interface, are not modelled.

**Two assumed parameters.** The microlens clear aperture and the receiving
fiber's mode field diameter are not reported in the primary paper. Both are
recorded as `ASSUMED` in `literature/extracted_parameters.csv` with the
reasoning stated. Sensitivity to both must be checked before any conclusion
rests on them.

**Nominal, not manufactured, geometry.** Surface figure error, etch roughness,
lens sag deviation and adhesive effects are absent. The lens is a perfect
sphere. Real etched lenses are not.

---

## 4. How the central simplification is tested, not assumed

Two experiments convert this limitation from a caveat into a measurement.

**Test A — source sensitivity.** Vary the assumed source waist by +/-20% and
record the change in predicted 1-dB tolerance. If tolerance is insensitive, the
exact source shape does not matter and the reduction is demonstrably safe. If it
is sensitive, that sensitivity is itself a reportable result and bounds the
confidence in every tolerance figure.

**Test B — real source substitution (Stage 8b).** Compute the actual grating
field in Lumerical FDTD, export it, and substitute it for the idealised Gaussian
with every other element of the validated model unchanged. Any shift in the
tolerance curves is the error the simplification was concealing, measured in
microns and degrees rather than described in prose.

Expected outcome of Test B: tolerance becomes anisotropic, and the asymmetry
should be traceable to the grating length. Confirming that would validate both
the model boundary and the explanation for it.

---

## 5. What may and may not be said

Per section 19 of the project brief:

- **"Experimentally demonstrated"** — reserved for results reported in the cited
  literature.
- **"Numerically determined"** — this project's simulation output.
- Conclusions about *relative* benefit between architectures are better supported
  than conclusions about *absolute* loss.
- No statement in this project constitutes a validated claim about grating
  coupler design, PIC-level electromagnetic behaviour, polarisation, or optical
  bandwidth.
