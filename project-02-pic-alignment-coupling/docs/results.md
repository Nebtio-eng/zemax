# Results — Stage 8 comparison and synthesis

Configurations A0 (conventional top-side), A1 (backside, flat Si exit) and B
(backside, monolithic Si micro-lens, R = -480 um), 1310 nm, 4.6 um source waist,
20 um nominal air gap. A0 uses a standard SMF (9.2 um MFD); A1 and B use a
34 um MFD TEC fiber. All values are **numerically determined** (Zemax POP)
unless marked otherwise. Paper values are **experimentally demonstrated**
(Mangal et al. 2021). How each number was obtained is in `docs/validation.md`;
the raw data is in `results/`.

---

## 1. Comparison table

"rel" = loss rises 1 dB above its value at the reference point (on-axis,
nominal gap). "abs" = total loss reaches 1.000 dB. Both are shown because the
nominal losses differ by 6 dB, and the relative number alone flatters whichever
configuration already has more loss.

| Quantity | A0 rel | A0 abs | A1 rel | A1 abs | B rel | B abs | Source |
|---|---|---|---|---|---|---|---|
| Nominal loss (dB) | 0.165 | | 6.255 | | 0.0004 | | M: Stage 3, Stage 7 |
| Lateral X 1-dB (um) | 2.250 | 2.055 | 6.188 | none | 8.140 | 8.138 | M: Stages 5, 7 |
| Lateral Y 1-dB (um) | 2.250 | 2.055 | 6.188 | none | 8.140 | 8.138 | M: Stages 6, 8 (radial cut at 90 deg) |
| Angular 1-dB (deg) | 2.446 | 2.234 | 1.318 | none | 0.676 | 0.676 | M: Stages 7, 8 |
| Longitudinal 1-dB from 20 um gap (um, gap increasing) | 36.3 | 31.6 | 240 | none | 713 | 713 | M: Stages 7, 8 |
| Longitudinal 1-dB from best gap (um) | 51.6 (best gap 0) | 51.6 | 249 (best gap 0) | none | 702 (best gap 30.9) | 702 | M: Stages 7, 8 |
| 2-D area A_1dB (um^2) | 15.90 | 13.27 | 120.3 | **0** | 208.2 | 208.1 | M: Stages 6, 8 (radial cuts) |
| Lateral x angular (um rad) | 0.0960 | 0.0801 | 0.1423 | - | 0.0960 | 0.0960 | D: from the measured rows |

M = measured in Zemax; D = derived from measured values. "none" = the
configuration never reaches 1 dB of loss or less, so no absolute window exists.
Every measured value agrees with the independent ABCD + overlap calculation to
within 0.2%.

Ratios, B against A0 (relative convention): lateral **3.62x**, angular
**1/3.62**, area **13.1x**, longitudinal **13.6x** from best gap (19.6x from the
20 um nominal).

---

## 2. Stage 8 predictions against measurements

| Prediction | Stated value | Measured | Difference |
|---|---|---|---|
| A0 angular 1-dB | 2.32 deg | **2.446 deg** (abs 2.234) | +5.4% |
| A0 longitudinal 1-dB | 52 um | **51.6 um** (from best gap) | -0.7% |
| B area A_1dB | ~208 um^2 | **208.2 um^2** | +0.07% |
| A1 area A_1dB | - | 120.3 um^2 rel, 0 abs | - |
| A0 invariant product | 0.0911 um rad (5% below ideal) | **0.0960 um rad** | +5.4%: the prediction is wrong |

**A0 angular (+5.4%).** The 2.32 deg prediction used the single-mode formula
with w = 4.94 um for both beam and fiber. The fiber mode is 4.6 um, and unequal
modes widen the angular window; the arriving wavefront curvature widens it a
little more. The exact overlap predicted 2.448 deg before measurement; Zemax
gives 2.446.

**The A0 invariant: A0 sits on the invariant, not 5% below it.** The ordering
asked for, A0 between B and A1, holds only trivially: A0 = 0.09602, B = 0.09601,
A1 = 0.1423. A0 and B are both on the ideal 0.0733 lambda = 0.0960 to four
figures. The premise that curvature pulls A0 *below* the ideal is wrong, for a
structural reason:

> 0.0733 lambda/n is the lateral x angular product of two **matched, flat**
> Gaussian modes (relative convention). The trade-off between offset
> tolerance, angular tolerance and coupling efficiency is due to Joyce &
> DeLoach (1984); the 1-dB form is in Gradkowski & O'Brien (2024). For
> unequal flat modes the product becomes 0.0733 lambda / sqrt(eta_0), with
> eta_0 the mode-mismatch efficiency: a mismatch buys a larger product only by
> losing efficiency (confirmed numerically to 0.012%, `validation.md`,
> "Reproduction of published results").

A0's curvature mismatch is small (R = 149 um against a 51 um Rayleigh range at
20 um). It shortens the lateral tolerance by about 2% and lengthens the angular
one by about as much, so the product barely moves. A1's large mismatch pushes it
to 0.142, 48% above the matched value. The absolute-convention product for A0 (0.080)
is lower only because its 0.17 dB nominal loss uses up part of the 1 dB budget.
That is a property of the convention, not of the optics.

---

## 3. B against the paper

| Quantity | B (Zemax) | Paper | Difference |
|---|---|---|---|
| Lateral 1-dB | 8.140 um | +/-7 um (X; +/-9 um along Y) | **+16.3%** (Y: -9.6%) |
| Angular 1-dB | 0.676 deg | +/-0.6 deg | **+12.6%** |
| Longitudinal 1-dB | 713 um (from 20 um), 702 um (from best gap) | 700 um (JSTQE 2020) | +1.9% / +0.3% |

Sources of the paper column: Mangal et al. 2021 (Opt. Express 29, 7601),
Fig. 10, for lateral (+/-7 um along the grating axis X, +/-9 um along Y) and
angular; the 1-dB longitudinal 700 um is from the companion JSTQE paper. The
2021 paper itself reports longitudinal differently (Fig. 11): **a 0.2 dB drop
over 300 um when either fibre alone is retracted**, and **1 dB at 400 um when
both fibres are retracted together** (two interfaces). Those two points are
compared in `validation.md`, "Reproduction of published results"
(single interface, B: 0.16 dB at 300 um).

The model reproduces its own analytic prediction to 0.2%, so these gaps sit
between the idealised model and the experiment, not inside the model. Candidate
explanations:

1. **Receiver MFD is ASSUMED** (34 um, matched to the beam). The paper does not
   report the TEC fiber's MFD numerically. A smaller real mode narrows the
   lateral tolerance (about 24 um MFD would give 7.0 um, at a 0.55 dB mismatch
   penalty) but *widens* the angular one. Untested.
2. **The real grating beam is elliptical**, not circular (`limitations.md`
   section 1). Measured tolerance is then axis-dependent; a circular model
   cannot reproduce that. **Tested in Stage 8b**: the FDTD beam's core is
   elliptical (3.54 x 4.23 um) and the tolerance is axis-dependent, but with X
   wider (8.91 um) than Y (8.54 um), the reverse of the paper's order, which
   follows from ASSUMED grating etch and width (`validation.md`, Stage 8b).
3. **Angular measurement pivot.** If the fiber was rotated about a point behind
   its facet, each tilt adds a lateral offset and the measured angular
   tolerance falls. The model tilts about the facet centre. Untested.
4. **Real lens imperfections** (sag error, roughness, decentre of the etched
   lens relative to the grating; the companion paper quotes a 2.5 um
   lens-to-grating 1-dB placement tolerance). Absent from the model.
5. **The two source papers are not self-consistent** (tested, below).

### Discrepancy 5: which literature pair was adopted

The 2021 paper reports a **32 um** expanded beam with the +/-7 um lateral
tolerance. The companion JSTQE paper reports the **480 um** lens radius with a
**630 um** substrate, which produces a **34 um** beam. The two cannot both hold:
a 32 um beam implies 592 um of silicon, which needs a 461 um lens to collimate.
B was built on the 630/480 pair (Stage 7). B was rebuilt on the 2021 pair,
**B592**: 592 um Si, R = -461 um, 32 um MFD TEC receiver
(`zemax/microlens/B592/`).

| | B (630 um, R 480, 34 um fiber) | B592 (592 um, R 461, 32 um fiber) | Paper |
|---|---|---|---|
| Beam radius at fiber | 16.92 um | 16.00 um | (32 um diameter, 2021) |
| Nominal loss | 0.0004 dB | 0.0010 dB | - |
| Lateral 1-dB | 8.140 um (+16.3%) | **7.678 um (+9.7%)** | +/-7 um |
| Angular 1-dB | 0.676 deg (+12.6%) | **0.717 deg (+19.4%)** | +/-0.6 deg |
| Lateral x angular | 0.0960 | 0.0960 | 0.0733 um rad (from +/-7 um x +/-0.6 deg) |

**Lateral lands at 7.68 um, as predicted (7.7 um).** So 6.6 of the 16.3
percentage points of lateral discrepancy come from which literature pair was
adopted, not from missing physics. The remaining 9.7 points are attributable to
explanations 1-4.

**This does not close the angular gap; it widens it, from +12.6% to +19.4%.**
A smaller beam has a wider angular tolerance; that is the lateral x angular
trade (Joyce & DeLoach 1984) again. The ideal matched Gaussian tolerances for
the paper's 32 um beam are **7.68 um and 0.72 deg**; the measured +/-7 um and
+/-0.6 deg are **91% and 84% of ideal**. Both measured values sit below the
ideal, so no choice of beam size reproduces both at once. Plausible reasons,
all outside the simulation: the angular pivot (candidate 3), real lens and
placement imperfections (candidate 4), the non-Gaussian grating beam
(candidate 2), or tolerances referenced to a fixed loss level rather than to
the peak (on the absolute convention a nominal excess loss shrinks both
windows; A0's absolute product is 0.080 against 0.096 relative). Recorded in
`limitations.md`.

### Status: closed as far as this model can close it (after Stage 9)

Two independent attempts to close the lateral gap both reach the paper's lateral
value, and both fail the same way:

| Attempt | Lateral 1-dB | vs paper 7 um | Angular 1-dB | vs paper 0.6 deg | Source |
|---|---|---|---|---|---|
| B as built (630 um, 34 um fibre) | 8.140 um | +16.3% | 0.676 deg | +12.6% | Stage 7 |
| B592: smaller beam (2021-paper pair) | 7.678 um | +9.7% | 0.717 deg | **+19.4%** | Stage 8 |
| B with a 24 um receiver (candidate 1) | 7.040 um | **+0.6%** | 0.828 deg | **+38.0%** | Stage 9, `receiver_mfd/` |

Lateral improves and angular worsens by the same factor each time: that is the
lateral x angular trade. With the paper's measured pair at 91% (lateral) and
84% (angular) of the ideal for its own 32 um beam, the remaining discrepancy is
attributed to the experiment's definitions, measurement geometry or hardware,
not to a parameter this Gaussian model omits. No further parameter hunting is
warranted.

**Stage 8b (the real grating beam, candidate 2) does not close it either.**
Against the paper's per-axis values the FDTD beam moves Y to within 5%
(8.54 vs 9 um) and angular to +4.5% / +7.3% (from +12.6%), but moves X further
away (8.91 vs 7 um, +27%), because the axis order comes out reversed
(`validation.md`, Stage 8b).

---

## 4. Headline claim

> For a 1310 nm grating-coupler output (4.6 um waist), a monolithic backside
> silicon micro-lens (B) with a matched 34 um MFD fiber widens the 1-dB lateral
> tolerance from **+/-2.25 um to +/-8.14 um (3.6x)**, the 1-dB X-Y area from
> **15.9 to 208 um^2 (13x)**, and the 1-dB working distance from **52 um to
> about 700 um (14x)**. **The cost:** angular tolerance falls from
> **+/-2.45 deg to +/-0.68 deg (3.6x)**, a large-mode TEC receiver becomes
> mandatory, the chip needs a backside lens process, and the silicon-air exit
> face must be anti-reflection coated. Without that coating, Fresnel reflection
> alone (31%, 1.6 dB at normal incidence; 1.85 dB expected and 2 dB recovered
> with a 170 nm SiN coating in Mangal et al. 2021) is ten times A0's entire
> nominal loss and is not counted by this model. The lens does not create alignment tolerance; it **converts
> positional tolerance into angular tolerance at a fixed exchange rate**. The
> only net gain is longitudinal. The substrate alone (A1) is not a working
> interface: 6.26 dB of loss and no 1-dB window at all. **The lens's
> contribution is phase, not size.**

No winner is declared on one metric. B is better if the assembly process holds
angle better than position (the usual case: flat-polished facets and planar
bonding control tilt mechanically, while lateral placement needs active
alignment). A0 is better if angle is the weaker axis, or if the TEC fiber and
backside processing are not available.

---

## 5. Two points an examiner will attack

### 5.1 "B's 3.6x lateral gain is cancelled by a 3.6x angular loss"

**Correct, and it is the central result, not a weakness to hide.** Measured:
lateral x3.618, angular /3.619. The lateral x angular product is 0.0960 um rad
for both A0 and B, to four figures. This is the invariant of
`methodology.md` section 6 appearing in the architecture comparison. It is not
a finding of this project: the trade between offset tolerance, angular
tolerance and coupling efficiency for Gaussian beams is due to Joyce & DeLoach
(Appl. Opt. 23, 4187, 1984), and its 1-dB form is given by Gradkowski &
O'Brien (Appl. Opt. 63, 8407, 2024). This project **confirms it numerically**
in POP. Gaussian coupling conserves the product of positional and angular
acceptance (a phase-space, or etendue, argument), so beam expansion cannot
enlarge it. It can only change its shape.

What survives the objection:

- **The longitudinal axis is not covered by the invariant.** It scales with the
  Rayleigh range, which goes as w^2, so it grows 13.6x while lateral grows only
  3.6x. This is a real, net gain.
- **Which axis matters is a property of the assembly process, not the optics.**
  A tolerance moved from a hard axis (sub-micron lateral placement) to an easier
  one (tilt held by flat reference surfaces) is an engineering gain even though
  the optical product is unchanged. The claim must be stated this way: "the
  lens moves the tolerance to the axis packaging controls best", not "the lens
  improves alignment tolerance".
- **The 2-D area (13x) does not contradict this.** A_1dB counts only the
  positional plane, so it grows as w^2 while angular acceptance shrinks.

### 5.2 "B's loss advantage over A0 is an artefact of the 20 um gap"

**Largely correct, and the comparison should say so.** A0's 0.165 dB at 20 um is
entirely beam spreading across the gap. At zero gap A0 is 0.000 dB (Stage 5),
the same as B. On nominal loss **there is no meaningful difference between A0
and B** in this model. What the 20 um gap actually exposes is the
**longitudinal** difference: A0 loses 1 dB by a 52 um gap, B not until about
700 um. So the defensible statement is "B tolerates the gap; A0 needs near
contact", not "B has lower loss".

Two further corrections both go against B, and must be stated alongside it:

- **Fresnel reflection is not modelled** (POP system efficiency S = 1.000
  everywhere). B's uncoated Si-air exit face reflects 30.9% (1.61 dB) at normal
  incidence. A0's oxide cladding exit reflects about 3.3% (0.15 dB). Uncoated,
  the real loss comparison favours **A0**. The paper confirms the size of the
  effect experimentally: of the extra loss of its backside interface, **1.85 dB
  was expected from Fresnel reflection at the silicon-air surface**, and a
  **170 nm SiN anti-reflection coating improved the measured coupling by
  2 dB** (Mangal et al. 2021, Section 5). Our 1.61 dB is for normal incidence
  on one face; the paper attributes its slightly larger value to off-normal
  incidence on the lens.
- **Grating directionality is excluded** for all three (`methodology.md`
  section 4). A backside configuration collects the *downward* grating
  emission, while A0 collects the upward one. How the two compare depends on
  the grating design; this is outside the Zemax boundary and is deferred to the
  FDTD stage.

---

## 6. Brief, section 20: the eleven questions

Questions as worded in `docs/brief-questions.md`. Questions 7 and 8 are answered
here from the Stage 9 parameter study (`docs/validation.md`, Stage 9;
data in `results/parameter_studies/`); question 11 from Stage 8b. The other
eight are answered in a later pass.

### Q7. Which lens parameters have the largest influence?

"Influence" splits into two questions with different answers: which parameters
you **choose** (design levers), and which you must **hold** (manufacturing
sensitivities). Ranking by one yardstick alone gives a misleading order.

**Design levers: what sets the tolerances you get.**

| Rank | Parameter | Effect over a plausible range | Source |
|---|---|---|---|
| 1 | **Substrate thickness, with the lens radius re-matched** (and the receiver matched) | Lateral 5.44 -> 12.61 um and angular 1.01 -> 0.44 deg for 400 -> 1000 um of Si, loss < 0.005 dB, product fixed at 0.0960. Sets *where* the tolerance sits, never how much there is in total | `thickness_rematched_receiver_matched/` |
| 2 | **Receiver MFD** | Must match the beam: +/-10% MFD costs 0.044 dB; 24 or 48 um costs 0.5 dB. A mismatched receiver widens the relative lateral window (6.67 -> 9.97 um for 20 -> 48 um) but shrinks the absolute one | `receiver_mfd/` |

**Manufacturing sensitivities: what must be controlled.**

| Rank | Parameter | Sensitivity | Source |
|---|---|---|---|
| 1 | **Lens radius** (and its twin, wafer thickness: 10 um of thickness acts like 6.6 um of radius) | +/-10% radius costs 0.19 / 0.10 dB at a 20 um gap. 1 dB at -21.7% / +40.2% (nearly symmetric in surface power, +28% / -29%; section 7.2; independently reproduces Mangal et al. 2021 Fig. 4, -22.7% / +40.9%). At a 700 um gap the absolute window is **+1.4%** (6.7 um) on the large side | `radius_alone/` |
| 2 | **Incidence angle** (grating emission angle in Si) | Optically benign once the fiber is re-pointed (0.033 dB at 6 deg; 4th by loss alone, section 7.1), but ranked here by its system consequence: the fiber must be tilted 3.5x the incidence angle (21.5 deg at 6 deg), and an un-re-pointed fiber passes 1 dB at 0.19 deg | `incidence/` |
| 3 | **Gap** | Negligible below 200 um (0.065 dB). The absolute lateral window falls from 8.1 to 2.6 um at 700 um | `gap/` |
| 4 | **Aperture** | Threshold: none at or above 80 um diameter (radius 2.4 w, not the 3 w rule of thumb); 0.1 dB at 50.7 um; 1 dB at 35.7 um (radius 1.05 w) | `aperture/` |
| 5 | **Conic constant** | Negligible: |k| <= 10 changes loss by <= 0.0015 dB and lateral by <= 0.005 um. **A spherical lens is sufficient; no aspheric correction is needed** | `conic/` |
| - | **Lens material** | Not a free parameter (monolithic silicon). A hypothetical +/-1% index moves lateral by +/-0.04 um | analytic, `validation.md` Stage 9 |

**Corrections to the predicted ranking** (thickness, MFD, radius, gap, aperture,
conic, material, incidence):

- Thickness with the radius re-matched is the top *design* lever but not a
  sensitivity: re-matching removes the loss penalty by construction. As a
  *tolerance*, thickness acts through the radius-matching condition and ranks
  with the radius.
- The lens radius is the top *manufacturing* sensitivity, ahead of the receiver
  MFD: +/-10% costs about 0.1-0.2 dB, against 0.04 dB for the MFD. Its
  criticality grows sharply with working distance.
- Incidence angle ranks higher than predicted, not for its optical penalty but
  because it forces a fiber tilt of 3.5x the incidence angle.
- The aperture threshold sits at about 36-51 um diameter, not 102 um.
- Conic is confirmed negligible.

### Q8. Is there a meaningful coupling-efficiency / alignment-tolerance trade-off?

**Only if the receiver is mismatched, and then it is an unfavourable trade.
The fundamental trade-off is between lateral and angular tolerance, not
between efficiency and tolerance.**

- **With the receiver matched, there is no efficiency cost to tolerance.**
  Scaling the design (thickness, radius and receiver together) moves lateral
  tolerance from 5.44 to 12.61 um at < 0.005 dB loss throughout
  (`thickness_rematched_receiver_matched/`). What it costs is angular
  tolerance: the lateral x angular product stays at 0.0960 um rad at every
  point.
- **Spending loss to buy tolerance works only on the relative convention.**
  Oversizing the receiver from 34 to 48 um MFD raises the relative lateral
  tolerance from 8.14 to 9.97 um (+22%) at a cost of 0.52 dB. On the absolute
  convention the same change *reduces* lateral tolerance from 8.14 to 6.91 um
  and angular from 0.676 to 0.406 deg (`receiver_mfd/`). The extra 1-dB window
  is spent on the loss that bought it.
- **Every deliberate mismatch tested raises the product above 0.0960 and adds
  loss**: receiver (0.110 at 20 um MFD, 1.15 dB), radius error (0.117 at 300
  um, 3.42 dB), aperture (0.132 at 20 um, 6.0 dB). There is no parameter
  setting that buys total tolerance with efficiency.
- **Across the architectures** (section 1): B's gain over A0 is also a
  redistribution (lateral x3.6, angular /3.6), not an efficiency trade. The
  only configuration with a large efficiency penalty, A1 (6.26 dB), is *less*
  tolerant on every absolute measure.

So the meaningful trade-off in this system is **positional against angular
tolerance**, set by the beam size, at a fixed product. Efficiency is a separate
axis that any mismatch degrades and nothing in the design space converts into
tolerance.

### Q11. What design would be worth investigating using Lumerical/FDTD next?

**Answered by doing it (Stage 8b, `docs/validation.md`), and the result names
the next design: an apodized grating whose emission is Gaussian, with a round
core.**

What was done: a 3-D FDTD of a uniform 20-period grating (all grating
parameters ASSUMED) supplied the real emitted field as the POP source of B,
with every Zemax number predicted first. Numerically determined:

- **Loss at the best fibre position: 1.20 dB** (B with a Gaussian source:
  0.0004 dB). Of this, **1.04 dB is the non-Gaussian shape** of the emitted
  field (overlap with its own best-fit elliptical Gaussian 0.787) and
  **0.16 dB is the ellipticity** of that Gaussian core (3.54 x 4.23 um).
  Walk-off and tilt cost 0.0004 dB; the lens position (on axis or centred on
  the beam) changes the loss by at most 0.026 dB.
- **Tolerances barely change**: lateral x 8.91, y 8.54 um; angular 0.63 /
  0.64 deg; longitudinal 663 um (relative to the 1.20 dB peak; no absolute
  1-dB window exists, because the peak is already above 1 dB). The 1-dB area
  is 14.8% larger than B's.

So the grating, not the micro-lens, now limits the interface, and it limits
**efficiency, not tolerance**. The design worth investigating next is
therefore:

1. **An apodized grating** (period or fill factor varying along x) that emits a
   Gaussian profile, to recover most of the 1.04 dB. Mangal et al. (2021)
   attribute their remaining ~1 dB to the same cause and project < 2 dB per
   interface with apodized gratings; this work supports that attribution
   numerically.
2. **Grating width and etch depth chosen for a round core**, to recover the
   0.16 dB and to set which axis gets the wider tolerance (here X; in the
   paper Y).
3. **A bottom reflector or directionality optimisation**, which this project
   does not model at all (POP sees only the downward beam; the paper reports
   -7.5 dB backside grating efficiency without a metal reflector against
   -2.3 dB with one).

The pipeline built for Stage 8b (`fdtd_source.py` -> .zbf -> B_fdtd ->
predictor) evaluates any such grating without changes to the Zemax model.

The remaining questions (1-6, 9-10) will be answered against their wording in a
later pass. The Stage 9 evidence behind Q7 and Q8 is set out in section 7.

---

## 7. Stage 9 — what each lens parameter does (configuration B)

One parameter changed at a time from B. Data: `results/parameter_studies/`;
method, convergence and full tables: `docs/validation.md`, Stage 9. All values
numerically determined.

### 7.1 Sensitivity ranking

Ranked by how far each parameter moves loss or tolerance over its plausible
range. "Design lever" = something you choose; "sensitivity" = something you must
hold in manufacture.

| Rank | Parameter | Kind | Measured effect | Prediction |
|---|---|---|---|---|
| 1 | Substrate thickness, ROC re-matched, receiver matched | design lever | 400 -> 1000 um: lateral 5.44 -> 12.61 um, angular 1.01 -> 0.44 deg, working distance 294 -> 1667 um, loss < 0.005 dB, product 0.0960 throughout | right |
| 2 | Receiver MFD | design lever | 20 um: 1.15 dB / 6.67 um; 34 um: 0.0004 dB / 8.14 um; 48 um: 0.52 dB / 9.97 um. +/-10% costs 0.044 dB | right, to 2 dp |
| 3 | Lens ROC alone (630 um fixed) | sensitivity | 1 dB at **-21.7% / +40.2%** at a 20 um gap; +/-10% costs 0.19 / 0.10 dB. Absolute window **+1.4%** at a 700 um gap | magnitude right, **symmetry wrong** |
| 4 | Incidence angle (in Si) | sensitivity | Residual 0.033 dB at 6 deg once the fibre is re-pointed; the fibre must tilt by the Snell angle (3.5x) | right on loss |
| 5 | Gap | sensitivity | 0.011 dB at 100 um, 0.065 dB at 200 um; absolute lateral window 8.1 -> 2.6 um at 700 um | right |
| 6 | Aperture diameter | sensitivity | Flat at >= 80 um (radius **2.4 w**); 0.11 dB at 50 um; 1 dB at 35.7 um | threshold right, **value too conservative** |
| 7 | Conic constant | sensitivity | -10 to +10: loss 0.0002 -> 0.0019 dB, lateral 8.142 -> 8.135 um | right |
| 8 | Lens material | not free | Monolithic Si; a hypothetical +/-1% index moves lateral by +/-0.04 um | - |

Incidence angle is optically 4th, but it carries a system cost the ranking by
loss hides: the fibre must be tilted by 3.5x the incidence angle, and an
un-re-pointed fibre passes 1 dB at only 0.19 deg (Q7 answer below).

### 7.2 Why each trend looks the way it does

- **Thickness (ROC re-matched).** A thicker wafer gives the beam more distance to
  spread before the lens collimates it, so the beam leaving the chip is wider.
  A wider matched beam forgives more sideways offset (lateral proportional to
  w), forgives less tilt (angular proportional to 1/w), and stays collimated
  over a longer distance (longitudinal proportional to w^2). Lateral x angular
  cannot change: 0.0960 at every thickness. With a *fixed* 34 um receiver the
  same sweep also adds size-mismatch loss away from 630 um (0.70 dB at 400 um,
  0.80 dB at 1000 um) and lifts the product to 0.104-0.105.
- **Receiver MFD.** The fibre mode is the target the beam must land on. A bigger
  target forgives more offset; a smaller one forgives more tilt. Any size
  mismatch costs loss as 4 w1^2 w2^2 / (w1^2 + w2^2)^2, which is flat at the
  match and rises on both sides (hence 0.044 dB for +/-10%).
- **ROC alone: why asymmetric.** With the substrate fixed, the wrong ROC leaves
  the beam converging (ROC too small: lens too strong) or diverging (ROC too
  large). The two errors are not mirror images. Surface power goes as 1/R, so
  -21.7% of radius is +28% of power, while +40.2% of radius is only -29% of
  power: in *power* the 1-dB window is nearly symmetric (+28% / -29%), and the
  asymmetry is an artefact of measuring a 1/R quantity in R. A too-strong lens
  also pulls the waist out to several hundred microns, so it *lengthens* the
  working distance (966 um at 420 um ROC) while a too-weak lens shortens it
  (293 um at 700 um). This is the origin of Stage 10's one-sided trade.
  **Both results are published.** Mangal et al. (2021) report the ROC window
  from their own OpticStudio model, -100 / +180 um about a 440 um nominal
  (-22.7% / +40.9%, fibre 100 um from the lens; Fig. 4), and state the
  one-sided rule in section 2(v): a larger ROC gives a divergent beam with no
  benefit, a smaller one a convergent beam usable for a longer working
  distance. This project's numbers are an **independent reproduction**: at the
  paper's conditions our model gives -94 / +186 um (`validation.md`,
  "Reproduction of published results").
- **Gap.** B's beam is collimated with a ~686 um Rayleigh range, so nothing
  changes until the gap becomes a sizeable fraction of that. At 700 um the beam
  has grown and curved enough that B starts at 0.92 dB, and the absolute window
  collapses.
- **Aperture: a threshold, not a slope.** A Gaussian carries almost no power
  beyond about 2w. An aperture is invisible until it cuts into the beam, then
  removes power (S falls) *and* distorts the mode (T falls) together:
  0.002 dB at 70 um, 0.11 dB at 50 um, 1.1 dB at 35 um, 6 dB at 20 um. The
  truncated-Gaussian overlap formula predicted every point within 0.015 dB.
- **Conic.** The beam (radius 17 um) uses only a small cap of a 480 um-radius
  lens. Over that cap a conic changes the surface height by about k x 0.0015 um
  at twice the beam radius, about lambda/346 of path error per unit k. Negative
  k very slightly helps (0.0002 dB at k = -10) by trimming the sphere's residual
  aberration.
- **Incidence.** Off-axis use of a sphere adds coma and astigmatism, visible as a
  slight left-right asymmetry in the lateral scan, but at 0.033 dB at 6 deg it
  is negligible. The large effect is geometric: Snell refraction from silicon
  (n = 3.5) multiplies the angle by about 3.5 on exit.

### 7.3 Predictions: where they were right and wrong

Right: the thickness and receiver trends (to 2 dp), the conic null, the gap
behaviour, the incidence penalty, and the existence of an aperture threshold.
Every Gaussian-physics point agreed with its ABCD prediction to within
0.025 dB.

Three corrections:

1. **The lens ROC tolerance is asymmetric: 1 dB at -21.7% / +40.2%**, not a
   symmetric ~30%. The asymmetry comes from measuring a 1/R quantity in R; in
   surface power the window is nearly symmetric (+28% / -29%). The asymmetry
   was a wrong prediction here but is not new: Mangal et al. (2021, Fig. 4)
   report -22.7% / +40.9%.
2. **The aperture threshold sits at 2.4x the beam radius** (flat above 80 um
   diameter, w = 16.9 um), **not the 3x rule of thumb** (102 um) quoted
   beforehand. 3x is safe but over-conservative. The 1-dB point is at 1.05x the
   beam radius (35.7 um diameter).
3. **The conic constant does nothing**, and that is a manufacturing
   recommendation, not a dull result: **a spherical lens is optically
   sufficient; no aspheric correction is required.** Lens-process effort should
   go into ROC control (rank 3), not surface shape.

How to read the product. For flat Gaussian modes the relative-convention
product is 0.0960 um rad / sqrt(eta_0), where eta_0 is the size-mismatch
efficiency (Joyce & DeLoach 1984; verified to 0.012% across the receiver
sweep, `validation.md`, "Reproduction of published results"). A larger
product is therefore paid for in efficiency. The relation is exact only for
Gaussian modes: in the aperture sweep a slightly truncated (non-Gaussian) beam
gave 0.0958 at 50-60 um diameter, 0.2% below the matched value.

---

## 8. Stage 10 — Pareto front and the passive-alignment design

Script: `python/optimization.py`. Outputs: `results/optimization/`
(`design_grid.csv`, `pareto_front.csv`, `run_config.json`, three plots). All
values numerically determined (POP); every point was predicted first (ABCD +
overlap; worst disagreement 0.010 dB in loss and 0.11 um in lateral tolerance).

### 8.1 Why two dimensions

Stage 9 showed conic, aperture (above its 80 um threshold), gap (below
~200 um), incidence angle (once the fibre is re-pointed) and material to be
flat. Searching them would add dimensions with no effect. The optimisation is
over the only two real freedoms:

- **Beam size**: substrate thickness t, with the receiver MFD matched (2 w(t))
  and the ROC scaled with it. This chooses *where on the lateral/angular
  invariant* the design sits.
- **rho = ROC / collimating ROC**: below 1 the lens is slightly focusing, which
  buys working distance at a loss (the convergent-beam option of Mangal et al.
  2021, section 2(v), quantified here).

Grid: t = 400, 500, 630, 700, 786, 800, 900, 1000, 1200 um x rho = 0.80-1.05
(54 designs) plus the two Stage 9 reference points. Objectives (absolute
convention): minimum loss; maximum lateral, angular and longitudinal 1-dB and
area A_1dB. The front uses epsilon-dominance (1% on tolerances, 0.01 dB on
loss) so that differences below measurement significance do not keep a design
alive. Robustness: worst loss under a +/-5% ROC error (5% is ASSUMED; no
manufacturing tolerance is reported in the sources). Convergence at the three
grid corners: worst change 0.0074%.

### 8.2 Validation gate

| Check | Result |
|---|---|
| Stage 9 point (630 um, ROC 420 um, 34 um fibre) reproduced | **0.2975 dB, 964.8 um** (Stage 9: 0.2975 dB, 964.8 um) |
| Below-collimation branch on the front | **yes: 18 of 25 front designs have rho < 1** |
| Above-collimation designs on the front | **none**: rho > 1 is dominated, as Stage 9 implied |

**The specific Stage 9 point is itself dominated** once thickness is free, by
t = 700 um, rho = 0.95: longer working distance (981 vs 864 um, absolute), at
one eighth of the loss (0.037 vs 0.298 dB), with more lateral (8.75 vs
6.57 um) and more angular tolerance (0.606 vs 0.59 deg). The branch is real;
it is cheaper to use it on a thicker wafer than on the 630 um one.

Working distance does not grow without limit as rho falls. It saturates or
peaks at rho = 0.80-0.90 (630 um: 976 um, plateau over 0.80-0.85; 800 um:
1557 um at 0.85; 1000 um: 2388 um at 0.90) and then falls, because the lens
starts focusing the beam too tightly. Past the peak a design pays more loss for
the same or *less* distance.

### 8.3 The invariant across the front

Every matched, collimated design sits on the matched value: lateral x angular
= 0.09600-0.09601 um rad against the exact (ln10/10) lambda/pi =
0.096015 um rad, i.e. -0.011% to -0.002%. That is within the 1e-4 precision of
POP's tan(theta) tilt convention (`docs/zosapi_gotchas.md` #30). The focusing
(rho < 1) front designs sit slightly above it: +0.003% to +0.22%, rising as rho
falls. Their curvature mismatch is real but small. The largest excess on the
whole grid is +10.3%, at a dominated rho = 0.80 design. Note:
`design_grid.csv` column `product_excess_pct` was computed against the rounded
0.0733 lambda (0.096023); the figures here use the exact value.

### 8.4 Engineering criterion: passive alignment

Passively aligned fibre interfaces (registration features, alignment detents)
are reported at about +/-10 um in all three axes, a few microns in the
best-constrained direction (prompt; not yet in `literature/papers.md`). The
criterion is therefore an **absolute** 1-dB lateral tolerance of at least
10 um.

| Design | t (um) | ROC (um) | Fibre MFD (um) | Loss (dB) | Lateral abs (um) | Angular abs (deg) | Working distance (um) | Worst loss, ROC +/-5% |
|---|---|---|---|---|---|---|---|---|
| A0 (top-side) | - | - | 9.2 | 0.165 | 2.06 | 2.23 | 32 (from 20 um) | - |
| B (as built) | 630 | 480 | 34 | 0.0004 | 8.14 | 0.676 | 713 | 0.046 dB |
| **Threshold design** | **786** | **590** | **41.7** | **0.0004** | **10.01** | **0.550** | **1042** | **0.057 dB** |
| Predicted design | 800 | 600 | 42.4 | 0.0004 | 10.17 | 0.541 | 1078 | 0.058 dB |
| Thick design | 1000 | 737 | 52.6 | 0.0002 | 12.61 | 0.436 | 1666 | 0.087 dB |
| Thickest tested | 1200 | 876 | 62.8 | 0.0001 | 15.06 | 0.365 | 2386 | 0.129 dB |

**Predictions tested.** 800 um, ~42 um MFD, ~600 um ROC: predicted ~10.2 um
lateral, ~0 dB, ~1.1 mm working distance, ~0.54 deg. Measured 10.17 um,
0.0004 dB, 1.078 mm, 0.541 deg: **confirmed**. 1000 um: predicted 12.6 um and
0.44 deg; measured 12.61 um and 0.436 deg: **confirmed**.

**A0 cannot be passively aligned** (2.06 um absolute, five times short).
**B as built is close** (8.14 um). The threshold is reached at **786 um of
silicon**, and every thicker design exceeds it.

**What it costs**, in order of size:

1. **Angular tolerance**: 0.676 -> 0.550 deg at the threshold (0.436 at
   1000 um). The passive-alignment features must therefore also hold fibre
   tilt to about +/-0.55 deg (+/-9.6 mrad). Whether a given scheme does is a
   requirement to check, not something this model knows.
2. **A larger-mode receiver**: 41.7 um MFD instead of 34 um, which must exist
   as a thermally expanded core fibre (not verified here).
3. **Tighter ROC control**: the worst loss for a +/-5% ROC error rises from
   0.046 dB (B) to 0.057 dB (786 um) and 0.129 dB (1200 um). A longer working
   distance amplifies de-collimation.
4. **A thicker chip** (786 vs 630 um), and the same anti-reflection coating
   requirement as B (Fresnel is not counted by POP).

**Not a cost: loss.** Every collimated, matched design on the thickness axis
is at <= 0.005 dB, so the lateral gain is not bought with efficiency.

**Headline.** Passive alignment at +/-10 um is reachable by using a thicker
wafer (>= 786 um, with ROC and fibre scaled to it). It is paid for entirely in
angular tolerance (about -19% against B) and in receiver mode size, not in
coupling loss. The invariant guarantees there is no design that improves
lateral and angular tolerance together.
