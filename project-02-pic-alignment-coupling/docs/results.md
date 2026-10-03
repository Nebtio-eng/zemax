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

> 0.0733 lambda/n is not a value every coupling obeys. It is the **minimum** of
> the lateral x angular product, reached only when the arriving beam and the
> fiber mode are matched and flat. Any mismatch (size, or wavefront curvature)
> raises the product. It never lowers it.

A0's curvature mismatch is small (R = 149 um against a 51 um Rayleigh range at
20 um). It shortens the lateral tolerance by about 2% and lengthens the angular
one by about as much, so the product barely moves. A1's large mismatch pushes it
to 0.142, 48% above the minimum. The absolute-convention product for A0 (0.080)
is lower only because its 0.17 dB nominal loss uses up part of the 1 dB budget.
That is a property of the convention, not of the optics.

---

## 3. B against the paper

| Quantity | B (Zemax) | Paper | Difference |
|---|---|---|---|
| Lateral 1-dB | 8.140 um | +/-7 um | **+16.3%** |
| Angular 1-dB | 0.676 deg | +/-0.6 deg | **+12.6%** |
| Longitudinal 1-dB | 713 um (from 20 um), 702 um (from best gap) | 700 um | +1.9% / +0.3% |

The model reproduces its own analytic prediction to 0.2%, so these gaps sit
between the idealised model and the experiment, not inside the model. Candidate
explanations:

1. **Receiver MFD is ASSUMED** (34 um, matched to the beam). The paper does not
   report the TEC fiber's MFD numerically. A smaller real mode narrows the
   lateral tolerance (about 24 um MFD would give 7.0 um, at a 0.55 dB mismatch
   penalty) but *widens* the angular one. Untested.
2. **The real grating beam is elliptical**, not circular (`limitations.md`
   section 1). Measured tolerance is then axis-dependent; a circular model
   cannot reproduce that. Stage 8b tests it.
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
A smaller beam has a wider angular tolerance; that is the invariant again. No
choice of beam size can close both gaps at once: the paper's own pair,
7 um x 0.6 deg, gives a product of 0.0733 um rad, **24% below the theoretical
minimum** of 0.0960 for a peak-referenced (relative) 1-dB tolerance between
Gaussian modes at 1310 nm. In this model class that pair cannot be reached on
the relative convention. A product below the minimum *is* possible on the
absolute convention when there is nominal excess loss (A0's absolute product is
0.080), so one candidate is that the paper's tolerances are referenced to a
fixed loss level rather than to the peak. Others: the angular pivot,
non-Gaussian modes, or a different definition of +/-. All lie outside the
simulation. This is the strongest single argument
that the remaining discrepancy is about the experiment's definitions or real
hardware, not about the simulation. It is recorded as a limitation of the
source material in `limitations.md`.

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
> alone (31%, 1.6 dB) is ten times A0's entire nominal loss and is not counted by
> this model. The lens does not create alignment tolerance; it **converts
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
for both A0 and B, to four figures. This is the invariant from
`methodology.md` section 6 appearing in the architecture comparison. Gaussian
coupling conserves the product of positional and angular acceptance (a
phase-space, or etendue, argument), so beam expansion cannot enlarge it. It can
only change its shape.

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
  the real loss comparison favours **A0**. Whether the paper's lens was
  AR-coated has not been checked.
- **Grating directionality is excluded** for all three (`methodology.md`
  section 4). A backside configuration collects the *downward* grating
  emission, while A0 collects the upward one. How the two compare depends on
  the grating design; this is outside the Zemax boundary and is deferred to the
  FDTD stage.

---

## 6. Brief, section 20: the eleven questions

Questions as worded in `docs/brief-questions.md`. Questions 7 and 8 are answered
here from the Stage 9 parameter study (`docs/validation.md`, Stage 9;
data in `results/parameter_studies/`). The other nine are answered in a later
pass.

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
| 1 | **Lens radius** (and its twin, wafer thickness: 10 um of thickness acts like 6.6 um of radius) | +/-10% radius costs 0.19 / 0.10 dB at a 20 um gap. 1 dB at -21.7% / +40.2%. At a 700 um gap the absolute window is **+1.4%** (6.7 um) on the large side | `radius_alone/` |
| 2 | **Incidence angle** (grating emission angle in Si) | Optically benign once the fiber is re-pointed (0.033 dB at 6 deg), but the fiber must then be tilted 3.5x the incidence angle (21.5 deg at 6 deg). An un-re-pointed fiber passes 1 dB at 0.19 deg | `incidence/` |
| 3 | **Gap** | Negligible below 200 um (0.065 dB). The absolute lateral window falls from 8.1 to 2.6 um at 700 um | `gap/` |
| 4 | **Aperture** | Threshold: none at or above 80 um diameter; 0.1 dB at 50.7 um (D = 3.0 w); 1 dB at 35.7 um (D = 2.1 w) | `aperture/` |
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

The remaining questions (1-6, 9-11) will be answered against their wording in a
later pass.
