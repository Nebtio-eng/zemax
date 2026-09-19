# LED Collimator & Uniform Illumination System — Methodology

**Status: FINALIZED — ARCHITECTURE-LIMITED, PROJECT STOPPED AFTER BUG-CORRECTED CAPTURE DIAGNOSTIC.**

## 1. Optical models

This repository contains the two Zemax OpticStudio NSC models developed for the project:

- `models/baseline_led_lens.zos` — baseline LED + single-lens + detector system.
- `models/two_lens_led_system.zos` — extended LED + two-lens + detector system.

The models are included specifically so the optical construction can be inspected directly in OpticStudio.

### Baseline single-lens system

- Object #1: Source Ellipse, Z = 0 mm
- Object #2: Detector Rectangle, Z = 50 mm
- Object #3: Standard Lens, Z = 20 mm, N-BK7
- Detector: 50 × 50 mm
- Detector sampling: 100 × 100
- Wavelength: 550 nm

### Two-lens system

- Object #1: Source Ellipse, Z = 0 mm
- Object #2: Detector Rectangle, Z = 60 mm
- Object #3: Standard Lens 1, N-BK7, Z = 15 mm
- Object #4: Standard Lens 2, N-BK7, Z = 30 mm
- Initial lens parameters were based on the baseline geometry.

## 2. Measurement and coverage definition

- Detector: **50 × 50 mm**
- Detector sampling: **100 × 100**
- Coverage threshold: **0.10 of peak**
- Coverage floor: **8.00%**

The corrected coverage metric uses the complete fixed detector grid. The active illuminated region is not allowed to shrink the evaluation region.

## 3. Wavelength audit

- Configured wavelength channels: **1**
- The ZRD capture denominator is the actual number of primary-ray records returned by the ZRD reader, not the requested `NumberOfAnalysisRays`.
- Captured-ray identity includes wavelength index plus ray number.

## 4. Corrected illumination metrics

Two CV definitions are retained for different purposes:

- **Fixed-region CV:** calculated over the complete 100 × 100 detector grid and used for feasibility/ranking.
- **Active-only CV:** calculated only over pixels above the coverage threshold and retained as a diagnostic.

The coverage peak is the **maximum of a 3 × 3 box-smoothed detector map**, rather than the raw maximum of one pixel.

Candidates below the **8.00% coverage floor** receive a fixed infeasibility penalty and are not called “best”. Feasible candidates are ranked by fixed-region CV, with centroid error as the secondary tiebreaker.

## 5. Ray-trace reproducibility

- Analysis rays: **50,000**
- 10× bare-source sanity rays: **500,000**
- Random seed: **12345**
- Wavelength channels: **1**
- Source seed property supported: **True**

The seed API is attempted where supported by the installed OpticStudio/ZOS-API version.

## 6. Geometry validation

The two-lens search includes explicit geometry checks:

- Lens 2 must be axially after Lens 1 plus Lens 1 thickness and a safety gap.
- The detector must be after Lens 2 plus Lens 2 thickness and a safety gap.
- `Edge1 >= abs(Clear1)` and `Edge2 >= abs(Clear2)`.
- Conic sag domain is checked.
- Front/rear axial separation is sampled over radius to reject self-intersecting lens geometries.
- Rejection reasons are logged.

## 7. Architecture capture diagnostic

Before running the expensive 21-variable optimization, the project imposed an architecture gate requiring **Lens-1 capture > 70.0%**.

Corrected results:

- Minimum Lens-1 capture: **8.9360%**
- Mean Lens-1 capture: **38.7463%**
- Capture-vs-coverage correlation: **0.5726335029337644**
- Required architecture gate: **> 70.0%**
- Conclusion: **capture/coverage-limited**

The corrected capture denominator uses the actual primary-ray population returned by the ZRD reader. The requested analysis-ray count is logged separately.

## 8. Bare-source sanity check

The normal bare-source sweep uses 50,000 rays. A separate 500,000-ray run was used only as a diagnostic.

An earlier implementation used the raw single-pixel maximum, which made the 10%-of-peak threshold overly sensitive to Monte Carlo shot noise. The corrected implementation uses the maximum of a 3 × 3 box-smoothed map.

## 9. Planned optimization that was gated

The planned search was:

- 150 LHS coarse proposals
- 150 multi-start fine evaluations
- fine refinement across the top 10 coarse candidates
- full 100 × 100 detector evaluation
- validation of the top 3 feasible candidates

This search was **not run** because the architecture failed the >70.0% Lens-1 capture gate.

## 10. Post-fix verification

### Bug 1 — capture-fraction denominator

**Symptom:** an intermediate implementation could produce physically impossible capture fractions above 100%.

**Cause:** the requested analysis-ray count was used as the denominator even when the ZRD reader returned a different number of primary-ray records.

**Fix:** normalize capture using the actual primary-ray records returned by the ZRD reader, while auditing wavelength channels and keeping requested ray count as a separate diagnostic.

**Verification:** corrected minimum and mean Lens-1 capture are **8.9360%** and **38.7463%**, respectively.

### Bug 2 — raw peak Monte Carlo shot noise

**Symptom:** several bare-source detector distances showed an identical raw maximum even while mean irradiance changed.

**Cause:** a single detector pixel controlled the 10%-of-peak threshold.

**Fix:** use the maximum of a **3 × 3 box-smoothed detector map**.

**Verification:** the 10× bare-source run remains diagnostic only, and the corrected peak definition is used for coverage decisions.

### Post-fix conclusion

After both fixes, the architectural conclusion did not change. The tested rotationally symmetric one-/two-lens standard/conic architecture is **capture/coverage-limited** at the tested working distances.

## 11. Interpretation

The project does not treat a low active-only CV over a small illuminated region as evidence of successful uniform illumination. The fixed-grid CV and coverage requirement prevent this shrink-to-uniform artifact.

The negative result is therefore an architectural finding rather than an optimizer failure.

## 12. Literature context

The result is consistent with the cited freeform LED illumination literature: freeform secondary optics and lens arrays provide substantially more control over LED emission redistribution than a small number of rotationally symmetric standard/conic surfaces.

See [`REFERENCES.md`](../REFERENCES.md).

## 13. Published result files

The repository includes:

- `results/bare_source_detector_map.png`
- `results/baseline_single_lens_map.png`
- `results/baseline_two_lens_map.png`
- `results/historical_8p58_single_lens_map.png`
- `results/old_flawed_4p02_map.png`
- `results/capture_fraction_vs_clear1.png`
- `results/bare_source_sweep.csv`
- `results/bare_source_sweep_sanity_10x.csv`
- `results/before_after_comparison.csv`
- `results/capture_diagnostic.csv`

The earlier 4.02% candidate is retained only as a documented flawed-objective case study and is not presented as a valid optimized design.
