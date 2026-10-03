# Project brief, section 20 — the eleven final questions

Reproduced verbatim from the original project specification, so they can be
answered against their exact wording. Answers belong in `docs/results.md`
section 6, each with a number and a pointer to the table row it came from.

The final report must answer:

1. What is the baseline fiber-to-PIC coupling loss?
2. How sensitive is the system to X/Y/Z misalignment?
3. What is the baseline 1-dB alignment tolerance?
4. Does the micro-lens improve alignment tolerance?
5. By how much?
6. What coupling-loss penalty, if any, is introduced?
7. Which lens parameters have the largest influence?
8. Is there a meaningful coupling-efficiency/alignment-tolerance trade-off?
9. Does the Zemax model reproduce the published behavior?
10. What limitations prevent us from claiming complete PIC-level
    electromagnetic accuracy?
11. What design would be worth investigating using Lumerical/FDTD next?

## Notes on scope

- Question 7 cannot be answered from Stage 8 data. It requires the Stage 9
  parameter study (radius of curvature, conic constant, aperture, substrate
  thickness, receiver MFD, gap, incidence angle), varying one at a time.
- Question 4 must be answered with the axis named. "Yes" is incomplete: lateral
  and area improve, angular degrades, longitudinal is the only net gain.
- Question 11 should name a specific design and the reason, not a general
  intention.
