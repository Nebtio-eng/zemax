# Zemax Optical Design & Simulation

A portfolio of optical design and simulation projects developed with Ansys Zemax OpticStudio, NSC/sequential modeling, ZOS-API, and Python automation.

## Research areas

### Integrated Photonics

Projects involving optical component design, beam shaping, collimation, optical interfaces, photonic systems, and computational optical design.

### Fusion Optics

Projects involving laser beam delivery, diagnostic optics, imaging, spectroscopy, and optical simulation relevant to fusion research.

## Projects

### Integrated Photonics

| Project | Description |
| --- | --- |
| [Project 01 — LED collimator](project-01-integrated-photonics/project-01-led-collimator/) | Zemax NSC model of a single-lens and two-lens LED collimator, with ZOS-API automation and detector analysis. |
| [Project 02 — PIC alignment coupling](project-02-pic-alignment-coupling/) | Alignment-tolerant fiber-to-PIC coupling through a backside silicon micro-lens, modelled in Zemax OpticStudio POP with a 3-D Lumerical FDTD grating beam as the source. Quantifies how beam expansion trades lateral for angular tolerance (lateral 2.25 → 8.14 um, angular 2.45 → 0.68 deg) and attributes 1.04 dB of the real-grating loss to its non-Gaussian profile; published results reproduced, with differences reported. Status: Stages 1-11 and 8b complete; final report pending. |

### Fusion Optics

Future projects will be added here as they are developed.

## Tools

- Ansys Zemax OpticStudio
- Physical Optics Propagation (POP) fiber coupling
- Ansys Lumerical FDTD (via lumapi)
- ZOS-API
- Python
- NSC ray tracing
- Optical modeling
- Detector irradiance analysis
- Monte Carlo ray tracing
- Optical optimization and validation
