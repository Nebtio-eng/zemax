# ZOS-API and POP gotchas found in this project

Every entry was found the hard way on OpticStudio 2026 R1.00 (build 260127,
Enterprise) through pythonnet 3.2.0 and Python 3.13 (64-bit). The stage where
each one surfaced is in brackets. Where code works around it, the place is named.

## Connection and session

1. **`ZOSAPI_NetHelper.dll` sits at the install root** in 2026 R1, not under
   `ZOS-API\Libraries\`. Probe both. [Stage 0] `coupling_analysis.load_zosapi`
2. **There is no `app.Edition`.** `app.LicenseStatus` *is* the edition
   (`EnterpriseEdition`, ...). [Stage 0]
3. **`ConnectAsExtension(0)` returns a stub** (`IsAlive = False`,
   `PrimarySystem = None`) unless Programming > Interactive Extension is armed.
   [Stage 3]
4. **Arming is consumed by the first connection.** A second script needs a
   fresh arm. Keep one process connected for a whole session. [Stage 3]
5. **An OpticStudio restart silently kills an extension connection**: calls fail
   with `RemotingException: Failed to connect to an IPC Port`. [Stage 7]
6. **Never run two standalone instances at once**: they collide over the licence
   seat (`FRU__delta_init()` errors). [Stage 0]
7. **Analysis objects go stale** after `LoadFile`/`SaveAs`: the old handle
   reports another window's title and `SaveTo` returns False. Create a fresh
   analysis. [Stage 5]

## Units, materials, model

8. **Lens units cannot be micrometres** (mm, cm, in, m only). Work in mm:
   4.6 um = 0.0046. [Stage 3]
9. **Out-of-range glass silently becomes n = 1.0.** SILICON (INFRARED.AGF) is
   defined for 1.36-11 um only. At 1.31 um, INDX returned 1.0 and POP failed
   with no message. Fix: project catalogue `zemax/glasscat/PROJECT02.AGF`
   (SILICON_1310, lower limit 1.30 um), copied into `app.GlassDir` by
   `build_model`. [Stage 7]
10. **Custom catalogues must live in `app.GlassDir`**
    (`Documents\Zemax\GLASSCAT`); then `MaterialCatalogs.AddCatalog(name)`.
    [Stage 7]
11. **Surface comments are truncated to 32 characters.** [Stage 3]
12. **Refractive index:** `MFE.GetOperandValue(MeritOperandType.INDX, surf, wave,
    0, 0, 0, 0, 0, 0)`. [Stage 7]
13. **Aperture:** `ApertureData.CreateApertureTypeSettings(SurfaceApertureTypes.
    CircularAperture)`, set `._S_CircularAperture.MaximumRadius`, then
    `ChangeApertureTypeSettings`. Remove with `getattr(SurfaceApertureTypes,
    "None")` (`None` is a Python keyword). [Stage 9]
14. **Coordinate-break tilt about X** is `GetSurfaceCell(SurfaceColumn.Par3)`
    (Par1/2 decenter X/Y, Par3/4/5 tilt X/Y/Z). [Stage 9]

## POP settings

15. **`ModifySettings(cfg, token, value)` returns True for any token**, real or
    invented. It cannot be used to validate a setting. [Stage 3]
16. **The real token names are UTF-16 strings in `ZemaxCore.dll`**, found by
    scanning the binary: `POP_START`, `POP_END`, `POP_WAVE`, `POP_FIELD`,
    `POP_BEAMTYPE` (0 = Gaussian Waist), `POP_PARAM1..8` (waist X/Y, decenter
    X/Y, aperture X/Y, order X/Y), `POP_SURFTOBEAM`, `POP_SAMPX/Y`, `POP_WIDEX/Y`,
    `POP_COMPUTE`, `POP_FIBERTYPE` (0 = Gaussian), `POP_FPARAM1..8` (same layout
    for the receiver: FPARAM3/4 = X/Y decenter), `POP_TILTX/Y` (receiver tilt,
    degrees). [Stage 3, 6]
17. **`POP_SAMPX/Y` is an index:** grid = 2^(index + 4) (5 -> 512, 6 -> 1024,
    7 -> 2048). [Stage 3]
18. **The .CFG is binary and there is no getter.** Read-back works by locating
    each token's byte offset (write two marker values, diff the files): integers
    as int32, doubles as float64. `PopSession._locate`, `verify`. [Stage 5]
19. **A write occasionally does not land** (once in ~10^4: the CFG kept the old
    value). Retry once, then verify strictly. `PopSession.set`. [Stage 9]
20. **OpticStudio's .NET parser is not always correctly rounded**: a decimal
    string can come back one ulp away. Doubles are verified to 2 ulp, integers
    exactly. [Stage 9]
21. **`repr(numpy.float64(x))` is `np.float64(x)`**, which OpticStudio cannot
    parse; it silently stores 0. Convert to a Python float first. [Stage 9]
22. **POPD reads the *saved default* POP settings**, not the open analysis.
    Before any POPD evaluation: `settings.LoadFrom(cfg)` then `settings.Save()`.
    Without this POPD returns 0. `PopSession.commit`. [Stage 3]
23. **POPD columns:** Param1 = surface, Param2 = wavelength, Param3 = field,
    Param4 = Data (0 total, 1 system S, 2 receiver T). [Stage 3]
24. **`GetResults().GetTextFile(path)` writes UTF-16** with a BOM. [Stage 3]
25. **The POP Prop Report tab is not exposed by the API.** `NumberOfMessages` was
    0 everywhere; the tab itself must be read in the GUI. [Stage 3, 5]

## POP physics behaviour

26. **The default propagator silently rescales its grid beyond the Rayleigh
    range** (A0: window 0.08 -> 0.84 mm at a 100 um gap; A1/B: 0.40 -> 0.67 /
    0.60 mm). eta is unchanged to 1e-12, but the stated window no longer
    describes the run. Use `UseAngularSpectrumPropagator = True`. [Stage 5, 7]
27. **Even the angular-spectrum window moves slightly after a curved surface**
    (0.4 -> 0.39999 mm in B). The report check allows 0.1%. [Stage 7]
28. **Resample After Refraction without Auto destroys the field**: it re-grids
    onto the surface's default 32 x 32, 1 mm grid (A1: eta 0.9996 instead of
    0.2369). With Auto on it agrees to 0.03%. Left off. [Stage 7]
29. **POP defines the fibre relative to the beam's chief ray.** With the lens
    tilted 4 deg, the optimal POP fibre tilt was -0.003 deg although the beam
    physically leaves at 14.2 deg. Incidence studies measure the re-pointed case;
    the physical fibre tilt is the Snell angle. [Stage 9]
30. **POP fibre tilt enters as tan(theta)**: 2 deg gives 0.862084 against
    0.862188 (small-angle) and 0.862240 (sin). Negligible below 2 deg. [Stage 6]
31. **POP system efficiency S excludes Fresnel reflection** (S = 1.000 in every
    configuration). An uncoated Si-air face would lose 30.9% (1.61 dB); the
    experiment expected 1.85 dB and recovered 2 dB with a 170 nm SiN coating
    (Mangal et al. 2021). [Stage 8]
32. **Grid-doubling convergence is weakly discriminating for ideal Gaussians**
    (0.0000% in A0). It needs a lens, an aperture or long propagation to bite;
    deliberately bad settings (64-point grid, 40 um window) were used to show the
    test can fail. [Stage 3, 7]

## File beams (.zbf) and Lumerical [Stage 8b]

33. **Lumerical and OpticStudio use opposite phase signs.** Lumerical fields are
    e^{+ikz}; OpticStudio reads .zbf phase with the opposite sign. A Gaussian
    written with a +1 deg tilt put the best fibre at -11.102 um instead of the
    predicted +11.091 um; after conjugating, +11.102 um. A pure decentre
    (+3 um) is blind to the sign (2.603 vs 2.604 um), so test with a tilt.
    `fdtd_source.write_zbf` conjugates on write, `read_zbf` conjugates back.
34. **POP launches a file beam on the air side of surface 1** and refracts it in
    at the field angle. At 11.202 deg in air (3.178 deg in Si) the x axis is
    stretched by cos(theta_Si)/cos(theta_air) = 1.01785: the launch window
    became 0.407140 mm instead of 0.4 mm. The x pitch is written pre-shrunk by
    that factor (`fdtd_source.launch_x_factor`); a compensated Gaussian then
    arrives round at the lens (X/Y 0.9983) with the window exactly as set. The
    remaining 0.17% (0.40644 mm at the lens) is real oblique refraction at the
    curved face, not an error.
35. **B's 0.4 mm window does not hold a real grating beam.** 4.3% of its power
    leaves at 11-29 deg in Si (side lobes), and 0.41% reached the window edge
    at the lens. B_fdtd uses 2048 points on 0.8 mm (same pixel as B). A file
    source fixes the pixel itself, so the grid-doubling convergence test must
    re-write the .zbf at the new pitch (`_g2`, `_g2w2` variants), or it tests
    nothing.
36. **Lumerical: save before `run`.** Otherwise a save dialog appears; cancelling
    it returns no data, without an error.
37. **Lumerical: monitor data are on the non-uniform simulation mesh** under a
    mesh override. The first beam-angle estimate, made before this was
    noticed, was 0.90 deg; the correct mean angle is 3.2 deg. Use the
    returned x, y vectors and weight sums by the cell width
    (`np.gradient(x)`, `fdtd_source`).
38. **Lumerical: `farfieldexact` takes a list of points; `farfieldexact3d`
    treats its inputs as grid axes.** Same arguments, different meaning.
39. **Lumerical: `getindex` returns a 2-D array** even for a 1-D query; index it
    explicitly.
