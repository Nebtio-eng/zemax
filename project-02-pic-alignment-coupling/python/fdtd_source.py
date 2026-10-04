"""Stage 8b: the grating-coupler beam from Lumerical FDTD, exported as a Zemax POP .zbf.

Physics boundary: FDTD computes the field the grating actually emits (period,
etch, fill factor, length, width and directionality all included). Zemax then
propagates that field through the substrate, lens and air gap exactly as in B.

Geometry (all grating parameters ASSUMED; the primary paper does not publish
them; see literature/extracted_parameters.csv):
  220 nm Si device layer on 2 um BOX on a Si substrate, SiO2 top cladding,
  uniform shallow-etch grating (70 nm, fill factor 0.5), period chosen for
  ~10 deg air-equivalent emission, 20 periods, 12 um wide, TE input.
The beam is recorded on a z-normal monitor inside the Si substrate. Its
propagation angle comes from the angular spectrum, and the field is projected
(farfieldexact) onto a plane normal to that direction, then embedded in a
POP-sized grid and written as .zbf.

No OpticStudio here: this runs in its own process with lumapi only.

Usage:  python fdtd_source.py [--coarse] [--periods 20]
"""
import argparse, json, math, os, struct, sys, time
from pathlib import Path

import numpy as np

LUMAPI_DIR = r"C:\Program Files\Lumerical\v261\api\python"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "fdtd"
BEAMFILES = Path(os.path.expanduser("~")) / "Documents" / "Zemax" / "POP" / "BEAMFILES"
PROGRESS_LOG = ROOT / "results" / "progress.log"

P = dict(
    wavelength_um=1.31, t_si_um=0.22, etch_um=0.07, box_um=2.0, fill_factor=0.5,
    period_um=0.4895, n_periods=20, width_um=12.0, clad_top_um=1.0, sub_monitor_depth_um=0.4,
    input_length_um=2.0, si_material="Si (Silicon) - Palik", ox_material="SiO2 (Glass) - Palik",
    mesh_accuracy=2, sim_time_fs=1500.0,
    # B's pixel (0.4 mm / 1024) on a 0.8 mm window: the real beam has 4.3% of its power at 11-29 deg in Si, which
    # reaches ~190-345 um off-axis at the lens; B's 0.4 mm window left 0.41% at its edge (Stage 8b)
    zbf_n=2048, zbf_dx_um=0.4 * 1000 / 1024, zbf_index=3.5039127043661025,
)


def progress(msg):
    print(msg, flush=True)
    with open(PROGRESS_LOG, "a") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + "[stage8b] " + msg + "\n")


def lumapi():
    sys.path.insert(0, LUMAPI_DIR)
    import lumapi as L
    return L


# ===================================================================== model build
def build(fd, p):
    """Grating coupler in FDTD. x along the grating (light travels +x), y across, z up."""
    um = 1e-6
    L = p["n_periods"] * p["period_um"]
    t, e, box = p["t_si_um"], p["etch_um"], p["box_um"]
    z_top, z_bot = t / 2, -t / 2
    x0, x1 = -p["input_length_um"], L + 2.0
    W = p["width_um"]
    ymax = W / 2 + 2.0
    zmin = z_bot - box - p["sub_monitor_depth_um"] - 0.6
    zmax = z_top + p["clad_top_um"]
    fd.switchtolayout()
    fd.deleteall()
    # cladding and BOX: one big oxide block; substrate below; device layer; trenches
    fd.addrect(name="oxide", x_min=(x0 - 5) * um, x_max=(x1 + 5) * um, y_min=-(ymax + 5) * um, y_max=(ymax + 5) * um,
               z_min=(z_bot - box) * um, z_max=(zmax + 5) * um, material=p["ox_material"])
    fd.addrect(name="substrate", x_min=(x0 - 5) * um, x_max=(x1 + 5) * um, y_min=-(ymax + 5) * um, y_max=(ymax + 5) * um,
               z_min=(zmin - 5) * um, z_max=(z_bot - box) * um, material=p["si_material"])
    fd.addrect(name="device", x_min=(x0 - 5) * um, x_max=(L + 0.5) * um, y_min=-W / 2 * um, y_max=W / 2 * um,
               z_min=z_bot * um, z_max=z_top * um, material=p["si_material"])
    w_trench = (1 - p["fill_factor"]) * p["period_um"]
    for i in range(p["n_periods"]):
        xs = i * p["period_um"] + p["fill_factor"] * p["period_um"]
        fd.addrect(name="trench%d" % i, x_min=xs * um, x_max=(xs + w_trench) * um, y_min=-W / 2 * um, y_max=W / 2 * um,
                   z_min=(z_top - e) * um, z_max=z_top * um, material=p["ox_material"])
    # simulation region
    fd.addfdtd(dimension="3D", x_min=x0 * um, x_max=x1 * um, y_min=-ymax * um, y_max=ymax * um,
               z_min=zmin * um, z_max=zmax * um, mesh_accuracy=p["mesh_accuracy"],
               simulation_time=p["sim_time_fs"] * 1e-15)
    # fine mesh through the device layer and the grating
    fd.addmesh(name="grating_mesh", x_min=-0.5 * um, x_max=(L + 0.5) * um, y_min=-ymax * um, y_max=ymax * um,
               z_min=(z_bot - 0.02) * um, z_max=(z_top + 0.02) * um, override_x_mesh=1, override_y_mesh=0,
               override_z_mesh=1, dx=0.02 * um, dz=0.01 * um)
    # TE fundamental mode injected along +x in the wide input section
    fd.addmode(name="source", injection_axis="x-axis", direction="Forward", x=(x0 + 0.5) * um,
               y=0, y_span=(W + 3.0) * um, z=0, z_span=(t + 2.0) * um,
               center_wavelength=p["wavelength_um"] * um, wavelength_span=0)
    fd.set("mode selection", "fundamental TE mode")
    # field monitor inside the Si substrate, below the BOX
    zmon = z_bot - box - p["sub_monitor_depth_um"]
    fd.addpower(name="near", monitor_type="2D Z-normal", x_min=x0 * um, x_max=x1 * um, y_min=-ymax * um,
                y_max=ymax * um, z=zmon * um)
    fd.setglobalmonitor("frequency points", 1)
    fd.setglobalsource("wavelength start", p["wavelength_um"] * um)
    fd.setglobalsource("wavelength stop", p["wavelength_um"] * um)
    return dict(grating_length_um=L, z_monitor_um=zmon, x_range_um=(x0, x1), y_half_um=ymax)


# =================================================================== field analysis
def near_field(fd):
    E = fd.getresult("near", "E")
    x, y, z = (np.asarray(E[k]).ravel() / 1e-6 for k in ("x", "y", "z"))
    F = np.squeeze(np.asarray(E["E"]))            # (nx, ny, 3) for a 2D monitor at one frequency
    return x, y, float(z[0]), F


def dominant_component(F):
    p = [float(np.sum(abs(F[..., i]) ** 2)) for i in range(3)]
    return int(np.argmax(p)), [v / sum(p) for v in p]


def propagation_angle(x, E, n, lam):
    """Beam direction in the substrate for a field E(x, y). The monitor grid is NOT uniform (FDTD mesh
    override), so everything is done on a uniform resample.
    Returns (mean angle, peak angle) in degrees: mean = power-weighted average kx (the direction the
    beam centroid travels: used for the plane and the Zemax field angle); peak = angular-spectrum
    maximum (what the grating equation predicts)."""
    xu = np.linspace(x[0], x[-1], 4 * len(x))
    Eu = np.array([np.interp(xu, x, c.real) + 1j * np.interp(xu, x, c.imag) for c in E.T]).T
    pad = 16 * len(xu)
    spec = sum(abs(np.fft.fftshift(np.fft.fft(Eu[:, j], pad))) ** 2 for j in range(Eu.shape[1]))
    kx = np.fft.fftshift(np.fft.fftfreq(pad, d=xu[1] - xu[0])) * 2 * math.pi
    k = 2 * math.pi * n / lam
    prop = abs(kx) < k
    kmean = (spec[prop] * kx[prop]).sum() / spec[prop].sum()
    return math.degrees(math.asin(kmean / k)), math.degrees(math.asin(kx[np.argmax(spec)] / k))


def second_moment_radii(u, v, I):
    """1/e^2-equivalent radii (2 sigma) along u and v of an intensity map I[u, v]."""
    I = I / I.sum()
    U, V = np.meshgrid(u, v, indexing="ij")
    cu, cv = (I * U).sum(), (I * V).sum()
    return 2 * math.sqrt((I * (U - cu) ** 2).sum()), 2 * math.sqrt((I * (V - cv) ** 2).sum()), cu, cv


# ======================================================================== ZBF output
def launch_x_factor(angle_si_deg, index=P["zbf_index"]):
    """POP launches a file beam on the AIR side of surface 1 and refracts it into the silicon at the field angle,
    which stretches x by cos(theta_Si)/cos(theta_air) (Stage 8b: window 0.4 -> 0.40714 mm at 11.20 deg, exactly this
    ratio; a Gaussian became 1.7% narrower in x at the lens). The FDTD field is already in silicon, normal to the
    beam, so its x pitch is written pre-shrunk by the inverse factor: POP's launch then restores the true size."""
    t = math.radians(angle_si_deg)
    a = math.asin(index * math.sin(t))
    return math.cos(a) / math.cos(t)


def write_zbf(path, E, dx_um, lam_um, waist_x_um, waist_y_um, index, x_factor=1.0):
    """OpticStudio ZBF (format 1, unpolarised): 9 ints, 20 doubles, then complex E[y][x]. Units mm.

    E is in Lumerical's phase convention (e^{+ikz}: a phase exp(+i k sin(a) x) tilts the beam towards +x).
    OpticStudio reads the opposite sign (Stage 8b convention control: a +1 deg tilt sent the best fibre
    to -11.10 um instead of +11.09 um, while a +3 um amplitude decentre went to +2.604 um as predicted).
    So the field is conjugated on writing, and read_zbf conjugates back.
    x_factor: written x pitch = dx_um * x_factor (launch_x_factor; y pitch is always the true pitch)."""
    E = np.conj(E)
    ny, nx = E.shape
    mm = 1e-3
    zx = math.pi * (waist_x_um * mm) ** 2 * index / (lam_um * mm)
    zy = math.pi * (waist_y_um * mm) ** 2 * index / (lam_um * mm)
    hdr = struct.pack("<9i", 1, nx, ny, 0, 0, 1, 2, 0, 0)
    hdr += struct.pack("<20d", dx_um * x_factor * mm, dx_um * mm, 0.0, zx, waist_x_um * mm, 0.0, zy, waist_y_um * mm,
                       lam_um * mm, index, 0.0, 0.0, *([0.0] * 8))
    d = np.empty(nx * ny * 2)
    d[0::2], d[1::2] = E.real.ravel(), E.imag.ravel()
    Path(path).write_bytes(hdr + d.astype("<f8").tobytes())


def read_zbf(path):
    b = Path(path).read_bytes()
    _, nx, ny, ispol, units = struct.unpack_from("<5i", b, 0)
    hd = struct.unpack_from("<20d", b, 36)
    d = np.frombuffer(b, "<f8", offset=196, count=nx * ny * 2)
    return dict(nx=nx, ny=ny, dx_mm=hd[0], dy_mm=hd[1], E=np.conj((d[0::2] + 1j * d[1::2]).reshape(ny, nx)))  # Lumerical convention


def gaussian_zbf(path, waist_um=4.6):
    """The analytic B source as a ZBF: the pipeline control (must reproduce B exactly)."""
    n, dx = P["zbf_n"], P["zbf_dx_um"]
    c = (np.arange(n) - n // 2) * dx
    X, Y = np.meshgrid(c, c)
    write_zbf(path, np.exp(-(X ** 2 + Y ** 2) / waist_um ** 2).astype(complex), dx, P["wavelength_um"],
              waist_um, waist_um, P["zbf_index"])


# ============================================================================= main
def geometry(p):
    """The numbers build() returns, without building (for re-analysing a saved run)."""
    L = p["n_periods"] * p["period_um"]
    return dict(grating_length_um=L, z_monitor_um=-p["t_si_um"] / 2 - p["box_um"] - p["sub_monitor_depth_um"],
                x_range_um=(-p["input_length_um"], L + 2.0), y_half_um=p["width_um"] / 2 + 2.0)


def run(p, coarse=False, tag=None, analyse_only=False):
    OUT.mkdir(parents=True, exist_ok=True)
    L = lumapi()
    p = dict(p)
    if coarse:
        p.update(mesh_accuracy=1, sim_time_fs=1000.0)
    tag = tag or ("grating_%dp%s" % (p["n_periods"], "_coarse" if coarse else ""))
    t0 = time.time()
    fsp = OUT / (tag + ".fsp")
    fd = L.FDTD(filename=str(fsp), hide=True) if analyse_only else L.FDTD(hide=True)
    try:
        geo = build(fd, p) if not analyse_only else geometry(p)
        if not analyse_only:
            fd.save(str(fsp))
            progress("%s: built (grating %.2f um); running" % (tag, geo["grating_length_um"]))
            fd.run()
            progress("%s: FDTD finished in %.0f s" % (tag, time.time() - t0))
        x, y, zmon, F = near_field(fd)
        n_si = float(np.real(np.ravel(fd.getindex(p["si_material"], 299792458.0 / (p["wavelength_um"] * 1e-6)))[0]))
        res = dict(tag=tag, params=p, geometry=geo, n_si_fdtd=n_si, wall_time_s=None)
        comp, frac = dominant_component(F)
        Ec = F[..., comp]                              # (nx, ny)
        th, th_peak = propagation_angle(x, Ec, n_si, p["wavelength_um"])
        wx, wy, cx, cy = second_moment_radii(x, y, abs(Ec) ** 2 * np.gradient(x)[:, None])   # non-uniform x
        progress("%s: dominant E%s (%.3f of |E|^2); angle in Si: mean %.3f deg, peak %.3f deg; "
                 "near-field 2-sigma radii x %.3f y %.3f um" % (tag, "xyz"[comp], frac[comp], th, th_peak, wx, wy))
        np.save(OUT / (tag + "_near.npy"), dict(x=x, y=y, E=Ec), allow_pickle=True)
        res.update(component="xyz"[comp], power_fractions=frac, angle_si_deg=th, angle_si_peak_deg=th_peak,
                   angle_air_equivalent_deg=math.degrees(math.asin(p["zbf_index"] * math.sin(math.radians(th)))),
                   near_radii_um=dict(x=wx, y=wy), centroid_um=dict(x=cx, y=cy), zbf_variants={})
        # the POP grid of B, plus the two grids of the project convergence test (grid x2; grid x2 + window x2)
        dx0 = p["zbf_dx_um"]
        for suffix, n, dx in (("", p["zbf_n"], dx0), ("_g2", 2 * p["zbf_n"], dx0 / 2), ("_g2w2", 2 * p["zbf_n"], dx0)):
            proj = project(fd, th, cx, cy, zmon, p, comp, dx)
            embed_and_write(proj, p, tag, res, n, dx, suffix)
        res["projection"] = proj["meta"]
        res["wall_time_s"] = time.time() - t0
        (OUT / (tag + ".json")).write_text(json.dumps(res, indent=2, default=float))
        return res
    finally:
        fd.close()


def project(fd, theta_deg, cx, cy, zmon, p, comp, dx, half_um=16.0, d_um=1.0):
    """farfieldexact onto the plane normal to the beam direction, d_um past the monitor."""
    um = 1e-6
    th = math.radians(theta_deg)
    u_dir = np.array([math.sin(th), 0.0, -math.cos(th)])          # propagation: down, tilted in x
    e1 = np.array([math.cos(th), 0.0, math.sin(th)])              # in-plane axis normal to u in x-z
    e2 = np.array([0.0, 1.0, 0.0])
    c = np.array([cx, cy, zmon]) + d_um * u_dir
    s = np.arange(-half_um, half_um + dx / 2, dx)
    S1, S2 = np.meshgrid(s, s, indexing="ij")
    pts = c[None, None, :] + S1[..., None] * e1 + S2[..., None] * e2
    X, Y, Z = (pts[..., i].ravel() * um for i in range(3))
    # farfieldexact takes a point list and returns (N, 3); farfieldexact3d treats x, y, z as grid axes
    E = np.asarray(fd.farfieldexact("near", X, Y, Z)).reshape(len(s), len(s), 3)
    # field component transverse to the rotated plane: project onto e1 (x-like) and e2 (y)
    Ee1 = E[..., 0] * e1[0] + E[..., 2] * e1[2]
    Ee2 = E[..., 1]
    use = Ee2 if comp == 1 else Ee1
    return dict(s=s, field=use, other=Ee1 if comp == 1 else Ee2,
                meta=dict(plane_centre_um=c.tolist(), normal=u_dir.tolist(), distance_um=d_um, half_width_um=half_um))


def embed_and_write(proj, p, tag, res, n, dx, suffix=""):
    """Embed the projected field in an n x n grid of pitch dx, centre it, write P02_<tag><suffix>.ZBF.
    Only the base grid (suffix "") updates the beam measurements in res."""
    s, F = proj["s"], proj["field"]                                # F[i along e1 (x'), j along y]
    G = np.zeros((n, n), complex)
    k = len(s)
    o = n // 2 - k // 2
    G[o:o + k, o:o + k] = F.T                                      # ZBF rows are y, columns x
    I = abs(G) ** 2
    c = (np.arange(n) - n // 2) * dx
    wx, wy, cx, cy = second_moment_radii(c, c, I.T)
    # check: on a plane normal to the mean direction the power-weighted phase tilt must be ~0
    kw = 2 * math.pi * p["zbf_index"] / p["wavelength_um"]
    tilt = []
    for ax in (1, 0):                                              # x' (columns), then y (rows)
        g = np.gradient(G, dx, axis=ax)
        tilt.append(math.degrees(math.asin(float(np.imag(np.conj(G) * g).sum() / I.sum()) / kw)))
    # centre the beam centroid on the grid (the chief ray); the field itself is not altered
    shift_x, shift_y = int(round(cx / dx)), int(round(cy / dx))
    G = np.roll(np.roll(G, -shift_y, axis=0), -shift_x, axis=1)
    name = "P02_%s%s.ZBF" % (tag, suffix)
    xf = launch_x_factor(res["angle_si_deg"], p["zbf_index"])
    write_zbf(BEAMFILES / name, G, dx, p["wavelength_um"], wx, wy, p["zbf_index"], x_factor=xf)
    res["zbf_variants"][suffix or "base"] = dict(file=name, n=n, dx_um=dx, x_pitch_written_um=dx * xf, width_mm=n * dx / 1000,
                                                 radii_um=[wx, wy], residual_tilt_deg=tilt)
    if suffix:
        return
    np.save(OUT / (tag + "_zbf_field.npy"), G[o - 8:o + k + 8, o - 8:o + k + 8])
    res.update(residual_tilt_deg=dict(x=tilt[0], y=tilt[1]))
    res.update(zbf=str(BEAMFILES / name), rotated_radii_um=dict(x=wx, y=wy),
               ellipticity_wx_over_wy=wx / wy, other_component_power_fraction=float(
                   (abs(proj["other"]) ** 2).sum() / ((abs(proj["other"]) ** 2).sum() + (abs(F) ** 2).sum())))
    progress("%s: rotated-plane 2-sigma radii x' %.3f um, y %.3f um (ratio %.3f); residual tilt x' %.3f y %.3f deg; "
             "ZBF written" % (tag, wx, wy, wx / wy, tilt[0], tilt[1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coarse", action="store_true")
    ap.add_argument("--periods", type=int, default=P["n_periods"])
    ap.add_argument("--analyse", action="store_true", help="re-analyse the saved .fsp without re-running FDTD")
    a = ap.parse_args()
    p = dict(P, n_periods=a.periods)
    print(json.dumps(run(p, coarse=a.coarse, analyse_only=a.analyse), indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
