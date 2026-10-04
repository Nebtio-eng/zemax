"""Stage 8b, OpticStudio half: the FDTD grating beam as the POP source of configuration B.

B_fdtd is an otherwise untouched copy of B: same surfaces, receiver, window and
propagators. Two changes only: the POP source is the .zbf written by
fdtd_source.py (POP_BEAMTYPE 4 = File), and field 1 is set to the beam's
air-equivalent angle so the chief ray leaves the grating at the FDTD angle in Si.

Every Zemax number is predicted first, without OpticStudio, by reciprocity: the
receiver's fibre mode is carried back to the source plane with Gaussian ABCD
optics (fibre decenter and tilt become source-plane offset and tilt) and
overlapped with the FDTD field. The predictor must reproduce B for B's own
Gaussian source before it is used.

Controls: an elliptical Gaussian .zbf (waist 3.5 x 4.6 um) shows the pipeline
resolves X/Y asymmetry, and fixes which POP tilt token is which axis.

Usage (inside a connected session):  import fdtd_coupling as FC; FC.run_all(app, zos, "grating_20p")
"""
import copy, csv, json, math, time
from pathlib import Path

import numpy as np

import analytic as AN
import coupling_analysis as C
import fdtd_source as FS
from coupling_analysis import ROOT, PopSession

B_CONFIG = ROOT / "zemax" / "microlens" / "B" / "run_config.json"
OUT = ROOT / "results" / "fdtd"
B_REF = dict(lateral_rel=8.1399, angular_rel=0.6758, longitudinal_rel=713.3)   # Stage 7, Gaussian source
TILT = {"x": "POP_TILTY", "y": "POP_TILTX"}     # receiver tilt in the x-z plane is a rotation about Y (control 1)


def log(m):
    C.progress("stage8b", m)


# ============================================================ independent prediction
def system_matrix(gap_um, t_um=630.0, radius_um=-480.0, n=AN.N_SI):
    """Source plane -> fibre, reduced coordinates (x, n*theta): Si, lens surface, air gap."""
    return (np.array([[1.0, gap_um], [0.0, 1.0]]) @ np.array([[1.0, 0.0], [(n - 1) / radius_um, 1.0]])
            @ np.array([[1.0, t_um / n], [0.0, 1.0]]))


def fibre_mode_at_source(dx=0.0, dy=0.0, tx=0.0, ty=0.0, gap=20.0, wf=17.0, n=AN.N_SI, lam=AN.LAM):
    """The receiver mode (decentred dx, dy um; tilted tx, ty deg, POP tan convention) carried back to the
    source plane: Gaussian radius w, wavefront R, centre (xg, yg), direction sines (sx, sy) in Si."""
    Mi = np.linalg.inv(system_matrix(gap, n=n))
    qf = 1j * math.pi * wf ** 2 / lam
    q = n * (Mi[0, 0] * qf + Mi[0, 1]) / (Mi[1, 0] * qf + Mi[1, 1])
    w, R = AN.w_R(q, n, lam)
    xg, ux = Mi @ [dx, math.tan(math.radians(tx))]
    yg, uy = Mi @ [dy, math.tan(math.radians(ty))]
    return w, R, xg, yg, ux / n, uy / n


def overlap_eta(E, c, mode, n=AN.N_SI, lam=AN.LAM):
    """Power overlap of a source-plane field E[y, x] (grid c, um; Lumerical phase convention e^{+ikz})
    with the back-propagated receiver mode."""
    w, R, xg, yg, sx, sy = mode
    X, Y = np.meshgrid(c, c)
    k = 2 * math.pi * n / lam
    r2 = (X - xg) ** 2 + (Y - yg) ** 2
    Eb = np.exp(-r2 / w ** 2 + (0 if math.isinf(R) else 1j * k * r2 / (2 * R)) + 1j * k * (sx * X + sy * Y))
    return float(abs(np.sum(E * np.conj(Eb))) ** 2 / (np.sum(abs(E) ** 2) * np.sum(abs(Eb) ** 2)))


def predicted_loss(E, c, **kw):
    return float(AN.loss_db(overlap_eta(E, c, fibre_mode_at_source(**kw))))


def predicted_optimum(E, c, start=(0.0, 0.0)):
    """Fibre X decenter and X-Z tilt that maximise coupling (gap held at its nominal 20 um)."""
    from scipy.optimize import minimize
    r = minimize(lambda v: predicted_loss(E, c, dx=v[0], tx=v[1]), list(start), method="Nelder-Mead",
                 options=dict(xatol=1e-4, fatol=1e-7))
    return dict(dx=float(r.x[0]), tx=float(r.x[1]))


def predict(E, c, base=None):
    """Loss and 1-dB tolerances (rel and abs, each side, each axis) of a source field in B, about the fibre
    position `base` (default: on the chief ray)."""
    base = dict(base or {})
    L = lambda **kw: predicted_loss(E, c, **{**base, **kw})
    l0 = L()
    out = dict(loss_dB=l0, base=base)
    axes = (("lateral_x", "dx", 6.0, 80.0), ("lateral_y", "dy", 6.0, 80.0),
            ("angular_x", "tx", 0.5, 10.0), ("angular_y", "ty", 0.5, 10.0))
    for conv, level in (("rel", l0 + 1.0), ("abs", 1.0)):
        for name, key, g0, hmax in axes:
            b0 = base.get(key, 0.0)
            sides = [AN.level_crossing(lambda v: L(**{key: b0 + s * v}), level, g0, hmax, l0) if level > l0 else float("nan")
                     for s in (1, -1)]
            out["%s_%s_sides" % (name, conv)] = sides
            out["%s_%s" % (name, conv)] = float(np.mean(sides))
        out["longitudinal_" + conv] = AN.level_crossing(lambda dz: L(gap=20.0 + dz), level, 300.0, 8000.0, l0) \
            if level > l0 else float("nan")
        for ax in "xy":
            out["product_%s_%s" % (ax, conv)] = out["lateral_%s_%s" % (ax, conv)] * math.radians(out["angular_%s_%s" % (ax, conv)])
    return out


def source_grid(n=1024, dx=FS.P["zbf_dx_um"], half=None):
    c = (np.arange(n) - n // 2) * dx
    return c if half is None else c[abs(c) <= half]


def gaussian_field(c, wx, wy):
    X, Y = np.meshgrid(c, c)
    return np.exp(-X ** 2 / wx ** 2 - Y ** 2 / wy ** 2).astype(complex)


def fdtd_field(tag):
    """The written .zbf field (base grid), cropped to where it lives. The y pitch is the true pitch (the x pitch is
    written pre-shrunk for POP's launch, fdtd_source.launch_x_factor)."""
    z = FS.read_zbf(FS.BEAMFILES / ("P02_%s.ZBF" % tag))
    c = (np.arange(z["nx"]) - z["nx"] // 2) * z["dy_mm"] * 1000
    keep = abs(c) <= 30.0
    return z["E"][np.ix_(keep, keep)], c[keep]


def predictor_selftest():
    """B's own Gaussian source must give B (Stage 7) through the reciprocity route."""
    c = source_grid(half=30.0)
    p = predict(gaussian_field(c, AN.W0, AN.W0), c)
    for k, v in B_REF.items():
        got = p[k.replace("lateral", "lateral_x").replace("angular", "angular_x")]
        if abs(got / v - 1) > 2e-3:
            raise RuntimeError("predictor self-test: %s %.4f vs B %.4f" % (k, got, v))
    return p


# ===================================================================== OpticStudio
def make_config(app, zos, tag, angle_air_deg, src):
    """Write zemax/microlens/B_fdtd/<tag>/run_config.json and the model: B plus field angle plus file source."""
    b = json.loads(B_CONFIG.read_text())
    c = copy.deepcopy(b)
    c["id"] = "B_fdtd_" + tag
    c["configuration"] = "B with the FDTD grating beam (%s) as POP source; field angle %.4f deg (air)" % (tag, angle_air_deg)
    c["stage"] = "8b"
    c["model"] = "zemax/microlens/B_fdtd/%s/B_fdtd.zmx" % tag
    c["pop_tokens"]["POP_BEAMTYPE"] = 4
    # numerical settings, not optics: B's 0.4 mm window does not hold the real beam's 11-29 deg skirt (0.41% at the
    # edge); same 0.39 um pixel on twice the window
    n = FS.P["zbf_n"]
    c["pop_tokens"].update(POP_SAMPX=int(round(math.log2(n))) - 4, POP_SAMPY=int(round(math.log2(n))) - 4,
                           POP_WIDEX=n * FS.P["zbf_dx_um"] / 1000, POP_WIDEY=n * FS.P["zbf_dx_um"] / 1000)
    c["pop_strings"] = {"_comment": "string tokens: not byte-verifiable; verified by eta", "POP_SOURCEFILE": "P02_%s.ZBF" % tag}
    c["field_x_deg"] = angle_air_deg
    c["source"] = src
    c.pop("results_nominal", None)
    sysm = app.PrimarySystem
    sysm.LoadFile(str(ROOT / b["model"]), False)
    C.check_model(sysm, c)
    f = sysm.SystemData.Fields
    if "Angle" not in str(f.GetFieldType()):
        raise C.SettingsMismatch("B field type is %s, expected Angle" % f.GetFieldType())
    if f.NumberOfFields != 1:
        raise C.SettingsMismatch("B has %d fields" % f.NumberOfFields)
    f.GetField(1).X = angle_air_deg
    f.GetField(1).Y = 0.0
    if abs(f.GetField(1).X - angle_air_deg) > 1e-12:
        raise C.SettingsMismatch("field angle did not read back")
    path = ROOT / c["model"]
    path.parent.mkdir(parents=True, exist_ok=True)
    sysm.SaveAs(str(path))
    cfg_path = path.parent / "run_config.json"
    cfg_path.write_text(json.dumps(c, indent=2))
    return c, cfg_path


class Bench8b:
    def __init__(self, app, zos, cfgd):
        self.app, self.zos, self.cfgd = app, zos, cfgd
        self.ps = PopSession(app, cfgd, zos)
        self.gap0 = cfgd["parameters"]["gap_um"]
        self.source = cfgd["pop_strings"]["POP_SOURCEFILE"]

    def set_source(self, name):
        self.ps.st.ModifySettings(self.ps.cfg, "POP_SOURCEFILE", name)
        self.source = name

    def popd(self, dx=0.0, dy=0.0, tx=0.0, ty=0.0, gap=None):
        ps = self.ps
        ps.set("POP_FPARAM3", dx / 1000)
        ps.set("POP_FPARAM4", dy / 1000)
        ps.set(TILT["x"], tx)
        ps.set(TILT["y"], ty)
        ps.set_gap_um(self.gap0 if gap is None else gap)
        return ps.popd()

    def loss(self, **kw):
        return float(AN.loss_db(self.popd(**kw)[0]))

    def chief_ray(self):
        """Chief ray of field 1 at the lens (surface 2) and the fibre (3): position and exit angle in air."""
        mfe, T = self.ps.sys.MFE, self.zos.Editors.MFE.MeritOperandType
        g = lambda op, s: mfe.GetOperandValue(op, s, 1, 1.0, 0.0, 0.0, 0.0, 0, 0)
        x2, z2, x3 = g(T.REAX, 2), g(T.REAZ, 2), g(T.REAX, 3)
        return dict(x_lens_um=x2 * 1000, x_fibre_um=x3 * 1000,
                    exit_angle_air_deg=math.degrees(math.atan2(x3 - x2, self.gap0 / 1000 - z2)),
                    stop_surface=self.ps.sys.LDE.StopSurface)

    def window(self):
        """POP report at the fibre and the lens: grid/window as configured, power at the window edge."""
        out = {}
        for end in (1, 2, 3):
            self.ps.set("POP_END", end)
            rep = self.ps.report(strict=False)      # width checked below: oblique refraction at the lens is real
            dg = self.ps.pop.GetResults().GetDataGrid(0)
            v = np.fromiter(dg.Values, float, count=dg.Nx * dg.Ny).reshape(dg.Ny, dg.Nx)
            b = max(4, dg.Nx // 64)
            edge = (v.sum() - v[b:-b, b:-b].sum()) / v.sum()
            yy, xx = np.indices(v.shape)
            cx = ((xx - dg.Nx / 2) * v).sum() / v.sum() * dg.Dx * 1000
            cy = ((yy - dg.Ny / 2) * v).sum() / v.sum() * dg.Dy * 1000
            out["surface_%d" % end] = dict(grid=rep["grid"], width_mm=rep["width_mm"], pilot_size_um=rep["pilot_size_um"],
                                           beam_width_um=rep["beam_width_um"], edge_power_fraction=float(edge),
                                           centroid_um=[float(cx), float(cy)], eff=rep["eff"])
        self.ps.set("POP_END", self.cfgd["pop_tokens"]["POP_END"])
        w0 = self.ps.tokens["POP_WIDEX"]
        if abs(out["surface_1"]["width_mm"] / w0 - 1) > 1e-4:
            raise C.SettingsMismatch("launch window %.6f mm != %.6f: launch stretch not compensated" % (
                out["surface_1"]["width_mm"], w0))
        for s in ("surface_2", "surface_3"):
            if abs(out[s]["width_mm"] / w0 - 1) > 0.01 or out[s]["edge_power_fraction"] > 1e-6:
                raise C.SettingsMismatch("POP window does not hold the beam at %s: %s" % (s, out[s]))
        return out

    def sweep(self, key, values, **fixed):
        rows = []
        for v in values:
            e = self.popd(**{**fixed, key: v})
            rows.append((v, *e, float(AN.loss_db(e[0]))))
        return np.array(rows)

    def convergence(self, tag, points):
        """Project test with the source swapped together with the grid: grid x2 (half pixel) and
        grid x2 + window x2 (same pixel). The .zbf variants are written by fdtd_source.py."""
        s0, w0 = self.ps.tokens["POP_SAMPX"], self.ps.tokens["POP_WIDEX"]
        res = []
        for label, kw in points:
            etas = {}
            for name, s_, w_, suf in (("base", s0, w0, ""), ("grid x2", s0 + 1, w0, "_g2"),
                                      ("grid x2, window x2", s0 + 1, 2 * w0, "_g2w2")):
                for k in ("POP_SAMPX", "POP_SAMPY"):
                    self.ps.set(k, s_)
                for k in ("POP_WIDEX", "POP_WIDEY"):
                    self.ps.set(k, w_)
                self.set_source("P02_%s%s.ZBF" % (tag, suf))
                etas[name] = self.popd(**kw)[0]
                if name != "base":
                    rep = self.ps.report(strict=False)
                    etas[name + " (report grid, width)"] = [rep["grid"], rep["width_mm"]]
            ch = {k: 100 * (v - etas["base"]) / etas["base"] for k, v in etas.items()
                  if k in ("grid x2", "grid x2, window x2")}
            res.append(dict(point=label, eta=etas, change_pct=ch, passed=max(abs(v) for v in ch.values()) < 0.5))
            for k in ("POP_SAMPX", "POP_SAMPY"):
                self.ps.set(k, s0)
            for k in ("POP_WIDEX", "POP_WIDEY"):
                self.ps.set(k, w0)
            self.set_source("P02_%s.ZBF" % tag)
        return res

    def close(self):
        self.popd()
        self.ps.close()


def write_sweep(path, arr, xname):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([xname, "eta", "S", "T", "loss_dB"])
        w.writerows([["%.6g" % r[0]] + ["%.9f" % v for v in r[1:]] for r in arr])


# ======================================================================== controls
def control_elliptical(app, zos):
    """Gaussian .zbf with waist 3.5 um (x) x 4.6 um (y) in B: the pipeline must show the predicted asymmetry,
    and the tilt-token assignment (TILT) must put the narrower angular tolerance on x."""
    c = source_grid()
    FS.write_zbf(FS.BEAMFILES / "P02_ctrl_ellipse.ZBF", gaussian_field(c, 3.5, 4.6), FS.P["zbf_dx_um"], AN.LAM,
                 3.5, 4.6, AN.N_SI)
    ck = source_grid(half=30.0)
    pred = predict(gaussian_field(ck, 3.5, 4.6), ck)
    log("control ellipse: PREDICTED loss %.4f dB; lateral x %.3f y %.3f um; angular x %.4f y %.4f deg (rel)" % (
        pred["loss_dB"], pred["lateral_x_rel"], pred["lateral_y_rel"], pred["angular_x_rel"], pred["angular_y_rel"]))
    cfgd = json.loads(B_CONFIG.read_text())
    cfgd["pop_tokens"]["POP_BEAMTYPE"] = 4
    cfgd["pop_strings"] = {"POP_SOURCEFILE": "P02_ctrl_ellipse.ZBF"}
    b = Bench8b(app, zos, cfgd)
    try:
        l0 = b.loss()
        m = {}
        for name, key, g in (("lateral_x", "dx", pred["lateral_x_rel"]), ("lateral_y", "dy", pred["lateral_y_rel"]),
                             ("angular_x", "tx", pred["angular_x_rel"]), ("angular_y", "ty", pred["angular_y_rel"])):
            m[name] = AN.level_crossing(lambda v: b.loss(**{key: v}), l0 + 1.0, 0.8 * g, 80.0, l0)
    finally:
        b.close()
    m["loss_dB"] = l0
    dev = {k: 100 * (m[k] / pred[k + ("_rel" if k != "loss_dB" else "")] - 1) for k in m if k != "loss_dB"}
    log("control ellipse: MEASURED loss %.4f dB; lateral x %.3f y %.3f um; angular x %.4f y %.4f deg; deviation %%: %s" % (
        l0, m["lateral_x"], m["lateral_y"], m["angular_x"], m["angular_y"], {k: round(v, 2) for k, v in dev.items()}))
    if not (m["angular_x"] < m["angular_y"] and m["lateral_x"] > m["lateral_y"]):
        raise RuntimeError("control ellipse: asymmetry not where predicted (tilt token assignment wrong?) %s" % m)
    if max(abs(v) for v in dev.values()) > 2.0:
        raise RuntimeError("control ellipse: Zemax and prediction disagree by more than 2%%: %s" % dev)
    return dict(predicted=pred, measured=m, deviation_pct=dev)


def control_convention(app, zos):
    """Does OpticStudio read the .zbf with the same axes and the same phase sign as Lumerical writes it?
    Two Gaussian sources: a +1 deg linear phase (Lumerical e^{+ikz}: tilted towards +x) and a +3 um amplitude
    decentre. The predicted best fibre X decenters (reciprocity model) must be reproduced in sign and size."""
    from scipy.optimize import minimize_scalar
    c = source_grid()
    X, Y = np.meshgrid(c, c)
    k = 2 * math.pi * AN.N_SI / AN.LAM
    tests = {"tilt_plus_1deg": np.exp(-(X ** 2 + Y ** 2) / AN.W0 ** 2 + 1j * k * math.sin(math.radians(1.0)) * X),
             "decentre_plus_3um": np.exp(-((X - 3.0) ** 2 + Y ** 2) / AN.W0 ** 2).astype(complex)}
    ck = source_grid(half=30.0)
    Xk, Yk = np.meshgrid(ck, ck)
    small = {"tilt_plus_1deg": np.exp(-(Xk ** 2 + Yk ** 2) / AN.W0 ** 2 + 1j * k * math.sin(math.radians(1.0)) * Xk),
             "decentre_plus_3um": np.exp(-((Xk - 3.0) ** 2 + Yk ** 2) / AN.W0 ** 2).astype(complex)}
    out = {}
    cfgd = json.loads(B_CONFIG.read_text())
    cfgd["pop_tokens"]["POP_BEAMTYPE"] = 4
    cfgd["pop_strings"] = {"POP_SOURCEFILE": "P02_ctrl_tilt_plus_1deg.ZBF"}
    for name, F in tests.items():
        FS.write_zbf(FS.BEAMFILES / ("P02_ctrl_%s.ZBF" % name), F, FS.P["zbf_dx_um"], AN.LAM, AN.W0, AN.W0, AN.N_SI)
    b = Bench8b(app, zos, cfgd)
    try:
        for name in tests:
            pr = minimize_scalar(lambda d: predicted_loss(small[name], ck, dx=d), bounds=(-30, 30), method="bounded",
                                 options=dict(xatol=1e-4)).x
            b.set_source("P02_ctrl_%s.ZBF" % name)
            zm = minimize_scalar(lambda d: b.loss(dx=d), bounds=(-30, 30), method="bounded", options=dict(xatol=1e-3)).x
            out[name] = dict(predicted_dx_um=float(pr), zemax_dx_um=float(zm))
            log("control convention %s: best fibre X predicted %.3f um, OpticStudio %.3f um" % (name, pr, zm))
    finally:
        b.close()
    for name, r in out.items():
        if r["predicted_dx_um"] * r["zemax_dx_um"] < 0 or abs(r["zemax_dx_um"] / r["predicted_dx_um"] - 1) > 0.02:
            raise RuntimeError("control convention %s failed: %s (axis or phase-sign mismatch)" % (name, r))
    return out


def control_field_angle(app, zos, angle_si_deg):
    """B's 4.6 um Gaussian as a .zbf with the launch-stretch compensation, at the FDTD field angle. It must come
    back round at the lens, with the launch window exactly as configured, and with B's tolerances apart from the
    (real) oblique refraction at the lens."""
    angle_air = math.degrees(math.asin(AN.N_SI * math.sin(math.radians(angle_si_deg))))
    c = source_grid(n=FS.P["zbf_n"])
    FS.write_zbf(FS.BEAMFILES / "P02_ctrl_field_angle.ZBF", gaussian_field(c, AN.W0, AN.W0), FS.P["zbf_dx_um"], AN.LAM,
                 AN.W0, AN.W0, AN.N_SI, x_factor=FS.launch_x_factor(angle_si_deg))
    cfgd, _ = make_config(app, zos, "ctrl_field_angle", angle_air, dict(control="B Gaussian, compensated, at field angle"))
    b = Bench8b(app, zos, cfgd)
    try:
        win = b.window()
        l0 = b.loss()
        m = dict(loss_dB=l0)
        for name, key, g in (("lateral_x", "dx", 8.14), ("lateral_y", "dy", 8.14), ("angular_x", "tx", 0.676),
                             ("angular_y", "ty", 0.676)):
            m[name] = float(np.mean([AN.level_crossing(lambda v: b.loss(**{key: s * v}), l0 + 1.0, 0.8 * g, 80.0, l0)
                                     for s in (1, -1)]))
    finally:
        b.close()
    bw = win["surface_2"]["beam_width_um"]
    out = dict(angle_si_deg=angle_si_deg, angle_air_deg=angle_air, window=win, measured=m,
               lens_beam_width_ratio_x_over_y=bw[0] / bw[1],
               vs_B_pct={k: 100 * (m[k] / B_REF[k.split("_")[0] + "_rel"] - 1) for k in m if k != "loss_dB"})
    log("control field angle (%.3f deg in Si, %.3f in air): launch window %.6f mm; lens beam X/Y %.3f/%.3f um "
        "(ratio %.4f); loss %.5f dB; lateral x %.4f y %.4f um, angular x %.4f y %.4f deg; vs B %%: %s" % (
            angle_si_deg, angle_air, win["surface_1"]["width_mm"], bw[0], bw[1], out["lens_beam_width_ratio_x_over_y"],
            l0, m["lateral_x"], m["lateral_y"], m["angular_x"], m["angular_y"],
            {k: round(v, 3) for k, v in out["vs_B_pct"].items()}))
    if abs(out["lens_beam_width_ratio_x_over_y"] - 1) > 0.005:
        raise RuntimeError("control field angle: beam not round at the lens after compensation: %s" % bw)
    return out


# ============================================================================= run
def fmt_pred(p):
    return ("loss %.4f dB; lateral x %.3f (+%.3f/-%.3f) y %.3f um; angular x %.4f (+%.4f/-%.4f) y %.4f deg; "
            "longitudinal +%.1f um (rel); abs lateral x %.3f y %.3f" % (
                p["loss_dB"], p["lateral_x_rel"], *p["lateral_x_rel_sides"], p["lateral_y_rel"], p["angular_x_rel"],
                *p["angular_x_rel_sides"], p["angular_y_rel"], p["longitudinal_rel"], p["lateral_x_abs"], p["lateral_y_abs"]))


def run_config_beam(app, zos, tag, full=True):
    """B_fdtd for one FDTD beam: chief-ray nominal, re-pointed optimum, sweeps about the optimum, convergence."""
    from scipy.optimize import minimize
    src = json.loads((OUT / (tag + ".json")).read_text())
    angle_air = src["angle_air_equivalent_deg"]
    src_summary = {k: src[k] for k in ("tag", "angle_si_deg", "angle_si_peak_deg", "angle_air_equivalent_deg",
                                       "rotated_radii_um", "ellipticity_wx_over_wy", "residual_tilt_deg",
                                       "n_si_fdtd", "zbf_variants")}
    E, c = fdtd_field(tag)
    pred_cr = predict(E, c)
    opt_p = predicted_optimum(E, c)
    pred = predict(E, c, base=opt_p)
    log("%s: PREDICTED on chief ray: %s" % (tag, fmt_pred(pred_cr)))
    log("%s: PREDICTED optimum fibre dx %.3f um, tilt %.4f deg: %s" % (tag, opt_p["dx"], opt_p["tx"], fmt_pred(pred)))
    cfgd, cfg_path = make_config(app, zos, tag, angle_air, src_summary)
    b = Bench8b(app, zos, cfgd)
    t0 = time.time()
    res = dict(tag=tag, predicted_chief_ray=pred_cr, predicted_optimum=opt_p, predicted=pred,
               config=str(cfg_path.relative_to(ROOT)))
    try:
        e0 = b.popd()
        res.update(chief_ray_nominal=dict(eta=e0[0], S=e0[1], T=e0[2], loss_dB=float(AN.loss_db(e0[0]))),
                   chief_ray=b.chief_ray(), window=b.window())
        log("%s: on chief ray eta %.6f (S %.6f T %.6f) = %.4f dB; chief ray at lens %.2f um, exit %.3f deg in air; "
            "stop surface %d; window edge power fibre %.1e lens %.1e" % (
                tag, e0[0], e0[1], e0[2], res["chief_ray_nominal"]["loss_dB"], res["chief_ray"]["x_lens_um"],
                res["chief_ray"]["exit_angle_air_deg"], res["chief_ray"]["stop_surface"],
                res["window"]["surface_3"]["edge_power_fraction"], res["window"]["surface_2"]["edge_power_fraction"]))
        r = minimize(lambda v: b.loss(dx=v[0], tx=v[1]), [opt_p["dx"], opt_p["tx"]], method="Nelder-Mead",
                     options=dict(xatol=2e-3, fatol=1e-6))
        ox, ot = float(r.x[0]), float(r.x[1])
        eo = b.popd(dx=ox, tx=ot)
        res["optimum"] = dict(dx=ox, tx=ot, eta=eo[0], S=eo[1], T=eo[2], loss_dB=float(AN.loss_db(eo[0])), nfev=r.nfev)
        log("%s: MEASURED optimum fibre dx %.3f um, tilt %.4f deg: eta %.6f = %.4f dB (%d evaluations)" % (
            tag, ox, ot, eo[0], res["optimum"]["loss_dB"], r.nfev))
        lat = list(np.round(np.arange(-24.0, 24.0001, 0.5 if full else 1.0), 6))
        ang = list(np.round(np.arange(-1.6, 1.6001, 0.04 if full else 0.08), 6))
        gap = list(np.round(np.arange(0.0, 2000.0001, 20.0 if full else 50.0), 6))
        sw = {"lateral_x": b.sweep("dx", [ox + v for v in lat], tx=ot),
              "lateral_y": b.sweep("dy", lat, dx=ox, tx=ot),
              "angular_x": b.sweep("tx", [ot + v for v in ang], dx=ox),
              "angular_y": b.sweep("ty", ang, dx=ox, tx=ot),
              "longitudinal": b.sweep("gap", gap, dx=ox, tx=ot)}
        centre = {"lateral_x": ox, "lateral_y": 0.0, "angular_x": ot, "angular_y": 0.0, "longitudinal": 20.0}
        tol = {}
        for k, arr in sw.items():
            write_sweep(OUT / ("%s_%s.csv" % (tag, k)), arr, {"longitudinal": "gap_um"}.get(k, "position"))
            tol[k] = AN.tolerance(arr[:, 0], arr[:, 4], centre[k])
        lp, ll = AN.peak(sw["longitudinal"][:, 0], sw["longitudinal"][:, 4])
        res.update(sweeps=tol, longitudinal_peak=dict(gap_um=lp, loss_dB=ll))
        for conv in ("rel", "abs"):
            for ax in "xy":
                res["product_%s_%s" % (ax, conv)] = tol["lateral_" + ax][conv] * math.radians(tol["angular_" + ax][conv])
        log("%s: MEASURED about the optimum: lateral x %.3f (+%.3f/-%.3f) y %.3f (+%.3f/-%.3f) um; angular x %.4f "
            "(+%.4f/-%.4f) y %.4f deg; longitudinal +%.1f um (rel; peak at gap %.1f um); abs: lateral x %.3f y %.3f, "
            "angular x %.4f y %.4f; product x %.5f y %.5f um rad" % (
                tag, tol["lateral_x"]["rel"], tol["lateral_x"]["rel_plus"], tol["lateral_x"]["rel_minus"],
                tol["lateral_y"]["rel"], tol["lateral_y"]["rel_plus"], tol["lateral_y"]["rel_minus"],
                tol["angular_x"]["rel"], tol["angular_x"]["rel_plus"], tol["angular_x"]["rel_minus"],
                tol["angular_y"]["rel"], tol["longitudinal"]["rel_plus"], lp, tol["lateral_x"]["abs"],
                tol["lateral_y"]["abs"], tol["angular_x"]["abs"], tol["angular_y"]["abs"],
                res["product_x_rel"], res["product_y_rel"]))
        conv_pts = [("optimum", dict(dx=ox, tx=ot)), ("lateral x end (+24 um)", dict(dx=ox + 24.0, tx=ot)),
                    ("lateral y end (+24 um)", dict(dx=ox, dy=24.0, tx=ot)),
                    ("angular x end (+1.6 deg)", dict(dx=ox, tx=ot + 1.6)), ("gap end (2000 um)", dict(dx=ox, tx=ot, gap=2000.0))]
        res["convergence"] = b.convergence(tag, conv_pts if full else conv_pts[:2])
        log("%s: convergence (source re-sampled with the grid): %s" % (
            tag, {c_["point"]: {k: round(v, 4) for k, v in c_["change_pct"].items()} for c_ in res["convergence"]}))
    finally:
        b.close()
    res["wall_time_s"] = time.time() - t0
    (OUT / ("%s_zemax.json" % tag)).write_text(json.dumps(res, indent=2, default=float))
    return res, cfg_path


def run_all(app, zos, tag, controls=True, full=True, tol_map=True):
    import tolerance_map as TM
    st = predictor_selftest()
    log("predictor self-test (B Gaussian via reciprocity): lateral %.4f um, angular %.4f deg, longitudinal %.1f um: OK"
        % (st["lateral_x_rel"], st["angular_x_rel"], st["longitudinal_rel"]))
    out = dict(selftest=st)
    if controls:
        out["control_convention"] = control_convention(app, zos)
        out["control_ellipse"] = control_elliptical(app, zos)
        out["control_field_angle"] = control_field_angle(app, zos, json.loads((OUT / (tag + ".json")).read_text())["angle_si_deg"])
    out["beam"], cfg_path = run_config_beam(app, zos, tag, full=full)
    if tol_map:
        o = out["beam"]["optimum"]
        out["map"] = TM.run_map(app, zos, cfg_path, tag="B_fdtd_" + tag, ref_radius_um=B_REF["lateral_rel"],
                                half=12, step=1.0, n_angles=16, r_max=12, r_step=0.5,
                                centre_um=(o["dx"], 0.0), fixed_tokens={TILT["x"]: o["tx"]})
        m = out["map"]
        log("%s: 2-D map about the optimum: 1-dB radius %.3f-%.3f um (x %.3f, y %.3f), circularity %.4f, "
            "A_1dB rel %.1f um^2 (%.2f%% vs B's pi r^2)" % (
                tag, m["r_rel_min_um"], m["r_rel_max_um"], m["x_rel_um"], m["y_rel_um"], m["circularity"],
                m["A1dB_rel_from_cuts_um2"], m["area_rel_vs_pi_r2_pct"]))
    (OUT / ("%s_stage8b.json" % tag)).write_text(json.dumps(out, indent=2, default=float))
    return out
