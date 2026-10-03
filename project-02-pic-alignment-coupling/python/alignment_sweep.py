"""Lateral (X) alignment sweep of configuration A0, with 1-dB tolerance extraction.

Sweeps the receiver-fiber X decenter in OpticStudio POP, logs POPD 0 (total eta),
1 (system S) and 2 (receiver T) at every point, and extracts the 1-dB lateral
tolerance by spline interpolation, not by nearest sample. Repeats for several
air gaps and writes CSVs, a tolerance-versus-gap plot and a run config.

Everything about the POP setup comes from zemax/baseline/run_config.json
("pop_tokens"). The POP .CFG is generated from that file in code, and every
setting is read back out of the CFG and checked; any mismatch raises.

Usage:
    python alignment_sweep.py --mode standalone      # long sweeps (headless)
    python alignment_sweep.py --mode extension       # needs Interactive Extension armed
or import and call run(app) with an already-connected application.
"""
import argparse, csv, glob, json, math, os, re, struct, sys, tempfile, time
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "zemax" / "baseline" / "run_config.json"
MODEL = ROOT / "zemax" / "baseline" / "A0_conventional.zmx"
OUT = ROOT / "results" / "coupling_curves"

INT_TOKENS = {"POP_WAVE", "POP_FIELD", "POP_START", "POP_END", "POP_BEAMTYPE",
              "POP_SAMPX", "POP_SAMPY", "POP_COMPUTE", "POP_FIBERTYPE"}
ONE_DB = math.log(10) / 10          # ln of 10^0.1: eta falls by 1 dB when its log drops by this


class SettingsMismatch(AssertionError):
    """A POP setting did not read back as written."""


# --------------------------------------------------------------------- analytics
def analytic(gap_um, w0=4.6, w2=4.6, lam=1.31):
    """Gaussian-to-Gaussian overlap for a beam of waist w0 after gap_um of air.

    Returns on-axis eta0, the exact lateral 1-dB tolerance (including wavefront
    curvature) and the simple 0.339*sqrt(w1^2+w2^2) estimate. Lengths in um.
    """
    zr = math.pi * w0 ** 2 / lam
    w1 = w0 * math.sqrt(1 + (gap_um / zr) ** 2)
    a1 = 1 / w1 ** 2 + (1j * (2 * math.pi / lam) / (2 * gap_um * (1 + (zr / gap_um) ** 2)) if gap_um else 0)
    a2 = 1 / w2 ** 2
    s = a1 + a2
    eta0 = 4 / (w1 ** 2 * w2 ** 2 * abs(s) ** 2)
    c = 2 * (a1 * a2 / s).real                       # eta(d) = eta0 * exp(-c d^2)
    return dict(w1=w1, eta0=eta0, c=c,
                d_exact=math.sqrt(ONE_DB / c),
                d_abs=math.sqrt(math.log(eta0 / 10 ** -0.1) / c) if eta0 > 10 ** -0.1 else float("nan"),
                d_simple=0.339 * math.sqrt(w1 ** 2 + w2 ** 2))


# --------------------------------------------------------------- tolerance extraction
def crossing(d, loss, level):
    """Interpolated d (>0) where loss first reaches `level`, by cubic spline + root find."""
    d, loss = np.asarray(d), np.asarray(loss)
    over = np.nonzero(loss >= level)[0]
    if len(over) == 0 or over[0] == 0:
        return float("nan")
    i = over[0]
    spl = CubicSpline(d, loss)
    return brentq(lambda x: float(spl(x)) - level, d[i - 1], d[i], xtol=1e-9)


def extract_tolerance(d_um, loss_db):
    """1-dB lateral tolerance from a symmetric sweep. Returns relative-to-peak and absolute.

    relative: loss rises 1 dB above the on-axis value L(0) (primary convention).
    absolute: total loss reaches 1.000 dB.
    Both sides are extracted separately; the mean is reported and the spread kept.
    """
    d_um, loss_db = np.asarray(d_um), np.asarray(loss_db)
    i0 = int(np.argmin(np.abs(d_um)))
    if abs(d_um[i0]) > 1e-9:
        raise ValueError("sweep has no d = 0 point")
    l0 = loss_db[i0]
    if int(np.argmin(loss_db)) != i0:
        raise RuntimeError("loss minimum is not at d = 0")
    out = {}
    for name, level in (("rel", l0 + 1.0), ("abs", 1.0)):
        sides = []
        for sign in (+1, -1):
            m = sign * d_um >= 0
            order = np.argsort(sign * d_um[m])
            sides.append(crossing(sign * d_um[m][order], loss_db[m][order], level))
        out[name] = float(np.nanmean(sides)) if not all(math.isnan(s) for s in sides) else float("nan")
        out[name + "_sides"] = sides
    out["L0"] = float(l0)
    return out


def selftest():
    """The extractor must recover the analytic tolerance from analytic data."""
    a = analytic(20.0)
    d = np.round(np.arange(-6, 6.0001, 0.1), 6)
    eta = a["eta0"] * np.exp(-a["c"] * d ** 2)
    got = extract_tolerance(d, -10 * np.log10(eta))["rel"]
    if abs(got - a["d_exact"]) > 1e-4:
        raise RuntimeError("tolerance extractor self-test failed: %.6f vs %.6f" % (got, a["d_exact"]))
    return got


# ------------------------------------------------------------------ OpticStudio side
def load_zosapi():
    try:
        import ZOSAPI
        return ZOSAPI
    except ImportError:
        import clr
        roots = sorted(glob.glob(os.path.join(os.environ["ProgramFiles"], "Ansys Zemax OpticStudio*")), reverse=True)
        helper = next(p for r in roots for p in (os.path.join(r, "ZOSAPI_NetHelper.dll"),
                      os.path.join(r, "ZOS-API", "Libraries", "ZOSAPI_NetHelper.dll")) if os.path.isfile(p))
        clr.AddReference(helper)
        import ZOSAPI_NetHelper
        ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize()
        clr.AddReference(os.path.join(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory(), "ZOSAPI.dll"))
        import ZOSAPI
        return ZOSAPI


class PopSession:
    """One POP analysis whose settings file is generated from run_config.json and verified."""

    def __init__(self, app, cfgd, zos):
        self.zos, self.cfgd, self.app = zos, cfgd, app
        self.tokens = {k: v for k, v in cfgd["pop_tokens"].items() if not k.startswith("_")}
        self.sys = app.PrimarySystem
        self.sys.LoadFile(str(MODEL), False)
        if self.sys.LDE.NumberOfSurfaces != 3:
            raise SettingsMismatch("A0 model should have 3 surfaces, found %d" % self.sys.LDE.NumberOfSurfaces)
        for idx, props in cfgd.get("surface_settings", {}).items():
            if idx.startswith("_"):
                continue
            pod = self.sys.LDE.GetSurfaceAt(int(idx)).PhysicalOpticsData
            for name, val in props.items():
                setattr(pod, name, val)
                if getattr(self.sys.LDE.GetSurfaceAt(int(idx)).PhysicalOpticsData, name) != val:
                    raise SettingsMismatch("surface %s %s did not read back as %r" % (idx, name, val))
        self.pop = self.sys.Analyses.New_Analysis(zos.Analysis.AnalysisIDM.PhysicalOpticsPropagation)
        self.st = self.pop.GetSettings()
        tmp = tempfile.gettempdir()
        self.cfg = os.path.join(tmp, "a0_pop_settings.CFG")       # the file POP reads
        self._scratch = os.path.join(tmp, "a0_pop_probe.CFG")      # used only to locate offsets
        self.nchecks = 0
        self.off = {t: self._locate(t) for t in self.tokens}
        self.st.SaveTo(self.cfg)
        for tok, val in self.tokens.items():
            self.st.ModifySettings(self.cfg, tok, repr(val))
        self.verify()
        self._operands()

    # -- settings file generation and read-back ------------------------------------
    def _locate(self, tok):
        """Find where `tok` lives in the CFG by writing two marker values and diffing."""
        for a, b in ((3, 5), (0, 1)) if tok in INT_TOKENS else ((0.1111111, 0.2222222),):
            blobs = []
            for v in (a, b):
                self.st.SaveTo(self._scratch)
                self.st.ModifySettings(self._scratch, tok, repr(v))
                blobs.append(Path(self._scratch).read_bytes())
            diff = [i for i in range(len(blobs[0])) if blobs[0][i] != blobs[1][i]]
            if diff:
                width = 4 if tok in INT_TOKENS else 8
                off = diff[0] // width * width
                fmt = "<i" if tok in INT_TOKENS else "<d"
                if [struct.unpack_from(fmt, x, off)[0] for x in blobs] == [a, b]:
                    return off
        raise SettingsMismatch("could not locate %s inside the CFG file" % tok)

    def _read(self, tok):
        fmt = "<i" if tok in INT_TOKENS else "<d"
        return struct.unpack_from(fmt, Path(self.cfg).read_bytes(), self.off[tok])[0]

    def verify(self):
        """Permanent assertion: every setting reads back from the CFG exactly as written."""
        for tok, val in self.tokens.items():
            got = self._read(tok)
            if got != (int(val) if tok in INT_TOKENS else float(val)):
                raise SettingsMismatch("%s: wrote %r, CFG holds %r" % (tok, val, got))
        self.nchecks += 1

    def set(self, tok, val):
        self.tokens[tok] = val
        self.st.ModifySettings(self.cfg, tok, repr(val))
        self.verify()

    def commit(self):
        """POPD reads the SAVED default settings, so load into the analysis and save."""
        self.st.LoadFrom(self.cfg)
        self.st.Save()

    # -- merit operands ------------------------------------------------------------
    def _operands(self):
        mfe, mc = self.sys.MFE, self.zos.Editors.MFE.MeritColumn
        while mfe.NumberOfOperands < 3:
            mfe.AddOperand()
        for i, data in enumerate((0, 1, 2), start=1):
            op = mfe.GetOperandAt(i)
            op.ChangeType(self.zos.Editors.MFE.MeritOperandType.POPD)
            for col, v in ((mc.Param1, self.tokens["POP_END"]), (mc.Param2, 1), (mc.Param3, 1), (mc.Param4, data)):
                op.GetOperandCell(col).IntegerValue = v
        self.mfe = mfe

    def popd(self):
        self.commit()
        self.mfe.CalculateMeritFunction()
        return tuple(self.mfe.GetOperandAt(i).Value for i in (1, 2, 3))

    # -- model and report checks ---------------------------------------------------
    def set_gap_um(self, gap_um):
        s = self.sys.LDE.GetSurfaceAt(1)
        s.Thickness = gap_um / 1000.0
        if abs(s.Thickness - gap_um / 1000.0) > 1e-12:
            raise SettingsMismatch("gap did not read back: %r" % s.Thickness)

    def report_checks(self, gap_um):
        """Read the POP text report and assert it shows the configuration we asked for."""
        self.commit()
        self.pop.ApplyAndWaitForCompletion()
        f = os.path.join(tempfile.gettempdir(), "a0_pop_report.txt")
        self.pop.GetResults().GetTextFile(f)
        raw = Path(f).read_bytes()
        txt = raw.decode("utf-16" if raw[:2] in (b"\xff\xfe", b"\xfe\xff") else "latin-1")
        t, a = self.tokens, analytic(gap_um, w0=self.tokens["POP_PARAM1"] * 1000, w2=self.tokens["POP_FPARAM1"] * 1000,
                                     lam=self.cfgd["wavelength_um"])
        g = lambda pat: re.search(pat, txt)

        def need(cond, msg):
            if not cond:
                raise SettingsMismatch("report: " + msg)

        n = int(g(r"Grid size \(X by Y\): (\d+) by (\d+)").group(1))
        need(n == 2 ** (t["POP_SAMPX"] + 4), "grid %d != 2^(SAMPX+4)" % n)
        need(abs(float(g(r"Display X Width = ([\d.Ee+-]+)").group(1)) - t["POP_WIDEX"]) < 1e-9, "window width")
        need(abs(float(g(r"Beam wavelength is ([\d.]+)").group(1)) - self.cfgd["wavelength_um"]) < 1e-6, "wavelength")
        need(int(g(r"Total Irradiance surface (\d+)").group(1)) == t["POP_END"], "end surface")
        size, waist, pos = (float(x) for x in g(r"Pilot: Size= ([\d.Ee+-]+), Waist= ([\d.Ee+-]+), Pos= ([\d.Ee+-]+)").groups())
        need(abs(waist - t["POP_PARAM1"]) / t["POP_PARAM1"] < 1e-4, "source waist")
        need(abs(pos - gap_um / 1000) < 1e-7, "pilot position != gap")
        need(abs(size * 1000 - a["w1"]) / a["w1"] < 1e-4, "beam radius at fiber != analytic")
        need(t["POP_WIDEX"] * 1000 >= 4 * a["w1"], "window narrower than 4x beam radius")
        return dict(grid=n, pilot_size_um=size * 1000, pilot_waist_um=waist * 1000, pilot_pos_um=pos * 1000)

    def close(self):
        self.pop.Close()


# ------------------------------------------------------------------------- sweeps
def sweep_gap(ps, gap_um, step_um, range_um, log):
    ps.set_gap_um(gap_um)
    rep = ps.report_checks(gap_um)
    ps.set("POP_FPARAM3", 0.0)
    a = analytic(gap_um)
    eta0 = ps.popd()
    if abs(eta0[0] - a["eta0"]) > 1e-4:
        raise SettingsMismatch("gap %g: on-axis eta %.6f vs analytic %.6f (receiver waist?)" % (gap_um, eta0[0], a["eta0"]))
    pts = []
    for d in np.round(np.arange(-range_um, range_um + step_um / 2, step_um), 6):
        ps.set("POP_FPARAM3", float(d) / 1000.0)
        eta, s, t = ps.popd()
        pts.append((float(d), eta, s, t, -10 * math.log10(eta)))
    ps.set("POP_FPARAM3", 0.0)
    arr = np.array(pts)
    # symmetry: +d and -d must agree
    asym = float(np.max(np.abs(arr[:, 1] - arr[::-1, 1])))
    if asym > 1e-6:
        raise RuntimeError("gap %g: eta(+d) != eta(-d), max diff %.2e" % (gap_um, asym))
    tol = extract_tolerance(arr[:, 0], arr[:, 4])
    # convergence: double the grid at d = 0 and d = predicted tolerance
    conv = {}
    for d in (0.0, round(a["d_exact"], 3)):
        ps.set("POP_FPARAM3", d / 1000.0)
        base = ps.popd()[0]
        for k in ("POP_SAMPX", "POP_SAMPY"):
            ps.set(k, ps.cfgd["sweep_lateral"]["convergence_samp_index"])
        fine = ps.popd()[0]
        for k in ("POP_SAMPX", "POP_SAMPY"):
            ps.set(k, ps.cfgd["pop_tokens"][k])
        conv[d] = 100 * (fine - base) / base
        if abs(conv[d]) > 0.5:
            raise RuntimeError("gap %g: doubling the grid changed eta by %.3f%% at d=%g" % (gap_um, conv[d], d))
    ps.set("POP_FPARAM3", 0.0)
    return arr, tol, conv, rep, asym, a


def run(app, out=OUT, zos=None):
    zos = zos or load_zosapi()
    log = lambda m: print(m, flush=True)
    cfgd = json.loads(CONFIG.read_text())
    sw = cfgd["sweep_lateral"]
    log("extractor self-test (analytic data, gap 20 um): %.6f um" % selftest())
    out.mkdir(parents=True, exist_ok=True)
    ps = PopSession(app, cfgd, zos)
    nominal_loss = None
    rows, t0 = [], time.time()
    try:
        # reproducibility: the generated CFG must reproduce the committed Stage 3 value
        ps.set_gap_um(sw["nominal_gap_um"])
        ref = float(list(csv.DictReader(open(ROOT / "zemax" / "baseline" / "a0_baseline_popd.csv")))[0]["popd0_eta_total"])
        got = ps.popd()[0]
        if abs(got - ref) > 1e-6:
            raise SettingsMismatch("generated CFG gives eta %.9f, Stage 3 committed %.9f" % (got, ref))
        log("generated CFG reproduces Stage 3 eta = %.9f" % got)
        for gap in sw["gaps_um"]:
            step = sw["nominal_step_um"] if gap == sw["nominal_gap_um"] else sw["step_um"]
            arr, tol, conv, rep, asym, a = sweep_gap(ps, gap, step, sw["range_um"], log)
            with open(out / ("a0_lateral_gap%03dum.csv" % gap), "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["x_decenter_um", "popd0_eta_total", "popd1_S_system", "popd2_T_receiver", "loss_dB"])
                for r in arr:
                    w.writerow(["%.4f" % r[0]] + ["%.9f" % v for v in r[1:4]] + ["%.6f" % r[4]])
            rows.append(dict(gap_um=gap, step_um=step, eta0=arr[len(arr) // 2, 1], L0_dB=tol["L0"],
                             d1dB_rel_um=tol["rel"], d1dB_rel_minus_um=tol["rel_sides"][1], d1dB_rel_plus_um=tol["rel_sides"][0],
                             d1dB_abs_um=tol["abs"], pred_exact_um=a["d_exact"], pred_simple_um=a["d_simple"],
                             pred_abs_um=a["d_abs"], beam_radius_um=a["w1"],
                             diff_vs_exact_pct=100 * (tol["rel"] - a["d_exact"]) / a["d_exact"],
                             diff_vs_simple_pct=100 * (tol["rel"] - a["d_simple"]) / a["d_simple"],
                             conv_pct_d0=conv[0.0], conv_pct_dtol=[v for k, v in conv.items() if k != 0.0][0],
                             asym=asym))
            log("gap %3g um: d1dB = %.4f um (exact %.4f, simple %.4f)  S min %.6f  conv %.4f%%/%.4f%%  [%.0fs]" % (
                gap, tol["rel"], a["d_exact"], a["d_simple"], arr[:, 2].min(), rows[-1]["conv_pct_d0"],
                rows[-1]["conv_pct_dtol"], time.time() - t0))
        # interpolation resolution check at the nominal gap: thin the 0.05 um data to 0.2 um
        nom = np.loadtxt(out / ("a0_lateral_gap%03dum.csv" % sw["nominal_gap_um"]), delimiter=",", skiprows=1)
        thin = nom[::4] if abs(nom[::4][len(nom[::4]) // 2, 0]) < 1e-9 else nom[::2]
        res_check = extract_tolerance(thin[:, 0], thin[:, 4])["rel"] - extract_tolerance(nom[:, 0], nom[:, 4])["rel"]
        log("resolution check (0.05 um vs thinned step): change %.2e um" % res_check)
    finally:
        ps.set_gap_um(sw["nominal_gap_um"])
        ps.set("POP_FPARAM3", 0.0)
        ps.commit()
        ps.close()
    keys = list(rows[0].keys())
    with open(out / "a0_tolerance_vs_gap.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (("%.6f" % v) if isinstance(v, float) else v) for k, v in r.items()})
    plot_tolerance(rows, out / "a0_tolerance_vs_gap.png")
    meta = dict(stage="5", description="A0 lateral X decenter of the receiver", date=time.strftime("%Y-%m-%d"),
                opticstudio=dict(build=str(app.OpticStudioVersion), license=str(app.LicenseStatus), mode=str(app.Mode)),
                source_config=str(CONFIG.relative_to(ROOT)), model=str(MODEL.relative_to(ROOT)),
                pop_tokens=dict(ps.tokens),
                sweep=sw, convergence_rule="double grid, <0.5% change in eta, checked at d=0 and d=tolerance for every gap",
                settings_readbacks_passed=ps.nchecks, resample_after_refraction="False on all surfaces (no refracting surface in A0)", surface_settings=cfgd.get("surface_settings"),
                tolerance_method="cubic spline of loss(d), brentq root; both sides, mean reported",
                tolerance_conventions={"rel": "loss rises 1 dB above on-axis value (primary)", "abs": "total loss reaches 1.000 dB"},
                interpolation_resolution_check_um=res_check, wall_time_s=time.time() - t0)
    (out / "run_config_a0_lateral.json").write_text(json.dumps(meta, indent=2, default=str))
    return rows


# --------------------------------------------------------------------------- plot
def plot_tolerance(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e4e3de"
    blue, orange = "#2a78d6", "#eb6834"
    g = np.linspace(0, 105, 300)
    ex = [analytic(x)["d_exact"] for x in g]
    si = [analytic(x)["d_simple"] for x in g]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    ax.plot(g, ex, color=orange, lw=1.8, label="Exact Gaussian overlap (with wavefront curvature)")
    ax.plot(g, si, color=ink2, lw=1.4, ls="--", label="0.339 sqrt(w1^2 + w2^2)")
    ax.plot([r["gap_um"] for r in rows], [r["d1dB_rel_um"] for r in rows], color=blue, lw=0, marker="o", ms=7,
            mec="white", mew=1.5, label="Zemax POP, interpolated", zorder=5)
    for r in rows:
        ax.annotate("%.2f" % r["d1dB_rel_um"], (r["gap_um"], r["d1dB_rel_um"]), textcoords="offset points",
                    xytext=(0, -15), ha="center", fontsize=8, color=ink2)
    ax.set_xlabel("Air gap, grating plane to fiber facet (um)", color=ink)
    ax.set_ylabel("Lateral 1-dB tolerance, one side (um)", color=ink)
    ax.set_title("A0 lateral 1-dB tolerance versus air gap", loc="left", fontsize=11, color=ink)
    ax.set_xlim(-3, 108)
    ax.set_ylim(2.0, None)
    ax.grid(color=grid, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(ink2)
    ax.tick_params(colors=ink2)
    ax.legend(frameon=False, loc="upper left", fontsize=8.5, labelcolor=ink)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=("standalone", "extension"), default="standalone")
    args = ap.parse_args()
    zos = load_zosapi()
    conn = zos.ZOSAPI_Connection()
    app = conn.CreateNewApplication() if args.mode == "standalone" else conn.ConnectAsExtension(0)
    if app is None or (args.mode == "extension" and not conn.IsAlive):
        raise SystemExit("could not connect (is the Interactive Extension armed?)")
    try:
        run(app, zos=zos)
    finally:
        if args.mode == "standalone":
            app.CloseApplication()


if __name__ == "__main__":
    main()
