"""Shared POP machinery: connection, verified POP settings, model building,
Gaussian-beam analytics and 1-dB tolerance extraction.

Used by alignment_sweep.py (1-D sweeps) and tolerance_map.py (2-D maps).

Every POP setting comes from a run_config.json ("pop_tokens", "surface_settings").
The binary POP .CFG is generated from it, and every write is read back out of the
CFG; any mismatch raises SettingsMismatch. POPD reads the *saved default* POP
settings, so commit() loads the CFG into the analysis and saves it before any
POPD evaluation.
"""
import glob, math, os, re, struct, tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
INT_TOKENS = {"POP_WAVE", "POP_FIELD", "POP_START", "POP_END", "POP_BEAMTYPE",
              "POP_SAMPX", "POP_SAMPY", "POP_COMPUTE", "POP_FIBERTYPE"}
ONE_DB = math.log(10) / 10          # ln(10^0.1): eta falls 1 dB when ln(eta) drops by this


class SettingsMismatch(AssertionError):
    """A setting did not read back as written, or POP reported something else."""


def loss_db(eta):
    return -10 * np.log10(eta)


# =========================================================== Gaussian-beam analytics
def q_param(w, R, n, lam):
    """Complex beam parameter from 1/e^2 radius w, wavefront radius R, index n."""
    return 1 / ((0 if math.isinf(R) else 1 / R) - 1j * lam / (math.pi * n * w * w))


def w_R(q, n, lam):
    iq = 1 / q
    return math.sqrt(-lam / (math.pi * n * iq.imag)), (math.inf if abs(iq.real) < 1e-18 else 1 / iq.real)


def refract(q, n1, n2, radius):
    """Refraction at a spherical surface (OpticStudio sign convention for radius)."""
    c = 0 if math.isinf(radius) else (n1 - n2) / (radius * n2)
    return q / (c * q + n1 / n2)


def beam_at(surfaces, gap_um, n_after, w0, lam, gap_surface):
    """Pilot-beam (w, R) in um at the last surface for a sequential surface list.

    surfaces: run_config "surfaces" list. The source waist sits on surface 1 and
    propagates through each thickness; refraction at each later surface uses
    n_after[i], the index of the medium after surface i (from OpticStudio).
    """
    q = q_param(w0, math.inf, n_after[1], lam)
    for i in range(1, len(surfaces) - 1):
        if i > 1:
            r = surfaces[i].get("radius_mm", "inf")
            q = refract(q, n_after[i - 1], n_after[i], math.inf if r == "inf" else float(r) * 1000)
        t = gap_um if i == gap_surface else float(surfaces[i]["thickness_mm"]) * 1000
        q = q + t
    return w_R(q, n_after[len(surfaces) - 2], lam)


def overlap(w1, R1, w2, lam, d=0.0, theta_deg=0.0, n=1.0, half=None, npts=40001):
    """Numerical power overlap of a curved Gaussian (w1, R1) with a flat Gaussian
    receiver (w2) decentred by d and tilted by theta. Independent of OpticStudio.
    Separable in x and y; the y factor has neither offset nor tilt."""
    half = half or 8 * max(w1, w2) + abs(d)
    x = np.linspace(-half, half, npts)
    k = 2 * math.pi * n / lam
    curv = 0 if math.isinf(R1) else k / (2 * R1)

    def one(dd, th):
        e1 = np.exp(-x ** 2 / w1 ** 2 - 1j * curv * x ** 2)
        e2 = np.exp(-(x - dd) ** 2 / w2 ** 2 + 1j * k * math.sin(math.radians(th)) * x)
        num = abs(np.trapezoid(e1 * np.conj(e2), x)) ** 2
        return num / (np.trapezoid(abs(e1) ** 2, x) * np.trapezoid(abs(e2) ** 2, x))
    return one(d, theta_deg) * one(0.0, 0.0)


# ================================================================ tolerance extraction
def tolerance(x, loss, x0):
    """1-dB tolerance either side of x0, both conventions, by spline interpolation.

    rel: loss rises 1 dB above its value at x0.   abs: total loss reaches 1.000 dB.
    Returns distances from x0 (positive), nan where the level is not reached inside
    the sweep or, for abs, where L(x0) is already at or above 1 dB.
    """
    x, loss = np.asarray(x, float), np.asarray(loss, float)
    order = np.argsort(x)
    x, loss = x[order], loss[order]
    spl = CubicSpline(x, loss)
    l0 = float(spl(x0))
    out = dict(x0=float(x0), L0=l0)
    for name, level in (("rel", l0 + 1.0), ("abs", 1.0)):
        for side in ("plus", "minus"):
            out[name + "_" + side] = float("nan")
            if level <= l0:
                continue
            xs = x[x > x0 + 1e-12] if side == "plus" else x[x < x0 - 1e-12][::-1]
            hit = np.nonzero(spl(xs) >= level)[0]
            if len(hit):
                j = hit[0]
                prev = x0 if j == 0 else xs[j - 1]
                root = brentq(lambda v: float(spl(v)) - level, min(prev, xs[j]), max(prev, xs[j]), xtol=1e-10)
                out[name + "_" + side] = abs(root - x0)
        vals = [out[name + "_plus"], out[name + "_minus"]]
        out[name] = float(np.nanmean(vals)) if not all(math.isnan(v) for v in vals) else float("nan")
    return out


def peak(x, loss):
    """Interpolated position and value of minimum loss."""
    x, loss = np.asarray(x, float), np.asarray(loss, float)
    order = np.argsort(x)
    x, loss = x[order], loss[order]
    i = int(np.argmin(loss))
    spl = CubicSpline(x, loss)
    lo, hi = x[max(i - 1, 0)], x[min(i + 1, len(x) - 1)]
    if lo == hi:
        return float(x[i]), float(loss[i])
    r = minimize_scalar(lambda v: float(spl(v)), bounds=(lo, hi), method="bounded", options=dict(xatol=1e-9))
    return float(r.x), float(r.fun)


# ======================================================================= OpticStudio
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


def connect(mode, zos):
    conn = zos.ZOSAPI_Connection()
    app = conn.CreateNewApplication() if mode == "standalone" else conn.ConnectAsExtension(0)
    if app is None or (mode == "extension" and not conn.IsAlive):
        raise SystemExit("could not connect (is the Interactive Extension armed?)")
    return conn, app


def build_model(app, cfgd, zos):
    """Create the sequential model described by cfgd["surfaces"] and save it."""
    sysm = app.PrimarySystem
    sysm.New(False)
    sysm.SystemData.Wavelengths.GetWavelength(1).Wavelength = cfgd["wavelength_um"]
    if cfgd.get("glass_catalog_file"):
        # OpticStudio only reads catalogues from its glass directory; install the repo copy there
        src = ROOT / cfgd["glass_catalog_file"]
        dst = Path(app.GlassDir) / src.name
        dst.write_bytes(src.read_bytes())
    for cat in cfgd.get("catalogs", []):
        sysm.SystemData.MaterialCatalogs.AddCatalog(cat)
    lde = sysm.LDE
    while lde.NumberOfSurfaces < len(cfgd["surfaces"]):
        lde.InsertNewSurfaceAt(lde.NumberOfSurfaces - 1)
    for i, spec in enumerate(cfgd["surfaces"]):
        s = lde.GetSurfaceAt(i)
        s.Comment = spec.get("comment", "")[:32]
        if spec.get("type") == "CoordinateBreak":
            s.ChangeType(s.GetSurfaceTypeSettings(zos.Editors.LDE.SurfaceType.CoordinateBreak))
            s.Thickness = float(spec.get("thickness_mm", 0.0))
            continue
        if spec.get("radius_mm", "inf") != "inf":
            s.Radius = float(spec["radius_mm"])
        if "thickness_mm" in spec and spec["thickness_mm"] != "inf":
            s.Thickness = float(spec["thickness_mm"])
        if spec.get("material"):
            s.Material = spec["material"]
    check_model(sysm, cfgd)
    path = ROOT / cfgd["model"]
    path.parent.mkdir(parents=True, exist_ok=True)
    sysm.SaveAs(str(path))
    return path


def check_model(sysm, cfgd):
    """Assert the open model matches cfgd["surfaces"] cell by cell."""
    lde = sysm.LDE
    if lde.NumberOfSurfaces != len(cfgd["surfaces"]):
        raise SettingsMismatch("model has %d surfaces, config %d" % (lde.NumberOfSurfaces, len(cfgd["surfaces"])))
    for i, spec in enumerate(cfgd["surfaces"]):
        s = lde.GetSurfaceAt(i)
        if spec.get("type") == "CoordinateBreak":
            if "Coordinate Break" not in s.TypeName:
                raise SettingsMismatch("surface %d should be a coordinate break, is %s" % (i, s.TypeName))
            continue
        r = spec.get("radius_mm", "inf")
        if (r == "inf") != (abs(s.Radius) > 1e9) or (r != "inf" and abs(s.Radius - float(r)) > 1e-12):
            raise SettingsMismatch("surface %d radius %r != %r" % (i, s.Radius, r))
        if "thickness_mm" in spec and spec["thickness_mm"] != "inf" and abs(s.Thickness - float(spec["thickness_mm"])) > 1e-12:
            raise SettingsMismatch("surface %d thickness %r != %r" % (i, s.Thickness, spec["thickness_mm"]))
        if (spec.get("material") or "").upper() != (s.Material or "").upper():
            raise SettingsMismatch("surface %d material %r != %r" % (i, s.Material, spec.get("material")))


class PopSession:
    """One POP analysis whose settings file is generated from a run_config and verified."""

    def __init__(self, app, cfgd, zos, load=True):
        self.zos, self.cfgd, self.app = zos, cfgd, app
        self.tokens = {k: v for k, v in cfgd["pop_tokens"].items() if not k.startswith("_")}
        self.sys = app.PrimarySystem
        if load:
            self.sys.LoadFile(str(ROOT / cfgd["model"]), False)
        if any("radius_mm" in s for s in cfgd.get("surfaces", [])):   # build-spec configs only (A1, B)
            check_model(self.sys, cfgd)
        self.gap_surface = cfgd.get("gap_surface", 1)
        self.apply_surface_settings(cfgd.get("surface_settings", {}))
        self.pop = self.sys.Analyses.New_Analysis(zos.Analysis.AnalysisIDM.PhysicalOpticsPropagation)
        self.st = self.pop.GetSettings()
        tmp = tempfile.gettempdir()
        self.cfg = os.path.join(tmp, "pop_%s.CFG" % cfgd.get("id", "A0"))
        self._scratch = os.path.join(tmp, "pop_probe.CFG")
        self.nchecks = 0
        self.off = {t: self._locate(t) for t in self.tokens}
        self.st.SaveTo(self.cfg)
        for tok, val in self.tokens.items():
            self.st.ModifySettings(self.cfg, tok, repr(val))
        self.verify()
        self._operands()

    def apply_surface_settings(self, settings):
        for idx, props in settings.items():
            if idx.startswith("_"):
                continue
            for name, val in props.items():
                setattr(self.sys.LDE.GetSurfaceAt(int(idx)).PhysicalOpticsData, name, val)
                if getattr(self.sys.LDE.GetSurfaceAt(int(idx)).PhysicalOpticsData, name) != val:
                    raise SettingsMismatch("surface %s %s did not read back as %r" % (idx, name, val))

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
        return struct.unpack_from("<i" if tok in INT_TOKENS else "<d", Path(self.cfg).read_bytes(), self.off[tok])[0]

    @staticmethod
    def _same(tok, got, val):
        """Integers exactly; doubles to 2 ulp (OpticStudio's .NET parser is not always correctly rounded)."""
        if tok in INT_TOKENS:
            return got == int(val)
        return abs(got - float(val)) <= 2 * math.ulp(float(val))

    def verify(self):
        """Permanent assertion: every setting reads back from the CFG as written (doubles to 2 ulp)."""
        for tok, val in self.tokens.items():
            got = self._read(tok)
            if not self._same(tok, got, val):
                raise SettingsMismatch("%s: wrote %r, CFG holds %r" % (tok, val, got))
        self.nchecks += 1

    def set(self, tok, val):
        val = int(val) if tok in INT_TOKENS else float(val)     # plain Python number: repr(np.float64) is not parseable
        if tok not in self.off:
            self.off[tok] = self._locate(tok)
        text = str(val) if tok in INT_TOKENS else "%.12g" % val
        val = int(text) if tok in INT_TOKENS else float(text)
        self.tokens[tok] = val
        self.st.ModifySettings(self.cfg, tok, text)
        if not self._same(tok, self._read(tok), val):
            # an occasional write does not land (seen once in ~10^4 writes); rewrite once, then verify strictly
            self.retries = getattr(self, "retries", 0) + 1
            self.st.ModifySettings(self.cfg, tok, text)
        self.verify()

    def commit(self):
        """POPD reads the SAVED default settings: load the CFG into the analysis and save."""
        self.st.LoadFrom(self.cfg)
        self.st.Save()

    # -- operands ------------------------------------------------------------------
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
        """(eta_total, S, T) from POPD Data 0, 1, 2."""
        self.commit()
        self.mfe.CalculateMeritFunction()
        return tuple(self.mfe.GetOperandAt(i).Value for i in (1, 2, 3))

    def index_after(self, surface):
        """Refractive index of the medium after `surface` at wavelength 1, from OpticStudio."""
        return self.mfe.GetOperandValue(self.zos.Editors.MFE.MeritOperandType.INDX, surface, 1, 0, 0, 0, 0, 0, 0)

    # -- model ---------------------------------------------------------------------
    def set_thickness_um(self, surface, um):
        s = self.sys.LDE.GetSurfaceAt(surface)
        s.Thickness = um / 1000.0
        if abs(s.Thickness - um / 1000.0) > 1e-12:
            raise SettingsMismatch("surface %d thickness did not read back: %r" % (surface, s.Thickness))

    def set_gap_um(self, um):
        self.set_thickness_um(self.gap_surface, um)

    # -- POP text report -------------------------------------------------------------
    def report(self, strict=True):
        """Run the POP window, parse its text report and assert the generic settings."""
        self.commit()
        self.pop.ApplyAndWaitForCompletion()
        f = os.path.join(tempfile.gettempdir(), "pop_report.txt")
        self.pop.GetResults().GetTextFile(f)
        raw = Path(f).read_bytes()
        txt = raw.decode("utf-16" if raw[:2] in (b"\xff\xfe", b"\xfe\xff") else "latin-1")
        g = lambda pat: re.search(pat, txt)
        t = self.tokens

        def need(cond, msg):
            if strict and not cond:
                raise SettingsMismatch("report: " + msg)
        n = int(g(r"Grid size \(X by Y\): (\d+) by (\d+)").group(1))
        need(n == 2 ** (t["POP_SAMPX"] + 4), "grid %d != 2^(SAMPX+4)" % n)
        width = float(g(r"Display X Width = ([\d.Ee+-]+)").group(1))
        need(abs(width - t["POP_WIDEX"]) / t["POP_WIDEX"] < 1e-3, "window width %.6g mm != configured %.6g mm" % (width, t["POP_WIDEX"]))
        need(abs(float(g(r"Beam wavelength is ([\d.]+)").group(1)) - self.cfgd["wavelength_um"]) < 1e-6, "wavelength")
        need(int(g(r"Total Irradiance surface (\d+)").group(1)) == t["POP_END"], "end surface")
        size, waist, pos = (float(x) for x in g(r"Pilot: Size= ([\d.Ee+-]+), Waist= ([\d.Ee+-]+), Pos= ([\d.Ee+-]+)").groups())
        eff = [float(x) for x in g(r"Fiber Efficiency: System ([\d.Ee+-]+), Receiver ([\d.Ee+-]+), Coupling ([\d.Ee+-]+)").groups()]
        return dict(grid=n, width_mm=width, pilot_size_um=size * 1000, pilot_waist_um=waist * 1000,
                    pilot_pos_um=pos * 1000, eff=eff, messages=self.pop.GetResults().NumberOfMessages, text=txt)

    def close(self):
        self.pop.Close()
