"""Shared OpticStudio machinery: connection, model building, verified POP settings,
convergence test, progress log and small I/O helpers.

The independent analytic model lives in analytic.py and never imports from here.
Its functions are re-exported below only for convenience of the Zemax-side scripts.

Every POP setting comes from a run_config.json ("pop_tokens", "surface_settings").
The binary POP .CFG is generated from it, and every write is read back out of the
CFG; any mismatch raises SettingsMismatch. POPD reads the *saved default* POP
settings, so commit() loads the CFG into the analysis and saves it before any
POPD evaluation.
"""
import csv, glob, math, os, re, struct, tempfile, time
from pathlib import Path

import numpy as np

from analytic import (ONE_DB, loss_db, q_param, w_R, refract, beam_at, overlap,   # noqa: F401 (re-export)
                      tolerance, peak, level_crossing)

ROOT = Path(__file__).resolve().parents[1]
PROGRESS_LOG = ROOT / "results" / "progress.log"
INT_TOKENS = {"POP_WAVE", "POP_FIELD", "POP_START", "POP_END", "POP_BEAMTYPE",
              "POP_SAMPX", "POP_SAMPY", "POP_COMPUTE", "POP_FIBERTYPE"}


class SettingsMismatch(AssertionError):
    """A setting did not read back as written, or POP reported something else."""


# ============================================================================ helpers
def progress(tag, msg):
    """Print and append to results/progress.log so a long run can be inspected in flight."""
    print(msg, flush=True)
    PROGRESS_LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(PROGRESS_LOG, "a") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + "[%s] %s\n" % (tag, msg))


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([("%.9g" % v) if isinstance(v, float) else v for v in r])


def frange(lo, hi, step):
    return [float(v) for v in np.round(np.arange(lo, hi + step / 2, step), 6)]


def convergence(ps, points, after=None, limit_pct=0.5):
    """The project convergence test: grid x2 (half pixel) and grid x2 + window x2 (same pixel).

    points: list of (label, setter) where setter() puts the system in the state to test.
    after(): restores the nominal state between points. Returns per-point % changes in eta.
    """
    s0, w0 = ps.tokens["POP_SAMPX"], ps.tokens["POP_WIDEX"]
    out = []
    for label, setter in points:
        setter()
        etas = {}
        for name, s_, w_ in (("base", s0, w0), ("grid x2", s0 + 1, w0), ("grid x2, window x2", s0 + 1, 2 * w0)):
            for k in ("POP_SAMPX", "POP_SAMPY"):
                ps.set(k, s_)
            for k in ("POP_WIDEX", "POP_WIDEY"):
                ps.set(k, w_)
            etas[name] = ps.popd()[0]
        ch = {k: 100 * (v - etas["base"]) / etas["base"] for k, v in etas.items() if k != "base"}
        out.append(dict(point=label, eta=etas, change_pct=ch, passed=max(abs(v) for v in ch.values()) < limit_pct))
        for k in ("POP_SAMPX", "POP_SAMPY"):
            ps.set(k, s0)
        for k in ("POP_WIDEX", "POP_WIDEY"):
            ps.set(k, w0)
        if after:
            after()
    return out


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
        # string settings (POP_SOURCEFILE): no byte read-back possible; verified by their effect on eta
        for tok, val in cfgd.get("pop_strings", {}).items():
            if not tok.startswith("_"):
                self.st.ModifySettings(self.cfg, tok, val)
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
        bw = g(r"Beam Width X = ([\d.Ee+-]+), Y = ([\d.Ee+-]+)")
        eff = [float(x) for x in g(r"Fiber Efficiency: System ([\d.Ee+-]+), Receiver ([\d.Ee+-]+), Coupling ([\d.Ee+-]+)").groups()]
        return dict(grid=n, width_mm=width, pilot_size_um=size * 1000, pilot_waist_um=waist * 1000,
                    pilot_pos_um=pos * 1000, eff=eff, beam_width_um=[float(v) * 1000 for v in bw.groups()] if bw else None, messages=self.pop.GetResults().NumberOfMessages, text=txt)

    def close(self):
        self.pop.Close()
