"""Regression tests for the project's validated numbers.

    python test_validated.py                         # analytic tests only (no OpticStudio)
    ZEMAX_MODE=standalone python test_validated.py   # also the OpticStudio checks (headless)

Inside an already-connected session:  import test_validated as T; T.zemax_checks(app, zos)
"""
import json, math, os, unittest

import analytic as AN

# Values validated in Stages 3-9 (docs/validation.md)
NULL = 1.000000
CONTROL_DIAMETER = 0.640000
CONTROL_DECENTRED = math.exp(-1)
A0_ETA = 0.962617754            # Stage 3, POPD 0
A0_LATERAL = 2.249773           # Stage 5, um
B_LATERAL = 8.1399              # Stage 7, um
B_ANGULAR = 0.6758              # Stage 7, deg
INVARIANT = 0.0960              # um rad


class AnalyticModel(unittest.TestCase):
    """The independent model on its own: no OpticStudio involved."""

    def test_a0_eta_closed_form(self):
        self.assertAlmostEqual(AN.a0_closed_form(20.0)["eta0"], A0_ETA, places=8)

    def test_a0_lateral_closed_form(self):
        self.assertAlmostEqual(AN.a0_closed_form(20.0)["d_exact"], A0_LATERAL, places=5)

    def test_tolerance_extractor(self):
        self.assertAlmostEqual(AN.selftest(), A0_LATERAL, places=4)

    def test_b_prediction(self):
        p = AN.gauss_prediction()
        self.assertAlmostEqual(p["lateral_rel"], B_LATERAL, delta=2e-3)
        self.assertAlmostEqual(p["angular_rel"], B_ANGULAR, delta=2e-3)

    def test_invariant_floor(self):
        self.assertAlmostEqual(AN.IDEAL_PRODUCT, INVARIANT, delta=1e-4)
        self.assertAlmostEqual(AN.gauss_prediction()["product_rel"], AN.IDEAL_PRODUCT, delta=2e-5)

    def test_collimation_pair(self):
        self.assertAlmostEqual(AN.r_match_um(630.0), 486.06, delta=0.05)

    def test_independence(self):
        """analytic.py must not import from the Zemax path."""
        import ast, pathlib
        src = (pathlib.Path(__file__).parent / "analytic.py").read_text()
        mods = {(n.names[0].name if isinstance(n, ast.Import) else n.module)
                for n in ast.walk(ast.parse(src)) if isinstance(n, (ast.Import, ast.ImportFrom))}
        self.assertTrue(mods <= {"math", "numpy", "scipy.interpolate", "scipy.optimize"}, mods)


def zemax_checks(app, zos):
    """The validated OpticStudio numbers. Raises AssertionError on any mismatch."""
    import run
    import parameter_study as P
    r3 = run.stage3(app, zos)                       # null test, both controls, A0 eta (asserts inside)
    assert abs(r3["null_identical_modes"] - NULL) < 1e-6
    assert abs(r3["control_receiver_as_diameter"] - CONTROL_DIAMETER) < 1e-6
    assert abs(r3["control_decentred_by_w"] - CONTROL_DECENTRED) < 1e-6
    assert abs(r3["A0_popd"][0] - A0_ETA) < 1e-9
    b = P.Bench(app, zos)
    try:
        m = b.metrics(guess=dict(lateral=B_LATERAL, angular=B_ANGULAR, longitudinal=713.3))
    finally:
        b.close()
    assert abs(m["lateral_rel"] - B_LATERAL) < 1e-3, m["lateral_rel"]
    assert abs(m["angular_rel"] - B_ANGULAR) < 1e-3, m["angular_rel"]
    assert abs(m["product_rel"] - INVARIANT) < 1e-4, m["product_rel"]
    return dict(stage3=r3, B=dict(lateral=m["lateral_rel"], angular=m["angular_rel"], product=m["product_rel"]))


@unittest.skipUnless(os.environ.get("ZEMAX_MODE"), "set ZEMAX_MODE=standalone|extension to run OpticStudio checks")
class OpticStudio(unittest.TestCase):

    def test_validated_numbers(self):
        import coupling_analysis as C
        zos = C.load_zosapi()
        conn, app = C.connect(os.environ["ZEMAX_MODE"], zos)
        try:
            print(json.dumps(zemax_checks(app, zos), indent=1, default=str))
        finally:
            if os.environ["ZEMAX_MODE"] == "standalone":
                app.CloseApplication()


if __name__ == "__main__":
    unittest.main(verbosity=2)
