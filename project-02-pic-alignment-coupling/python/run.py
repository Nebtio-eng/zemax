"""Single entry point: reproduce any stage of Project 02 from its run_config.json.

    python run.py 3            # null test + A0 nominal, checked against the committed values
    python run.py 5            # A0 lateral sweep vs gap              -> results/coupling_curves/
    python run.py 6            # A0 X-Y tolerance map                 -> results/tolerance_maps/
    python run.py 7            # build A1 and B, all 1-D sweeps        -> results/coupling_curves/
    python run.py 8            # A0 angular/longitudinal, B592, A1/B maps
    python run.py 8b           # FDTD grating beam in B (first: python fdtd_source.py, Lumerical) -> results/fdtd/
    python run.py 9            # parameter study on B                  -> results/parameter_studies/
    python run.py 10           # Pareto front                          -> results/optimization/
    python run.py all          # 3, 5, 6, 7, 8, 9, 10 in order (8b is separate: it needs Lumerical)
    add --mode extension to drive an open OpticStudio (Interactive Extension armed);
    the default is a headless standalone instance (never run two at once).

Inside an already-connected session:  import run; run.stage("7", app, zos)
Tests: python test_validated.py   (analytic always; Zemax checks with ZEMAX_MODE=standalone)
"""
import argparse, json

import coupling_analysis as C


def stage3(app, zos):
    """Null test, its two controls, and the A0 nominal POPD, against the Stage 3 record."""
    from coupling_analysis import ROOT, PopSession
    cfgd = json.loads((ROOT / "zemax" / "baseline" / "run_config.json").read_text())
    ps = PopSession(app, cfgd, zos)
    w = cfgd["pop_tokens"]["POP_PARAM1"]
    out = {}
    try:
        ps.set_gap_um(0.0)
        out["null_identical_modes"] = ps.popd()[0]
        for k in ("POP_FPARAM1", "POP_FPARAM2"):
            ps.set(k, 2 * w)
        out["control_receiver_as_diameter"] = ps.popd()[0]
        for k in ("POP_FPARAM1", "POP_FPARAM2"):
            ps.set(k, w)
        ps.set("POP_FPARAM3", w)
        out["control_decentred_by_w"] = ps.popd()[0]
        ps.set("POP_FPARAM3", 0.0)
        ps.set_gap_um(cfgd["sweep_lateral"]["nominal_gap_um"])
        out["A0_popd"] = ps.popd()
    finally:
        ps.set_gap_um(cfgd["sweep_lateral"]["nominal_gap_um"])
        ps.set("POP_FPARAM3", 0.0)
        ps.commit()
        ps.close()
    ref = cfgd["null_test"]
    checks = [("null", out["null_identical_modes"], ref["identical_modes_zero_distance_eta"], 1e-6),
              ("control 1", out["control_receiver_as_diameter"], ref["control_receiver_entered_as_diameter_eta"], 1e-6),
              ("control 2", out["control_decentred_by_w"], ref["control_receiver_decentred_by_w_eta"], 1e-6),
              ("A0 eta", out["A0_popd"][0], cfgd["popd_raw"]["data0_total_eta"], 1e-9)]
    for name, got, want, tol in checks:
        C.progress("stage3", "%-10s %.9f (record %.9f) %s" % (name, got, want, "ok" if abs(got - want) <= tol else "MISMATCH"))
        if abs(got - want) > tol:
            raise AssertionError("Stage 3 %s not reproduced: %.9f vs %.9f" % (name, got, want))
    return out


def stage(name, app, zos):
    import alignment_sweep as A
    import tolerance_map as T
    if name == "3":
        return stage3(app, zos)
    if name == "5":
        return A.run(app, zos=zos)
    if name == "6":
        return T.run_map(app, zos)
    if name == "7":
        return A.run_stage7(app, zos)
    if name == "8":
        r = A.run_stage8(app, zos)
        for cid, ref in (("B", 8.1399), ("A1", 6.1875)):
            T.run_map(app, zos, C.ROOT / "zemax" / "microlens" / cid / "run_config.json", tag=cid, ref_radius_um=ref,
                      half=12, step=1.0, n_angles=16, r_max=12, r_step=0.5)
        return r
    if name == "8b":
        import fdtd_coupling as FC
        return FC.run_all(app, zos, "grating_20p")
    if name == "9":
        import parameter_study as P
        return P.run_all(app, zos)
    if name == "10":
        import optimization as O
        return O.run(app, zos)
    raise ValueError("unknown stage %r" % name)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=("3", "5", "6", "7", "8", "8b", "9", "10", "all"))
    ap.add_argument("--mode", choices=("standalone", "extension"), default="standalone")
    args = ap.parse_args()
    zos = C.load_zosapi()
    conn, app = C.connect(args.mode, zos)
    try:
        for s in (("3", "5", "6", "7", "8", "9", "10") if args.stage == "all" else (args.stage,)):
            C.progress("run", "=== stage %s ===" % s)
            stage(s, app, zos)
    finally:
        if args.mode == "standalone":
            app.CloseApplication()


if __name__ == "__main__":
    main()
