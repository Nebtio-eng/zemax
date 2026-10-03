"""2-D X-Y lateral tolerance map and enclosed 1-dB area A_1dB.

For a circular source and receiver the 1-dB contour must be a circle. The radius
is measured precisely along radial cuts (spline-interpolated), and the area is
computed two ways: from the radial cuts (integral of r^2/2 over angle) and from
the contour of an X-Y grid. Both conventions are reported: rel (1 dB above the
on-axis loss) and abs (total loss 1.000 dB; zero area if on-axis loss >= 1 dB).

Stage 6: A0 at 20 um (machinery check).   Stage 8: A1 and B.

Usage:
    python tolerance_map.py --config zemax/microlens/B/run_config.json --mode standalone
"""
import argparse, csv, json, math, time
from pathlib import Path

import numpy as np

import analytic as AN
import coupling_analysis as C
from analytic import loss_db
from coupling_analysis import ROOT, PopSession

A0_CONFIG = ROOT / "zemax" / "baseline" / "run_config.json"
OUT = ROOT / "results" / "tolerance_maps"
CIRCULARITY_LIMIT = 0.005         # r_max/r_min - 1 above this means something is broken
STAGE5_TOL_UM = 2.249773          # A0, 20 um gap (results/coupling_curves/a0_tolerance_vs_gap.csv)


def nominal_gap(cfgd):
    return cfgd["parameters"]["gap_um"] if "parameters" in cfgd else cfgd["sweep_lateral"]["nominal_gap_um"]


def eta_at(ps, x_um, y_um):
    ps.set("POP_FPARAM3", x_um / 1000)
    ps.set("POP_FPARAM4", y_um / 1000)
    return ps.popd()


def cut_area(r):
    """Area of a star-shaped region from equally spaced radial cuts: integral of r^2/2."""
    r = np.nan_to_num(np.asarray(r, float))
    return float(0.5 * np.sum(r ** 2) * (2 * math.pi / len(r)))


def run_map(app, zos, cfg_path=A0_CONFIG, tag="A0", ref_radius_um=STAGE5_TOL_UM,
            half=3.5, step=0.25, n_angles=24, r_max=3.2, r_step=0.1):
    cfgd = json.loads(Path(cfg_path).read_text())
    gap_um = nominal_gap(cfgd)
    OUT.mkdir(parents=True, exist_ok=True)
    ps = PopSession(app, cfgd, zos)
    t0 = time.time()
    try:
        ps.set_gap_um(gap_um)
        e0 = eta_at(ps, 0.0, 0.0)
        l0 = float(loss_db(e0[0]))
        radii = [float(r) for r in np.round(np.arange(0, r_max + r_step / 2, r_step), 6)]
        cuts = []
        for k in range(n_angles):
            phi = 2 * math.pi * k / n_angles
            rows = [(0.0, *e0, l0)]
            for r in radii[1:]:
                e = eta_at(ps, r * math.cos(phi), r * math.sin(phi))
                rows.append((r, *e, float(loss_db(e[0]))))
            arr = np.array(rows)
            t = AN.tolerance(arr[:, 0], arr[:, 4], 0.0)
            cuts.append((math.degrees(phi), t["rel_plus"], t["abs_plus"]))
        xs = [float(v) for v in np.round(np.arange(-half, half + step / 2, step), 6)]
        grid = np.zeros((len(xs), len(xs)))
        for i, y in enumerate(xs):
            for j, x in enumerate(xs):
                grid[i, j] = float(loss_db(eta_at(ps, x, y)[0]))
    finally:
        ps.set("POP_FPARAM3", 0.0)
        ps.set("POP_FPARAM4", 0.0)
        ps.set_gap_um(gap_um)
        ps.commit()
        nchecks = ps.nchecks
        ps.close()

    r = np.array([c[1] for c in cuts])
    ra = np.array([c[2] for c in cuts])
    circularity = float(r.max() / r.min() - 1)
    by_angle = {round(c[0], 3): c for c in cuts}
    res = dict(config=tag, gap_um=gap_um, eta0=e0[0], L0_dB=l0,
               r_rel_mean_um=float(r.mean()), r_rel_min_um=float(r.min()), r_rel_max_um=float(r.max()),
               x_rel_um=by_angle[0.0][1], y_rel_um=by_angle[90.0][1] if 90.0 in by_angle else float("nan"),
               x_abs_um=by_angle[0.0][2], y_abs_um=by_angle[90.0][2] if 90.0 in by_angle else float("nan"),
               circularity=circularity, circular=circularity < CIRCULARITY_LIMIT,
               A1dB_rel_from_cuts_um2=cut_area(r), A1dB_abs_from_cuts_um2=cut_area(ra),
               A1dB_rel_from_grid_contour_um2=contour_area(xs, grid, l0 + 1.0),
               A1dB_abs_from_grid_contour_um2=contour_area(xs, grid, 1.0) if l0 < 1.0 else 0.0,
               reference_radius_um=ref_radius_um, pi_r2_reference_um2=math.pi * ref_radius_um ** 2)
    res["area_rel_vs_pi_r2_pct"] = 100 * (res["A1dB_rel_from_cuts_um2"] / res["pi_r2_reference_um2"] - 1)
    with open(OUT / ("%s_radial_cuts.csv" % tag), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["angle_deg", "d1dB_rel_um", "d1dB_abs_um"])
        w.writerows([["%.1f" % a, "%.6f" % b, "%.6f" % c] for a, b, c in cuts])
    with open(OUT / ("%s_xy_loss_map.csv" % tag), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["y_um \\ x_um"] + ["%.2f" % x for x in xs])
        for y, row in zip(xs, grid):
            w.writerow(["%.2f" % y] + ["%.6f" % v for v in row])
    plot_map(xs, grid, l0, cuts, res, OUT / ("%s_xy_tolerance_map.png" % tag))
    (OUT / ("run_config_%s_map.json" % tag)).write_text(json.dumps(dict(
        configuration=cfgd.get("configuration", tag), source_config=str(Path(cfg_path).resolve().relative_to(ROOT)),
        opticstudio=dict(build=str(app.OpticStudioVersion), license=str(app.LicenseStatus), mode=str(app.Mode)),
        pop_tokens=cfgd["pop_tokens"], surface_settings=cfgd.get("surface_settings"),
        grid=dict(half_width_um=half, step_um=step, points=len(xs) ** 2),
        radial_cuts=dict(n_angles=n_angles, r_max_um=r_max, r_step_um=r_step),
        settings_readbacks_passed=nchecks, results=res, wall_time_s=time.time() - t0), indent=2))
    return res


def contour_area(xs, grid, level):
    """Area enclosed by the `level` contour of the grid (shoelace on matplotlib's contour path)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    cs = ax.contour(xs, xs, grid, levels=[level])
    paths = cs.allsegs[0]
    plt.close(fig)
    if not paths:
        return 0.0
    seg = max(paths, key=len)
    x, y = seg[:, 0], seg[:, 1]
    return float(0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))))


def plot_map(xs, grid, l0, cuts, res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2 = "#0b0b0b", "#52514e"
    fig, ax = plt.subplots(figsize=(5.6, 4.8), dpi=150)
    im = ax.imshow(grid - l0, extent=[xs[0], xs[-1], xs[0], xs[-1]], origin="lower", cmap="Blues", vmin=0, vmax=3)
    ax.contour(xs, xs, grid - l0, levels=[1.0], colors=["#eb6834"], linewidths=2)
    th = np.linspace(0, 2 * math.pi, 300)
    rr = res["reference_radius_um"]
    ax.plot(rr * np.cos(th), rr * np.sin(th), color=ink, lw=1, ls="--")
    ax.plot([c[1] * math.cos(math.radians(c[0])) for c in cuts], [c[1] * math.sin(math.radians(c[0])) for c in cuts],
            lw=0, marker="o", ms=4, color=ink)
    ax.set_xlabel("Receiver X decenter (um)", color=ink)
    ax.set_ylabel("Receiver Y decenter (um)", color=ink)
    ax.set_title("%s: 1-dB contour (orange) vs circle r = %.3f um (dashed)" % (res["config"], rr),
                 loc="left", fontsize=9, color=ink)
    ax.set_aspect("equal")
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("Loss above on-axis (dB)", color=ink2)
    ax.text(0.02, 0.02, "A_1dB (rel) = %.2f um^2 (pi r^2 = %.2f)\ncircularity %.1e" % (
        res["A1dB_rel_from_cuts_um2"], res["pi_r2_reference_um2"], res["circularity"]),
        transform=ax.transAxes, fontsize=8, color=ink, va="bottom")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(A0_CONFIG))
    ap.add_argument("--mode", choices=("standalone", "extension"), default="standalone")
    args = ap.parse_args()
    zos = C.load_zosapi()
    conn, app = C.connect(args.mode, zos)
    cfgd = json.loads(Path(args.config).read_text())
    try:
        if cfgd.get("id", "A0") == "A0":
            res = run_map(app, zos)
        else:
            res = run_map(app, zos, args.config, tag=cfgd["id"], ref_radius_um=float("nan"),
                          half=12, step=1.0, n_angles=16, r_max=12, r_step=0.5)
        print(json.dumps(res, indent=2))
    finally:
        if args.mode == "standalone":
            app.CloseApplication()


if __name__ == "__main__":
    main()
