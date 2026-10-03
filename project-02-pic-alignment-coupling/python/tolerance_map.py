"""2-D X-Y lateral tolerance map and enclosed 1-dB area A_1dB.

Stage 6 runs it on A0 at the 20 um gap as a machinery check: with a circular
source and receiver the 1-dB contour must be a circle of radius equal to the
Stage 5 1-D tolerance. The radius is measured precisely along radial cuts
(spline-interpolated), and the area is computed two ways: from the radial cuts
(integral of r^2/2 over angle) and from the contour of the X-Y grid.

Usage:
    python tolerance_map.py --mode standalone
"""
import argparse, csv, json, math, time

import numpy as np

import coupling_analysis as C
from coupling_analysis import ROOT, PopSession, loss_db

CONFIG = ROOT / "zemax" / "baseline" / "run_config.json"
OUT = ROOT / "results" / "tolerance_maps"
STAGE5_TOL_UM = 2.249773          # results/coupling_curves/a0_tolerance_vs_gap.csv, gap 20 um
CIRCULARITY_LIMIT = 0.005         # r_max/r_min - 1 above this means something is broken


def eta_at(ps, x_um, y_um):
    ps.set("POP_FPARAM3", x_um / 1000)
    ps.set("POP_FPARAM4", y_um / 1000)
    return ps.popd()


def run_map(app, zos, cfg_path=CONFIG, gap_um=20.0, half=3.5, step=0.25, n_angles=24, r_max=3.2, r_step=0.1):
    cfgd = json.loads(cfg_path.read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    ps = PopSession(app, cfgd, zos)
    t0 = time.time()
    try:
        ps.set_gap_um(gap_um)
        e0 = eta_at(ps, 0.0, 0.0)
        l0 = float(loss_db(e0[0]))
        # radial cuts
        radii = [float(r) for r in np.round(np.arange(0, r_max + r_step / 2, r_step), 6)]
        cuts = []
        for k in range(n_angles):
            phi = 2 * math.pi * k / n_angles
            rows = [(0.0, *e0, l0)]
            for r in radii[1:]:
                e = eta_at(ps, r * math.cos(phi), r * math.sin(phi))
                rows.append((r, *e, float(loss_db(e[0]))))
            arr = np.array(rows)
            t = C.tolerance(arr[:, 0], arr[:, 4], 0.0)
            cuts.append((math.degrees(phi), t["rel_plus"], t["abs_plus"]))
        # grid
        xs = [float(v) for v in np.round(np.arange(-half, half + step / 2, step), 6)]
        grid = np.zeros((len(xs), len(xs)))
        for i, y in enumerate(xs):
            for j, x in enumerate(xs):
                grid[i, j] = float(loss_db(eta_at(ps, x, y)[0]))
    finally:
        ps.set("POP_FPARAM3", 0.0)
        ps.set("POP_FPARAM4", 0.0)
        ps.set_gap_um(cfgd["sweep_lateral"]["nominal_gap_um"])
        ps.commit()
        nchecks = ps.nchecks
        ps.close()

    r = np.array([c[1] for c in cuts])
    phis = np.radians([c[0] for c in cuts])
    circularity = float(r.max() / r.min() - 1)
    # area from radial cuts: integral of r(phi)^2 / 2, periodic trapezoid
    area_cuts = float(0.5 * np.sum(r ** 2) * (2 * math.pi / len(r)))
    area_contour = contour_area(xs, grid, l0 + 1.0)
    res = dict(gap_um=gap_um, eta0=e0[0], L0_dB=l0, r_mean_um=float(r.mean()), r_min_um=float(r.min()),
               r_max_um=float(r.max()), circularity=circularity, stage5_tolerance_um=STAGE5_TOL_UM,
               A1dB_from_cuts_um2=area_cuts, A1dB_from_grid_contour_um2=area_contour,
               pi_r2_stage5_um2=math.pi * STAGE5_TOL_UM ** 2,
               area_cuts_vs_pi_r2_pct=100 * (area_cuts / (math.pi * STAGE5_TOL_UM ** 2) - 1),
               area_contour_vs_pi_r2_pct=100 * (area_contour / (math.pi * STAGE5_TOL_UM ** 2) - 1),
               circular=circularity < CIRCULARITY_LIMIT)
    with open(OUT / "A0_radial_cuts.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["angle_deg", "d1dB_rel_um", "d1dB_abs_um"])
        w.writerows([["%.1f" % a, "%.6f" % b, "%.6f" % c] for a, b, c in cuts])
    with open(OUT / "A0_xy_loss_map.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["y_um \\ x_um"] + ["%.2f" % x for x in xs])
        for y, row in zip(xs, grid):
            w.writerow(["%.2f" % y] + ["%.6f" % v for v in row])
    plot_map(xs, grid, l0, cuts, res, OUT / "A0_xy_tolerance_map.png")
    (OUT / "run_config_A0_map.json").write_text(json.dumps(dict(
        stage="6", configuration="A0 conventional, 20 um gap", source_config=str(cfg_path.relative_to(ROOT)),
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
    ax.plot(STAGE5_TOL_UM * np.cos(th), STAGE5_TOL_UM * np.sin(th), color=ink, lw=1, ls="--")
    ax.plot([c[1] * math.cos(math.radians(c[0])) for c in cuts], [c[1] * math.sin(math.radians(c[0])) for c in cuts],
            lw=0, marker="o", ms=4, color=ink)
    ax.set_xlabel("Receiver X decenter (um)", color=ink)
    ax.set_ylabel("Receiver Y decenter (um)", color=ink)
    ax.set_title("A0, 20 um gap: 1-dB contour (orange) vs circle r = %.3f um (dashed)" % STAGE5_TOL_UM,
                 loc="left", fontsize=9, color=ink)
    ax.set_aspect("equal")
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("Loss above on-axis (dB)", color=ink2)
    ax.text(0.02, 0.02, "A_1dB = %.3f um^2 (pi r^2 = %.3f)\ncircularity %.1e" % (
        res["A1dB_from_cuts_um2"], res["pi_r2_stage5_um2"], res["circularity"]),
        transform=ax.transAxes, fontsize=8, color=ink, va="bottom")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=("standalone", "extension"), default="standalone")
    args = ap.parse_args()
    zos = C.load_zosapi()
    conn, app = C.connect(args.mode, zos)
    try:
        print(json.dumps(run_map(app, zos), indent=2))
    finally:
        if args.mode == "standalone":
            app.CloseApplication()


if __name__ == "__main__":
    main()
