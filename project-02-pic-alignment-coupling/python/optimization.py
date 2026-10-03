"""Stage 10: robust optimisation and Pareto front for the backside micro-lens (family of B).

Two-dimensional by design. Stage 9 showed conic, aperture (above its threshold),
gap, incidence angle and material flat, so the only real freedoms are:

  t    substrate thickness, with the receiver MFD matched to the beam (2 w(t)):
       chooses WHERE on the lateral/angular invariant curve the design sits
  rho  lens ROC / collimating ROC for that thickness:
       chooses how much LOSS is paid for WORKING DISTANCE

Each grid point: predicted (ABCD + overlap), then measured in POP (tolerances by
root-finding on POPD, both conventions), invariant product recorded, plus the
worst-case loss under a +/-5% ROC error (ASSUMED manufacturing tolerance).

Objectives (Pareto, absolute convention): minimise loss; maximise lateral,
angular, longitudinal 1-dB and area A_1dB = pi * lateral^2 (every B-family
contour was measured circular to 1e-9, Stage 8).

Validation gate: the Stage 9 point (630 um, ROC 420 um, 34 um fibre: 0.30 dB,
965 um working distance) must be reproduced, and the front must contain the
ROC-below-collimation branch. If not, the search is declared inadequate.

Usage:  python optimization.py --mode standalone
"""
import argparse, csv, json, math, time
from pathlib import Path

import numpy as np

import analytic as AN
import coupling_analysis as C
import parameter_study as P
from coupling_analysis import ROOT

OUT = ROOT / "results" / "optimization"
T_GRID = [400.0, 500.0, 630.0, 700.0, 786.0, 800.0, 900.0, 1000.0, 1200.0]
RHO_GRID = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05]
ROC_ERR = 0.05                      # ASSUMED +/-5% lens ROC manufacturing tolerance, for the robustness metric
PASSIVE_UM = 10.0                   # passive-alignment target, +/- um lateral (prompt; literature ~+/-10 um)
EPS_REL, EPS_LOSS_DB = 0.01, 0.01   # epsilon-dominance: 1% on tolerances, 0.01 dB on loss


def log(m):
    C.progress("stage10", m)


def design(t, rho, mfd=None):
    return dict(t_um=t, rho=rho, roc_um=rho * AN.r_match_um(t), mfd_um=mfd or 2 * AN.w_exit_um(t))


def apply(b, d):
    b.set_thickness_um(d["t_um"])
    b.set_radius_um(-d["roc_um"])
    b.set_receiver_mfd(d["mfd_um"])


def evaluate(b, d):
    b.restore()
    apply(b, d)
    p = AN.gauss_prediction(d["t_um"], -d["roc_um"], mfd_um=d["mfd_um"])
    guess = {k: p[k + "_rel"] for k in ("lateral", "angular", "longitudinal") if not math.isnan(p[k + "_rel"])}
    m = b.metrics(guess=guess)
    # robustness: worst nominal loss under a +/-5% ROC error (thickness and fibre unchanged)
    rob = []
    for s in (1 - ROC_ERR, 1 + ROC_ERR):
        b.set_radius_um(-d["roc_um"] * s)
        rob.append(b.loss())
    b.set_radius_um(-d["roc_um"])
    row = dict(d, **{k: m[k] for k in ("loss_dB", "eta", "S", "T", "lateral_rel", "lateral_abs", "angular_rel",
                                         "angular_abs", "longitudinal_rel", "longitudinal_abs", "product_rel", "product_abs")})
    row.update(area_rel_um2=math.pi * m["lateral_rel"] ** 2,
               area_abs_um2=(math.pi * m["lateral_abs"] ** 2) if not math.isnan(m["lateral_abs"]) else 0.0,
               product_excess_pct=100 * (m["product_rel"] / AN.IDEAL_PRODUCT - 1),
               worst_loss_roc5_dB=max(rob), pred_loss_dB=p["loss_dB"], pred_lateral_abs=p["lateral_abs"],
               pred_angular_abs=p["angular_abs"], pred_longitudinal_abs=p["longitudinal_abs"],
               pred_product_rel=p["product_rel"])
    # window rule: 4x the largest beam met in the sweeps (beam at the far longitudinal point)
    surf = [{"thickness_mm": "inf"}, {"thickness_mm": d["t_um"] / 1000}, {"radius_mm": -d["roc_um"] / 1000}, {}]
    zfar = 20 + (m["longitudinal_rel"] if not math.isnan(m["longitudinal_rel"]) else 0)
    wfar = AN.beam_at(surf, zfar, {0: 1.0, 1: AN.N_SI, 2: 1.0}, 4.6, AN.LAM, 2)[0]
    row["largest_beam_um"] = max(wfar, AN.w_exit_um(d["t_um"]), d["mfd_um"] / 2)
    row["window_ok"] = b.ps.tokens["POP_WIDEX"] * 1000 >= 4 * row["largest_beam_um"]
    return row


def objectives(r):
    """Loss (min), then lateral, angular, longitudinal, area (max); absolute convention, nan -> 0."""
    z = lambda v: 0.0 if v is None or (isinstance(v, float) and math.isnan(v)) else v
    return (r["loss_dB"], z(r["lateral_abs"]), z(r["angular_abs"]), z(r["longitudinal_abs"]), z(r["area_abs_um2"]))


def dominates(a, b, eps_rel=0.0, eps_loss=0.0):
    fa, fb = objectives(a), objectives(b)
    ge = [fa[0] <= fb[0] + eps_loss] + [x >= y - eps_rel * abs(y) for x, y in zip(fa[1:], fb[1:])]
    gt = [fa[0] < fb[0] - eps_loss] + [x > y + eps_rel * abs(y) for x, y in zip(fa[1:], fb[1:])]
    return all(ge) and any(gt)


def pareto(rows, eps_rel=0.0, eps_loss=0.0):
    return [r for r in rows if not any(dominates(o, r, eps_rel, eps_loss) for o in rows if o is not r)]


def run(app, zos):
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    b = P.Bench(app, zos)
    rows, conv = [], []
    try:
        # --- gate part 1: reproduce the Stage 9 below-collimation point and B itself (34 um fibre) ---
        gate_pts = [dict(design(630.0, 420.0 / AN.r_match_um(630.0), 34.0), role="gate: Stage 9 ROC 420 um"),
                    dict(design(630.0, 480.0 / AN.r_match_um(630.0), 34.0), role="gate: B as built")]
        for d in gate_pts:
            r = evaluate(b, d)
            r["role"] = d["role"]
            rows.append(r)
            log("%s: loss %.4f dB, long rel %.1f um, lat abs %.3f um" % (d["role"], r["loss_dB"], r["longitudinal_rel"], r["lateral_abs"]))
        g = rows[0]
        if abs(g["loss_dB"] - 0.2975) > 0.005 or abs(g["longitudinal_rel"] - 964.8) > 2.0:
            raise RuntimeError("gate failed: Stage 9 ROC-420 point not reproduced (%.4f dB, %.1f um)" % (g["loss_dB"], g["longitudinal_rel"]))
        # --- grid ---
        for t in T_GRID:
            for rho in RHO_GRID:
                r = evaluate(b, design(t, rho))
                r["role"] = "grid"
                rows.append(r)
                log("t %6.1f rho %.2f (ROC %6.1f, MFD %5.1f): loss %.4f | lat %.3f/%.3f | ang %.4f/%.4f | long %.0f/%.0f | "
                    "prod %.5f (%+.2f%%) | worst +/-5%% ROC %.3f dB | window %s  [%.0fs]" % (
                        t, rho, r["roc_um"], r["mfd_um"], r["loss_dB"], r["lateral_rel"], r["lateral_abs"],
                        r["angular_rel"], r["angular_abs"], r["longitudinal_rel"], r["longitudinal_abs"],
                        r["product_rel"], r["product_excess_pct"], r["worst_loss_roc5_dB"],
                        "ok" if r["window_ok"] else "TOO SMALL", time.time() - t0))
                if not r["window_ok"]:
                    raise RuntimeError("POP window smaller than 4x the largest beam at t %g rho %g" % (t, rho))
        # --- convergence at the extremes of the grid ---
        for t, rho in ((400.0, 0.80), (1200.0, 0.80), (1200.0, 1.05)):
            b.restore()
            apply(b, design(t, rho))
            m = [r for r in rows if r["t_um"] == t and r["rho"] == rho and r["role"] == "grid"][0]
            c = b.convergence("t %g rho %g" % (t, rho), extra=dict(x=m["lateral_rel"]))
            conv.append(c)
            if not c["passed"]:
                raise RuntimeError("convergence failed at t %g rho %g" % (t, rho))
        b.restore()
    finally:
        b.close()

    grid = [r for r in rows if r["role"] == "grid"]
    strict = pareto(rows)
    front = pareto(rows, EPS_REL, EPS_LOSS_DB)
    for r in rows:
        r["pareto_strict"] = any(r is s for s in strict)
        r["pareto_eps"] = any(r is s for s in front)
    # --- gate part 2: the below-collimation branch must be on the front ---
    branch = [r for r in front if r["rho"] < 0.999]
    longer = [r for r in branch if r["longitudinal_abs"] > max(
        (q["longitudinal_abs"] for q in grid if q["t_um"] == r["t_um"] and abs(q["rho"] - 1) < 1e-9), default=0)]
    gate = dict(stage9_point_reproduced=True, below_collimation_on_front=len(longer) > 0,
                n_front=len(front), n_front_below_collimation=len(branch),
                above_collimation_on_front=[(r["t_um"], r["rho"]) for r in front if r["rho"] > 1.001])
    log("gate: %s" % gate)
    if not gate["below_collimation_on_front"]:
        raise RuntimeError("gate failed: the front contains no ROC-below-collimation design; the search is inadequate")
    # --- passive alignment ---
    passive = sorted([r for r in grid if abs(r["rho"] - 1) < 1e-9], key=lambda r: r["t_um"])
    reach = [r for r in passive if r["lateral_abs"] >= PASSIVE_UM]
    write(rows, conv, gate, reach, time.time() - t0)
    plots(rows, front)
    return rows, front, gate


def write(rows, conv, gate, reach, wall):
    keys = ["role", "t_um", "rho", "roc_um", "mfd_um", "loss_dB", "S", "T", "lateral_rel", "lateral_abs", "angular_rel",
            "angular_abs", "longitudinal_rel", "longitudinal_abs", "area_rel_um2", "area_abs_um2", "product_rel",
            "product_abs", "product_excess_pct", "worst_loss_roc5_dB", "pareto_strict", "pareto_eps", "pred_loss_dB",
            "pred_lateral_abs", "pred_angular_abs", "pred_longitudinal_abs", "pred_product_rel", "largest_beam_um", "window_ok"]
    with open(OUT / "design_grid.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (("%.6g" % v) if isinstance(v, float) else v) for k, v in r.items() if k in keys})
    with open(OUT / "pareto_front.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in sorted([r for r in rows if r["pareto_eps"]], key=lambda r: (r["t_um"], r["rho"])):
            w.writerow({k: (("%.6g" % v) if isinstance(v, float) else v) for k, v in r.items() if k in keys})
    (OUT / "run_config.json").write_text(json.dumps(dict(
        stage="10", base_config="zemax/microlens/B/run_config.json",
        free_parameters=dict(t_um=T_GRID, rho=RHO_GRID, receiver="MFD = 2 w(t), matched"),
        held_flat_from_stage9=dict(conic=0, aperture="none (above threshold)", gap_um=20, incidence_deg=0,
                                   material="SILICON_1310"),
        objectives="min loss; max lateral, angular, longitudinal 1-dB and area (absolute convention)",
        pareto=dict(strict="plain dominance", eps=dict(rel_tolerance=EPS_REL, loss_dB=EPS_LOSS_DB)),
        robustness=dict(roc_error=ROC_ERR, provenance="ASSUMED; no manufacturing tolerance reported in the sources"),
        passive_alignment=dict(target_um=PASSIVE_UM, reached_by=[(r["t_um"], round(r["lateral_abs"], 3)) for r in reach]),
        gate=gate, convergence=conv, results=rows, wall_time_s=wall), indent=2, default=str))


def plots(rows, front):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, grid_c = "#0b0b0b", "#52514e", "#e4e3de"
    g = [r for r in rows if r["role"] == "grid"]

    def style(ax):
        ax.grid(color=grid_c, lw=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(colors=ink2, labelsize=8)

    # 1. lateral vs angular, colour = loss, invariant curve
    fig, ax = plt.subplots(figsize=(7, 5), dpi=150)
    d = np.linspace(3, 14, 200)
    ax.plot(d, np.degrees(AN.IDEAL_PRODUCT / d), color=ink2, ls="--", lw=1, label="invariant floor, 0.0960 um rad")
    sc = ax.scatter([r["lateral_abs"] for r in g], [r["angular_abs"] for r in g], c=[r["loss_dB"] for r in g],
                    cmap="Blues", vmin=0, vmax=0.8, s=36, edgecolors=ink2, linewidths=0.4, zorder=3)
    fr = [r for r in front if r["role"] == "grid"]
    ax.scatter([r["lateral_abs"] for r in fr], [r["angular_abs"] for r in fr], facecolors="none", edgecolors="#eb6834",
               s=110, linewidths=1.5, label="Pareto front (1% / 0.01 dB)", zorder=4)
    ax.axvline(PASSIVE_UM, color="#eb6834", lw=1, ls=":")
    ax.annotate("passive alignment, +/-10 um", (PASSIVE_UM, 1.05), fontsize=8, color=ink, ha="left", xytext=(4, 0),
                textcoords="offset points")
    fig.colorbar(sc, ax=ax, shrink=0.85).set_label("Nominal loss (dB)", color=ink2)
    ax.set_xlabel("Lateral 1-dB, absolute (um)", color=ink)
    ax.set_ylabel("Angular 1-dB, absolute (deg)", color=ink)
    ax.set_title("Stage 10: lateral vs angular tolerance of every design", loc="left", fontsize=10, color=ink)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    style(ax)
    fig.tight_layout()
    fig.savefig(OUT / "pareto_lateral_angular.png")
    plt.close(fig)

    # 2. loss vs working distance along rho, per thickness
    fig, ax = plt.subplots(figsize=(7, 4.6), dpi=150)
    cols = {630.0: "#2a78d6", 800.0: "#eb6834", 1000.0: "#1baf7a"}
    for t, col in cols.items():
        s = sorted([r for r in g if r["t_um"] == t], key=lambda r: r["rho"])
        ax.plot([r["longitudinal_rel"] for r in s], [r["loss_dB"] for r in s], color=col, marker="o", ms=4, lw=1.6,
                label="t = %g um" % t)
        for r in s:
            ax.annotate("%.2f" % r["rho"], (r["longitudinal_rel"], r["loss_dB"]), fontsize=7, color=ink2,
                        xytext=(3, 3), textcoords="offset points")
    gp = [r for r in rows if r["role"].startswith("gate: Stage 9")][0]
    ax.scatter([gp["longitudinal_rel"]], [gp["loss_dB"]], marker="*", s=140, color=ink, zorder=5,
               label="Stage 9 point (ROC 420 um)")
    ax.set_xlabel("Working distance: longitudinal 1-dB from the 20 um gap, relative (um)", color=ink)
    ax.set_ylabel("Nominal loss (dB)", color=ink)
    ax.set_title("ROC below collimation buys working distance; above it loses both (labels: rho)",
                 loc="left", fontsize=9, color=ink)
    ax.legend(frameon=False, fontsize=8)
    style(ax)
    fig.tight_layout()
    fig.savefig(OUT / "loss_vs_working_distance.png")
    plt.close(fig)

    # 3. passive alignment: lateral and angular against thickness (two panels, no dual axis)
    s = sorted([r for r in g if abs(r["rho"] - 1) < 1e-9], key=lambda r: r["t_um"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), dpi=150)
    axes[0].plot([r["t_um"] for r in s], [r["lateral_abs"] for r in s], color="#2a78d6", marker="o", ms=4, lw=1.8)
    axes[0].axhline(PASSIVE_UM, color="#eb6834", ls=":", lw=1)
    axes[0].axhline(2.055, color=ink2, ls="--", lw=0.8)
    axes[0].annotate("A0 (2.06 um)", (s[0]["t_um"], 2.055), fontsize=8, color=ink2, xytext=(0, 3), textcoords="offset points")
    axes[0].annotate("passive target +/-10 um", (s[0]["t_um"], PASSIVE_UM), fontsize=8, color=ink, xytext=(0, 3),
                     textcoords="offset points")
    axes[0].set_ylabel("Lateral 1-dB, absolute (um)", color=ink)
    axes[1].plot([r["t_um"] for r in s], [r["angular_abs"] for r in s], color="#2a78d6", marker="o", ms=4, lw=1.8)
    axes[1].set_ylabel("Angular 1-dB, absolute (deg)", color=ink)
    for ax in axes:
        ax.set_xlabel("Substrate thickness (um); ROC collimating, fibre matched", color=ink)
        style(ax)
    fig.suptitle("Passive alignment: lateral and angular tolerance against substrate thickness",
                 x=0.01, ha="left", fontsize=10, color=ink)
    fig.tight_layout()
    fig.savefig(OUT / "passive_alignment.png")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=("standalone", "extension"), default="standalone")
    args = ap.parse_args()
    zos = C.load_zosapi()
    conn, app = C.connect(args.mode, zos)
    try:
        run(app, zos)
    finally:
        if args.mode == "standalone":
            app.CloseApplication()


if __name__ == "__main__":
    main()
