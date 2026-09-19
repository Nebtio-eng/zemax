# ============================================================
# ZEMAX LED COLLIMATOR & UNIFORM ILLUMINATION SYSTEM
# Corrected Master's-level optimization workflow
#
# Single Python file:
#   zemax_led_optimizer.py
#
# The script:
#   1. Loads baseline_led_lens.zos without saving it.
#   2. Captures baseline single-lens results.
#   3. Reconstructs the documented two-lens architecture in memory.
#   4. Reconstructs the previous flawed 4.02% candidate.
#   5. Uses fixed-region CV; active-only CV is diagnostic only.
#   6. Applies TARGET_COVERAGE as a hard feasibility floor.
#   7. Audits configured wavelength channels and runs a source/lens
#      capture-fraction diagnostic using ZRD data.
#   8. Runs a bare-source detector-Z sweep plus a 10x-ray sanity check.
#   9. Runs 150-point LHS coarse search if capture is adequate.
#  10. Refines the top 10 coarse candidates with 150 multi-start
#      fine evaluations.
#  11. Computes a Pareto front over coverage, fixed-region CV,
#      and centroid error.
#  12. Validates the top 3 at full 100x100 resolution.
#  13. Saves detector maps, CSVs, and a before/after comparison.
#
# IMPORTANT:
#   - The uploaded baseline file is never saved or overwritten.
#   - The two-lens architecture is constructed in memory from the
#     uploaded single-lens baseline if object #4 is absent.
#   - All search evaluations use the same fixed number of source
#     analysis rays and the same NSC random seed where supported.
#   - Coverage peak estimation uses a 3x3 box-smoothed detector map
#     to reduce single-pixel Monte Carlo shot-noise sensitivity.
# ============================================================

import os
import csv
import math
import traceback
import textwrap
from pathlib import Path

import numpy as np
import clr


# ============================================================
# USER SETTINGS
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
ZOS_FILE_ENV = os.environ.get("ZEMAX_LED_BASELINE")
ZOS_FILE = str(
    Path(ZOS_FILE_ENV).expanduser().resolve()
    if ZOS_FILE_ENV
    else REPO_ROOT / "models" / "baseline_led_lens.zos"
)

OPTICSTUDIO_PATH = (
    r"C:\Program Files\Ansys Zemax OpticStudio 2024 R1.00"
)

SOURCE_NUMBER = 1
DETECTOR_NUMBER = 2
LENS1_NUMBER = 3
LENS2_NUMBER = 4

DETECTOR_SIZE = 50.0
DETECTOR_PIXELS = 100

# The previous project baseline was approximately 8% coverage.
# This is a FEASIBILITY FLOOR, not a weighted term.
TARGET_COVERAGE = 8.0

# A candidate below the floor is never considered a feasible best.
INFEASIBLE_PENALTY = 1.0e6

# Capture-fraction gate requested for the architecture decision.
CAPTURE_THRESHOLD = 70.0

# Fixed source sampling for reproducibility.
ANALYSIS_RAYS = 50000
SANITY_ANALYSIS_RAYS = 500000  # 10x search-ray count for Monte Carlo peak audit
RANDOM_SEED = 12345

# Search sizes.
COARSE_RUNS = 150
FINE_RUNS = 150
FINE_STARTS = 10
VALIDATION_COUNT = 3

# LHS replacement limit when a geometry-filtered candidate is rejected.
MAX_LHS_REPLACEMENTS = 500

# Capture diagnostic grid.
CAPTURE_L1_Z = [8.0, 16.5, 25.0]
CAPTURE_CLEAR1 = [8.0, 14.0, 20.0]
CAPTURE_DETECTOR_Z = 60.0

# Bare-source detector sweep.
BARE_SOURCE_DETECTOR_Z = [45.0, 60.0, 75.0, 90.0, 100.0]

# Detector map threshold used only to define coverage / active-only CV.
COVERAGE_THRESHOLD_FRACTION = 0.10

# Geometry safety margin.
AXIAL_SAFETY_GAP = 1.0
GEOMETRY_TOL = 1.0e-7

# Search variables: 21 total.
VARIABLES = {
    "L1_Z": (8.0, 25.0),
    "L1_R1": (20.0, 100.0),
    "L1_R2": (-100.0, 100.0),
    "L1_C1": (-1.0, 1.0),
    "L1_C2": (-1.0, 1.0),
    "L1_Clear1": (8.0, 20.0),
    "L1_Edge1": (5.0, 15.0),
    "L1_Clear2": (8.0, 20.0),
    "L1_Edge2": (5.0, 15.0),
    "L1_Thickness": (2.0, 8.0),

    "L2_Z": (25.0, 50.0),
    "L2_R1": (20.0, 100.0),
    "L2_R2": (-100.0, 100.0),
    "L2_C1": (-1.0, 1.0),
    "L2_C2": (-1.0, 1.0),
    "L2_Clear1": (8.0, 20.0),
    "L2_Edge1": (5.0, 15.0),
    "L2_Clear2": (8.0, 20.0),
    "L2_Edge2": (5.0, 15.0),
    "L2_Thickness": (2.0, 8.0),

    "Detector_Z": (45.0, 100.0),
}

VARIABLE_NAMES = list(VARIABLES.keys())

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RESULTS_CSV = OUTPUT_DIR / "optimization_results.csv"
REJECTIONS_CSV = OUTPUT_DIR / "geometry_rejections.csv"
CAPTURE_CSV = OUTPUT_DIR / "capture_diagnostic.csv"
BARE_SOURCE_CSV = OUTPUT_DIR / "bare_source_sweep.csv"
BARE_SOURCE_SANITY_CSV = OUTPUT_DIR / "bare_source_sweep_sanity_10x.csv"
PARETO_CSV = OUTPUT_DIR / "pareto_front.csv"
VALIDATION_CSV = OUTPUT_DIR / "validation_top3.csv"
COMPARISON_CSV = OUTPUT_DIR / "before_after_comparison.csv"
METHODOLOGY_TXT = OUTPUT_DIR / "methodology.txt"


# ============================================================
# IMPORT ZOS-API
# ============================================================

os.add_dll_directory(OPTICSTUDIO_PATH)

clr.AddReference(os.path.join(OPTICSTUDIO_PATH, "ZOSAPI.dll"))
clr.AddReference(
    os.path.join(OPTICSTUDIO_PATH, "ZOSAPI_Interfaces.dll")
)

import ZOSAPI
from System import Enum


# ============================================================
# CONNECTION
# ============================================================

connection = ZOSAPI.ZOSAPI_Connection()
application = connection.CreateNewApplication()

if not Path(ZOS_FILE).is_file():
    raise FileNotFoundError(
        "Baseline Zemax model not found. Expected "
        f"{ZOS_FILE}. The repository includes the model under models/."
    )


if application is None:
    raise RuntimeError("Could not start OpticStudio.")

system = application.PrimarySystem

if system is None:
    raise RuntimeError("Could not obtain PrimarySystem.")

if not system.LoadFile(ZOS_FILE, False):
    raise RuntimeError(f"Could not load {ZOS_FILE}")

nce = system.NCE


# ============================================================
# WAVELENGTH AUDIT
# ============================================================

def audit_wavelength_channels():
    """
    Inspect the system wavelength configuration before any capture
    calculation. Source Ellipse can launch NumberOfAnalysisRays per
    wavelength, so an unexpected second wavelength can explain a ZRD
    record count that is an integer multiple of the requested rays.
    """
    wavelengths = system.SystemData.Wavelengths

    try:
        count = int(wavelengths.NumberOfWavelengths)
    except Exception:
        count = -1

    entries = []

    if count > 0:
        for index in range(1, count + 1):
            try:
                item = wavelengths.GetWavelength(index)
                wavelength_um = float(item.Wavelength)
                weight = float(item.Weight)
                entries.append(
                    {
                        "index": index,
                        "wavelength_um": wavelength_um,
                        "wavelength_nm": 1000.0 * wavelength_um,
                        "weight": weight,
                    }
                )
            except Exception as exc:
                entries.append(
                    {
                        "index": index,
                        "wavelength_um": np.nan,
                        "wavelength_nm": np.nan,
                        "weight": np.nan,
                        "error": str(exc),
                    }
                )

    print("\n" + "=" * 70)
    print("WAVELENGTH CHANNEL AUDIT")
    print("=" * 70)
    print(f"Number of wavelength channels: {count}")

    for item in entries:
        if np.isfinite(item.get("wavelength_nm", np.nan)):
            print(
                f"  #{item['index']}: "
                f"{item['wavelength_nm']:.6f} nm | "
                f"weight={item['weight']:.6g}"
            )
        else:
            print(
                f"  #{item['index']}: "
                f"could not read wavelength details: "
                f"{item.get('error', 'unknown error')}"
            )

    if count > 1:
        print(
            "WARNING: More than one wavelength channel is configured. "
            "NumberOfAnalysisRays may be launched per wavelength, so "
            "a ZRD record count near N x NumberOfAnalysisRays is expected."
        )
    elif count == 1:
        print(
            "One wavelength channel is configured. If ZRD primary records "
            "remain near 2x the requested source-ray count, the excess is "
            "not attributable to a second wavelength channel; the script "
            "therefore reports the ZRD/requested ratio for further API/source "
            "diagnosis without using the requested count as a denominator."
        )

    return count, entries


WAVELENGTH_COUNT, WAVELENGTH_ENTRIES = audit_wavelength_channels()


# ============================================================
# OBJECT ACCESS
# ============================================================

def get_object(number):
    obj = nce.GetObjectAt(number)
    if obj is None:
        raise RuntimeError(f"NCE object #{number} was not found.")
    return obj


source = get_object(SOURCE_NUMBER)
detector = get_object(DETECTOR_NUMBER)
lens1 = get_object(LENS1_NUMBER)


# ============================================================
# NCE PARAMETER ACCESS
# ============================================================

ObjectColumn = ZOSAPI.Editors.NCE.ObjectColumn


def get_cell(obj, parameter):
    column = Enum.Parse(ObjectColumn, "Par" + str(parameter))
    return obj.GetObjectCell(column)


def get_parameter(obj, parameter):
    return float(get_cell(obj, parameter).DoubleValue)


def set_parameter(obj, parameter, value):
    get_cell(obj, parameter).DoubleValue = float(value)


def get_object_type_name(obj):
    try:
        return str(obj.TypeName)
    except Exception:
        try:
            return str(obj.Type)
        except Exception:
            return "UNKNOWN"


# ============================================================
# BASELINE SNAPSHOT
# ============================================================

single_baseline = {
    "lens_z": float(lens1.ZPosition),
    "detector_z": float(detector.ZPosition),
    "material": str(lens1.Material),
}

for p in range(1, 10):
    single_baseline[f"p{p}"] = get_parameter(lens1, p)

try:
    source_data = source.ObjectData
    original_analysis_rays = int(source_data.NumberOfAnalysisRays)
except Exception:
    source_data = None
    original_analysis_rays = -1

try:
    original_source_seed = int(source_data.RandomSeed)
except Exception:
    original_source_seed = -1


# ============================================================
# OUTPUT / PLOTTING
# ============================================================

def save_detector_map(data, filename, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    extent = [
        -DETECTOR_SIZE / 2.0,
        DETECTOR_SIZE / 2.0,
        -DETECTOR_SIZE / 2.0,
        DETECTOR_SIZE / 2.0,
    ]
    im = ax.imshow(
        data,
        origin="lower",
        extent=extent,
        aspect="equal",
    )
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="Irradiance [model units]")
    fig.tight_layout()
    fig.savefig(filename, dpi=180)
    plt.close(fig)


def save_capture_plot(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    zvals = sorted(set(float(r["lens1_z"]) for r in rows))
    for z in zvals:
        sub = [
            r for r in rows
            if abs(float(r["lens1_z"]) - z) < 1e-9
        ]
        sub.sort(key=lambda r: float(r["clear1"]))
        x = [float(r["clear1"]) for r in sub]
        y = [float(r["capture_fraction_pct"]) for r in sub]
        plt.plot(x, y, marker="o", label=f"L1 Z={z:g} mm")

    plt.axhline(CAPTURE_THRESHOLD, linestyle="--")
    plt.xlabel("Lens 1 Clear1 [mm]")
    plt.ylabel("Capture fraction [%]")
    plt.title("Lens-1 source-ray capture diagnostic")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "capture_fraction_vs_clear1.png",
        dpi=180,
    )
    plt.close()


# ============================================================
# SOURCE SAMPLING CONTROL
# ============================================================

def set_source_sampling():
    global source_data

    source_data = source.ObjectData

    source_data.NumberOfAnalysisRays = int(ANALYSIS_RAYS)

    seed_supported = False
    try:
        source_data.RandomSeed = int(RANDOM_SEED)
        seed_supported = True
    except Exception:
        seed_supported = False

    return seed_supported


SOURCE_SEED_SUPPORTED = set_source_sampling()


# ============================================================
# LENS HELPERS
# ============================================================

def configure_standard_lens(obj, z, params, material="N-BK7"):
    obj.ZPosition = float(z)
    obj.Material = material

    for p, key in enumerate(
        [
            "R1", "C1", "Clear1", "Edge1",
            "Thickness", "R2", "C2",
            "Clear2", "Edge2",
        ],
        start=1,
    ):
        set_parameter(obj, p, params[key])


DEFAULT_LENS_PARAMS = {
    "R1": 50.0,
    "C1": 0.0,
    "Clear1": 10.0,
    "Edge1": 10.0,
    "Thickness": 5.0,
    "R2": 0.0,
    "C2": 0.0,
    "Clear2": 10.0,
    "Edge2": 10.0,
}


def make_two_lens_architecture():
    """
    Build the documented two-lens architecture in memory.

    Uploaded baseline:
        #1 Source Ellipse
        #2 Detector Rectangle
        #3 Standard Lens

    If #4 is absent, insert #4 as Standard Lens. The uploaded
    baseline file itself is never saved.
    """
    global lens1, lens2, detector, nce

    if nce.NumberOfObjects < LENS2_NUMBER:
        nce.InsertNewObjectAt(LENS2_NUMBER)

    lens2 = get_object(LENS2_NUMBER)
    detector = get_object(DETECTOR_NUMBER)
    lens1 = get_object(LENS1_NUMBER)

    type_name = get_object_type_name(lens2)

    if "StandardLens" not in type_name and "Standard Lens" not in type_name:
        settings = lens2.GetObjectTypeSettings(
            ZOSAPI.Editors.NCE.ObjectType.StandardLens
        )
        lens2.ChangeType(settings)

    configure_standard_lens(
        lens2,
        30.0,
        DEFAULT_LENS_PARAMS,
        material="N-BK7",
    )

    # Documented two-lens baseline architecture.
    lens1.ZPosition = 15.0
    detector.ZPosition = 60.0

    # Reset lens 1 to the baseline Standard Lens prescription.
    configure_standard_lens(
        lens1,
        15.0,
        DEFAULT_LENS_PARAMS,
        material=single_baseline["material"],
    )


# ============================================================
# DETECTOR DATA
# ============================================================

def get_detector_data():
    """
    Data type 2 = incoherent irradiance in the existing workflow.
    """
    total_pixels = DETECTOR_PIXELS * DETECTOR_PIXELS

    # Use the safe bulk API when available; otherwise fall back to
    # the proven pixel-by-pixel method from the original script.
    try:
        raw = nce.GetAllDetectorDataSafe(
            DETECTOR_NUMBER,
            2,
        )
        arr = np.array(raw, dtype=float)

        if arr.size == total_pixels:
            return arr.reshape(
                DETECTOR_PIXELS,
                DETECTOR_PIXELS,
            )
    except Exception:
        pass

    data = np.zeros(total_pixels)

    for pixel in range(1, total_pixels + 1):
        result = nce.GetDetectorData(
            DETECTOR_NUMBER,
            pixel,
            2,
            0.0,
        )
        data[pixel - 1] = float(result[1])

    return data.reshape(
        DETECTOR_PIXELS,
        DETECTOR_PIXELS,
    )


# ============================================================
# REPRODUCIBLE NSC RAY TRACE
# ============================================================

def run_nsc_ray_trace(
    save_zrd=False,
    zrd_name="capture_diagnostic.ZRD",
):
    ray_trace = system.Tools.OpenNSCRayTrace()

    if ray_trace is None:
        raise RuntimeError("Could not open NSC Ray Trace.")

    try:
        ray_trace.ClearDetectors(0)
        ray_trace.SplitNSCRays = False
        ray_trace.ScatterNSCRays = False
        ray_trace.UsePolarization = False
        ray_trace.IgnoreErrors = True
        ray_trace.SaveRays = bool(save_zrd)

        if save_zrd:
            # OpticStudio saves this in the system directory.
            ray_trace.SaveRaysFile = str(zrd_name)

        # OpticStudio exposes an explicit random-seed API for NSC
        # ray tracing in recent ZOS-API versions.  2024 R1 may not
        # expose it; fall back to the fixed Source Ellipse seed.
        trace_seed_supported = False
        try:
            ray_trace.SetRandomSeed(int(RANDOM_SEED))
            trace_seed_supported = True
        except Exception:
            trace_seed_supported = False

        ray_trace.RunAndWaitForCompletion()

        if not ray_trace.Succeeded:
            raise RuntimeError(
                f"NSC ray trace failed: {ray_trace.ErrorMessage}"
            )

        return {
            "trace_seed_supported": trace_seed_supported,
            "total_ray_energy": float(
                ray_trace.GetTotalRayEnergy()
            ),
        }

    finally:
        ray_trace.Close()


# ============================================================
# METRICS
# ============================================================

def estimate_smoothed_peak(data):
    """
    Estimate the detector peak from a 3x3 box-smoothed irradiance map.

    The raw detector maximum can be dominated by a single Monte Carlo
    ray landing in one pixel. A 3x3 uniform filter suppresses that
    single-pixel shot-noise spike while preserving the beam-scale peak.
    """
    data = np.asarray(data, dtype=float)

    try:
        from scipy.ndimage import uniform_filter

        smoothed = uniform_filter(
            data,
            size=3,
            mode="nearest",
        )
    except Exception:
        # Dependency-free 3x3 uniform-filter fallback.
        padded = np.pad(
            data,
            pad_width=1,
            mode="edge",
        )
        smoothed = (
            padded[:-2, :-2]
            + padded[:-2, 1:-1]
            + padded[:-2, 2:]
            + padded[1:-1, :-2]
            + padded[1:-1, 1:-1]
            + padded[1:-1, 2:]
            + padded[2:, :-2]
            + padded[2:, 1:-1]
            + padded[2:, 2:]
        ) / 9.0

    return float(np.max(smoothed))


def compute_metrics(data):
    data = np.asarray(data, dtype=float)

    peak = estimate_smoothed_peak(data)
    mean = float(np.mean(data))
    std_full = float(np.std(data))

    if peak <= 0.0:
        return {
            "coverage": 0.0,
            "active_cv": 999.0,
            "fixed_cv": 999.0,
            "centroid_x": 999.0,
            "centroid_y": 999.0,
            "centroid_error": 999.0,
            "mean": mean,
            "peak": peak,
            "std_full": std_full,
            "active_pixels": 0,
        }

    threshold = COVERAGE_THRESHOLD_FRACTION * peak
    active = data >= threshold
    active_values = data[active]

    coverage = (
        100.0
        * float(np.count_nonzero(active))
        / float(data.size)
    )

    if active_values.size > 1:
        active_mean = float(np.mean(active_values))
        active_std = float(np.std(active_values))
        active_cv = (
            100.0 * active_std / active_mean
            if active_mean > 0.0
            else 999.0
        )
    else:
        active_cv = 999.0

    # Corrected CV:
    # full fixed 100x100 reference region.
    # This is the metric used for ranking feasible candidates.
    fixed_mean = float(np.mean(data))
    fixed_std = float(np.std(data))
    fixed_cv = (
        100.0 * fixed_std / fixed_mean
        if fixed_mean > 0.0
        else 999.0
    )

    coords = np.linspace(
        -DETECTOR_SIZE / 2.0,
        DETECTOR_SIZE / 2.0,
        DETECTOR_PIXELS,
    )
    X, Y = np.meshgrid(coords, coords)

    total = float(np.sum(data))

    if total > 0.0:
        centroid_x = float(np.sum(X * data) / total)
        centroid_y = float(np.sum(Y * data) / total)
    else:
        centroid_x = 999.0
        centroid_y = 999.0

    centroid_error = math.hypot(
        centroid_x,
        centroid_y,
    )

    return {
        "coverage": coverage,
        "active_cv": active_cv,
        "fixed_cv": fixed_cv,
        "centroid_x": centroid_x,
        "centroid_y": centroid_y,
        "centroid_error": centroid_error,
        "mean": mean,
        "peak": peak,
        "std_full": std_full,
        "active_pixels": int(np.count_nonzero(active)),
    }


def trace_and_measure():
    run_nsc_ray_trace()
    data = get_detector_data()
    return compute_metrics(data), data


# ============================================================
# OBJECTIVE / RANKING
# ============================================================

def rank_record(metrics):
    """
    Coverage gate:
      coverage < TARGET_COVERAGE:
          infeasible, fixed penalty, never a best candidate.

      coverage >= TARGET_COVERAGE:
          rank by fixed-region CV.
          centroid error is the secondary tiebreaker.

    Active-only CV is never used here.
    """
    if metrics["coverage"] < TARGET_COVERAGE:
        return {
            "feasible": False,
            "objective_score": INFEASIBLE_PENALTY,
            "ranking_rule": (
                f"INFEASIBLE: coverage < {TARGET_COVERAGE:.2f}% "
                f"(fixed penalty; excluded from best)"
            ),
            "rank_key": (
                1,
                float("inf"),
                float("inf"),
            ),
        }

    return {
        "feasible": True,
        "objective_score": metrics["fixed_cv"],
        "ranking_rule": (
            "FEASIBLE: fixed-region CV, then centroid error"
        ),
        "rank_key": (
            0,
            metrics["fixed_cv"],
            metrics["centroid_error"],
        ),
    }


# ============================================================
# STANDARD-LENS GEOMETRY VALIDATION
# ============================================================

def sag(radius, conic, r):
    """
    OpticStudio Standard-surface sag:
        z = c*r^2 / (1 + sqrt(1-(1+k)c^2*r^2))

    Radius=0 is the flat-surface special case.
    """
    radius = float(radius)
    conic = float(conic)
    r = float(r)

    if abs(radius) < GEOMETRY_TOL:
        return 0.0

    c = 1.0 / radius
    q = 1.0 - (1.0 + conic) * (c * r) ** 2

    if q < -GEOMETRY_TOL:
        raise ValueError("invalid conic sag domain")

    q = max(q, 0.0)
    denominator = 1.0 + math.sqrt(q)

    if denominator <= GEOMETRY_TOL:
        raise ValueError("singular sag denominator")

    return (c * r * r) / denominator


def validate_lens_geometry(x, prefix):
    reasons = []

    clear1 = abs(float(x[f"{prefix}_Clear1"]))
    clear2 = abs(float(x[f"{prefix}_Clear2"]))
    edge1 = abs(float(x[f"{prefix}_Edge1"]))
    edge2 = abs(float(x[f"{prefix}_Edge2"]))

    r1 = float(x[f"{prefix}_R1"])
    r2 = float(x[f"{prefix}_R2"])
    c1 = float(x[f"{prefix}_C1"])
    c2 = float(x[f"{prefix}_C2"])
    thickness = float(x[f"{prefix}_Thickness"])

    # Explicit restrictions from the Standard Lens definition.
    if edge1 + GEOMETRY_TOL < clear1:
        reasons.append("Edge1<abs(Clear1)")

    if edge2 + GEOMETRY_TOL < clear2:
        reasons.append("Edge2<abs(Clear2)")

    if thickness <= 0.0:
        reasons.append("nonpositive thickness")

    # A flat radius with a nonzero conic constant is not a useful
    # finite-conic representation; force the flat special case.
    if abs(r1) < GEOMETRY_TOL and abs(c1) > GEOMETRY_TOL:
        reasons.append("R1=0 with nonzero C1")

    if abs(r2) < GEOMETRY_TOL and abs(c2) > GEOMETRY_TOL:
        reasons.append("R2=0 with nonzero C2")

    # Validate sag domain at both clear and edge radii.
    for label, radius, conic, clear, edge in [
        ("front", r1, c1, clear1, edge1),
        ("rear", r2, c2, clear2, edge2),
    ]:
        for radius_test, label2 in [
            (clear, "clear"),
            (edge, "edge"),
        ]:
            try:
                sag(radius, conic, radius_test)
            except ValueError as exc:
                reasons.append(
                    f"{label} {label2}: {str(exc)}"
                )

    # Sample the front/rear axial separation over the common
    # radial extent. A negative gap means the surfaces intersect.
    sample_radius = min(
        edge1,
        edge2,
        0.5 * (edge1 + edge2),
    )

    if sample_radius <= 0.0:
        reasons.append("nonpositive common radial extent")
    else:
        rs = np.linspace(
            0.0,
            sample_radius,
            101,
        )

        min_gap = float("inf")

        for rr in rs:
            try:
                z1 = sag(r1, c1, rr)
                z2 = sag(r2, c2, rr)
                gap = thickness + z2 - z1
                min_gap = min(min_gap, gap)
            except ValueError:
                min_gap = -1.0
                break

        if min_gap <= GEOMETRY_TOL:
            reasons.append(
                f"front/back self-intersection "
                f"(minimum axial gap={min_gap:.6g} mm)"
            )

    return reasons


def valid_geometry(x):
    reasons = []

    l2_min = (
        x["L1_Z"]
        + x["L1_Thickness"]
        + AXIAL_SAFETY_GAP
    )

    if x["L2_Z"] <= l2_min:
        reasons.append(
            "L2_Z<=L1_Z+L1_Thickness+safety_gap"
        )

    detector_min = (
        x["L2_Z"]
        + x["L2_Thickness"]
        + AXIAL_SAFETY_GAP
    )

    if x["Detector_Z"] <= detector_min:
        reasons.append(
            "Detector_Z<=L2_Z+L2_Thickness+safety_gap"
        )

    reasons.extend(
        validate_lens_geometry(x, "L1")
    )
    reasons.extend(
        validate_lens_geometry(x, "L2")
    )

    return len(reasons) == 0, "; ".join(reasons)


# ============================================================
# DESIGN APPLICATION
# ============================================================

def apply_design(x):
    l1_params = {
        "R1": x["L1_R1"],
        "C1": x["L1_C1"],
        "Clear1": x["L1_Clear1"],
        "Edge1": x["L1_Edge1"],
        "Thickness": x["L1_Thickness"],
        "R2": x["L1_R2"],
        "C2": x["L1_C2"],
        "Clear2": x["L1_Clear2"],
        "Edge2": x["L1_Edge2"],
    }

    l2_params = {
        "R1": x["L2_R1"],
        "C1": x["L2_C1"],
        "Clear1": x["L2_Clear1"],
        "Edge1": x["L2_Edge1"],
        "Thickness": x["L2_Thickness"],
        "R2": x["L2_R2"],
        "C2": x["L2_C2"],
        "Clear2": x["L2_Clear2"],
        "Edge2": x["L2_Edge2"],
    }

    configure_standard_lens(
        lens1,
        x["L1_Z"],
        l1_params,
        material="N-BK7",
    )

    configure_standard_lens(
        lens2,
        x["L2_Z"],
        l2_params,
        material="N-BK7",
    )

    detector.ZPosition = float(x["Detector_Z"])


# ============================================================
# LHS SAMPLING
# ============================================================

def lhs_samples(n, names, rng):
    dim = len(names)

    result = np.zeros((n, dim))

    for j, name in enumerate(names):
        permutation = rng.permutation(n)
        result[:, j] = (
            (permutation + rng.random(n))
            / float(n)
        )

    return result


def unit_to_design(unit_row, names):
    x = {}

    for value, name in zip(unit_row, names):
        low, high = VARIABLES[name]
        x[name] = float(
            low + value * (high - low)
        )

    return x


# ============================================================
# REJECTION LOGGING
# ============================================================

rejection_rows = []
global_rejection_count = 0


def log_rejection(stage, run, attempt, x, reason):
    global global_rejection_count
    global_rejection_count += 1

    row = {
        "stage": stage,
        "run": run,
        "attempt": attempt,
        "reason": reason,
    }
    row.update(x)
    rejection_rows.append(row)


# ============================================================
# EVALUATION
# ============================================================

def evaluate_design(
    x,
    stage,
    run,
    rejection_count_before_accept=0,
    last_rejection_reason="",
):
    apply_design(x)

    metrics, data = trace_and_measure()
    ranking = rank_record(metrics)

    result = {
        **x,
        **metrics,
        **ranking,
        "stage": stage,
        "run": run,
        "capture_fraction_pct": np.nan,
        "rejection_count_before_accept": (
            rejection_count_before_accept
        ),
        "last_rejection_reason": (
            last_rejection_reason
        ),
    }

    return result, data


def get_valid_lhs_design(
    stage,
    run,
    unit_row,
    names,
    rng,
):
    """
    Use the requested LHS point first. If geometry-invalid,
    record the exact reason and generate deterministic LHS
    replacements until a valid point is obtained.
    """
    x = unit_to_design(unit_row, names)
    ok, reason = valid_geometry(x)

    if ok:
        return x, 0, ""

    local_rejections = 1
    log_rejection(
        stage,
        run,
        local_rejections,
        x,
        reason,
    )

    for replacement in range(
        1,
        MAX_LHS_REPLACEMENTS + 1,
    ):
        replacement_row = lhs_samples(
            1,
            names,
            rng,
        )[0]
        x = unit_to_design(
            replacement_row,
            names,
        )
        ok, reason = valid_geometry(x)

        if ok:
            return (
                x,
                local_rejections,
                (
                    " | ".join(
                        [
                            r["reason"]
                            for r in rejection_rows[-local_rejections:]
                        ]
                    )
                ),
            )

        local_rejections += 1
        log_rejection(
            stage,
            run,
            local_rejections,
            x,
            reason,
        )

    raise RuntimeError(
        f"Could not generate a valid {stage} design "
        f"after {MAX_LHS_REPLACEMENTS} replacements."
    )


# ============================================================
# ZRD CAPTURE FRACTION
# ============================================================

def find_zrd_path(filename):
    system_dir = Path(str(system.SystemFile)).parent
    return system_dir / filename


def read_capture_fraction(
    zrd_filename,
    lens_number,
    clear_radius,
):
    """
    Count primary ZRD ray records that have a front-face segment
    on Lens 1 inside Clear1.

    ZRD provides ray-by-ray and segment-by-segment information.
    This is a direct ray-path diagnostic rather than an inference
    from detector power. The denominator is the actual number of
    primary ray records read from the ZRD file, so the result cannot
    exceed 100%.
    """
    zrd_path = find_zrd_path(zrd_filename)

    if not zrd_path.exists():
        raise RuntimeError(
            f"Expected ZRD file was not created: {zrd_path}"
        )

    reader = system.Tools.OpenRayDatabaseReader()

    if reader is None:
        raise RuntimeError(
            "Could not open ZRD reader."
        )

    try:
        reader.ZRDFile = str(zrd_path)
        reader.RunAndWaitForCompletion()

        if not reader.Succeeded:
            raise RuntimeError(
                f"ZRD reader failed: {reader.ErrorMessage}"
            )

        results = reader.GetResults()

        total_results = 0
        captured_ray_keys = set()
        all_primary_ray_keys = set()

        ok = True
        while ok:
            result = results.ReadNextResult()

            if not result:
                break

            # Python.NET normally returns:
            # success, rayNumber, waveIndex, wlUM, numSegments
            ok = bool(result[0])
            if not ok:
                break

            ray_number = int(result[1])
            wave_index = int(result[2])
            num_segments = int(result[4])
            total_results += 1
            all_primary_ray_keys.add((wave_index, ray_number))

            captured = False

            for _ in range(num_segments):
                seg = results.ReadNextSegment()

                if not seg:
                    break

                # success, level, parent, hitObj, hitFace, insideOf,
                # status, x, y, z, l, m, n, ..., intensity, path
                if not bool(seg[0]):
                    continue

                hit_obj = int(seg[3])
                hit_face = int(seg[4])
                x = float(seg[7])
                y = float(seg[8])

                if (
                    hit_obj == int(lens_number)
                    and hit_face == 1
                    and math.hypot(x, y)
                    <= abs(float(clear_radius)) + 1.0e-6
                ):
                    captured = True

            if captured:
                # Ray numbers can restart for each wavelength channel.
                # Include wavelength index so a captured ray from one
                # wavelength cannot collide with the same ray number
                # from another wavelength.
                captured_ray_keys.add((wave_index, ray_number))

        # The ZRD contains the actual primary ray records produced by
        # this trace.  Do NOT assume that Source Ellipse's
        # NumberOfAnalysisRays is the denominator: depending on the
        # ZOS-API/source configuration, the ZRD can contain a different
        # number of primary records.  Using ANALYSIS_RAYS here can
        # produce the impossible result capture > 100%.
        denominator = float(total_results)

        if denominator <= 0.0:
            return np.nan, total_results, len(captured_ray_keys)

        # A capture fraction is bounded by 0--100%.
        fraction = (
            100.0
            * len(captured_ray_keys)
            / denominator
        )
        fraction = min(100.0, max(0.0, fraction))

        return (
            fraction,
            total_results,
            len(captured_ray_keys),
        )

    finally:
        reader.Close()


def run_capture_measurement(
    lens1_z,
    clear1,
    detector_z,
    index,
):
    # Keep the rest of the two-lens baseline fixed.
    lens1.ZPosition = float(lens1_z)

    # Standard Lens requires Edge1 >= abs(Clear1).  The original
    # capture sweep changed Clear1 up to 20 mm while leaving Edge1 at
    # the 10 mm baseline, creating an invalid aperture definition.
    set_parameter(lens1, 3, clear1)
    set_parameter(lens1, 4, max(abs(float(clear1)), 10.0))

    detector.ZPosition = float(detector_z)

    zrd_name = (
        f"capture_{index:03d}.ZRD"
    )

    trace_info = run_nsc_ray_trace(
        save_zrd=True,
        zrd_name=zrd_name,
    )

    data = get_detector_data()
    metrics = compute_metrics(data)

    capture_fraction, total_results, captured = (
        read_capture_fraction(
            zrd_name,
            LENS1_NUMBER,
            clear1,
        )
    )

    try:
        find_zrd_path(zrd_name).unlink()
    except Exception:
        pass

    return {
        "lens1_z": lens1_z,
        "clear1": clear1,
        "detector_z": detector_z,
        "capture_fraction_pct": capture_fraction,
        "zrd_ray_results": total_results,
        "requested_analysis_rays": ANALYSIS_RAYS,
        "zrd_to_requested_ratio": (
            float(total_results) / float(ANALYSIS_RAYS)
            if ANALYSIS_RAYS > 0 else np.nan
        ),
        "captured_ray_count": captured,
        "wavelength_channel_count": WAVELENGTH_COUNT,
        "zrd_primary_record_audit": (
            "unique (wave_index, ray_number) keys counted during read"
        ),
        "coverage": metrics["coverage"],
        "fixed_cv": metrics["fixed_cv"],
        "active_cv": metrics["active_cv"],
        "centroid_error": metrics["centroid_error"],
        "mean": metrics["mean"],
        "peak": metrics["peak"],
        "trace_seed_supported": trace_info[
            "trace_seed_supported"
        ],
    }


# ============================================================
# BASELINE / HISTORICAL DESIGNS
# ============================================================

def run_and_save_map(title, filename):
    metrics, data = trace_and_measure()
    save_detector_map(
        data,
        OUTPUT_DIR / filename,
        title,
    )
    return metrics, data


def configure_single_baseline():
    global lens1, detector
    lens1 = get_object(LENS1_NUMBER)
    detector = get_object(DETECTOR_NUMBER)

    configure_standard_lens(
        lens1,
        single_baseline["lens_z"],
        {
            "R1": single_baseline["p1"],
            "C1": single_baseline["p2"],
            "Clear1": single_baseline["p3"],
            "Edge1": single_baseline["p4"],
            "Thickness": single_baseline["p5"],
            "R2": single_baseline["p6"],
            "C2": single_baseline["p7"],
            "Clear2": single_baseline["p8"],
            "Edge2": single_baseline["p9"],
        },
        material=single_baseline["material"],
    )

    detector.ZPosition = single_baseline["detector_z"]


# Historical 8.58% single-lens design from the previous
# parameter study.
HISTORICAL_8P58 = {
    "R1": 40.0,
    "C1": 0.0,
    "Clear1": 10.0,
    "Edge1": 10.0,
    "Thickness": 4.0,
    "R2": -60.0,
    "C2": 0.0,
    "Clear2": 10.0,
    "Edge2": 10.0,
}


# Previous flawed two-lens candidate, reconstructed from the
# recorded Phase 5E result.
OLD_FLAWED_DESIGN = {
    "L1_Z": 14.0269,
    "L1_R1": 97.6558,
    "L1_R2": 78.6242,
    "L1_C1": 0.5568,
    "L1_C2": -0.6107,
    "L1_Clear1": 13.6007,
    "L1_Edge1": 5.4380,
    "L1_Clear2": 9.8515,
    "L1_Edge2": 11.8305,
    "L1_Thickness": 6.4686,

    "L2_Z": 49.1877,
    "L2_R1": 46.0660,
    "L2_R2": -25.9081,
    "L2_C1": -0.0609,
    "L2_C2": -0.6211,
    "L2_Clear1": 9.5591,
    "L2_Edge1": 9.7570,
    "L2_Clear2": 10.7229,
    "L2_Edge2": 11.6981,
    "L2_Thickness": 4.6229,

    "Detector_Z": 90.7973,
}


# ============================================================
# CAPTURE DIAGNOSTIC
# ============================================================

def run_capture_diagnostic():
    print("\n" + "=" * 70)
    print("CAPTURE-FRACTION DIAGNOSTIC")
    print("=" * 70)

    rows = []
    index = 1

    # Start from documented two-lens baseline.
    make_two_lens_architecture()

    for z in CAPTURE_L1_Z:
        for clear1 in CAPTURE_CLEAR1:
            print(
                f"Capture test {index:02d}: "
                f"L1 Z={z:.2f} mm, "
                f"Clear1={clear1:.2f} mm"
            )

            row = run_capture_measurement(
                z,
                clear1,
                CAPTURE_DETECTOR_Z,
                index,
            )
            rows.append(row)

            print(
                f"  Capture={row['capture_fraction_pct']:.2f}% | "
                f"ZRD/requested={row['zrd_to_requested_ratio']:.3f}x | "
                f"Coverage={row['coverage']:.2f}% | "
                f"Fixed CV={row['fixed_cv']:.2f}% | "
                f"Active CV={row['active_cv']:.2f}%"
            )

            index += 1

    with open(
        CAPTURE_CSV,
        "w",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)

    save_capture_plot(rows)

    return rows


# ============================================================
# BARE-SOURCE SWEEP
# ============================================================

def set_analysis_rays_and_seed(ray_count):
    """Set the Source Ellipse analysis-ray count and fixed seed."""
    global source_data

    source_data = source.ObjectData
    source_data.NumberOfAnalysisRays = int(ray_count)

    seed_supported = False
    try:
        source_data.RandomSeed = int(RANDOM_SEED)
        seed_supported = True
    except Exception:
        pass

    return seed_supported


def run_bare_source_sweep(
    analysis_rays,
    csv_path,
    label,
    save_map=True,
):
    print("\n" + "=" * 70)
    print(f"BARE-SOURCE / NO-LENS DETECTOR-Z SWEEP — {label}")
    print("=" * 70)
    print(f"Analysis rays: {analysis_rays}")

    rows = []

    # Move both lenses safely beyond the detector range.
    # This removes lens interaction without changing the
    # uploaded model permanently.
    make_two_lens_architecture()

    original_l1_z = float(lens1.ZPosition)
    original_l2_z = float(lens2.ZPosition)
    original_detector_z = float(detector.ZPosition)
    # Do not re-read NumberOfAnalysisRays from IObject here.
    # In this OpticStudio/pythonnet build, source.ObjectData can expose
    # NumberOfAnalysisRays during initial setup but later return IObject
    # without that property. The original value was already captured at
    # startup in the baseline snapshot.
    original_analysis_rays = int(globals().get("original_analysis_rays", ANALYSIS_RAYS))

    try:
        set_analysis_rays_and_seed(analysis_rays)

        lens1.ZPosition = 500.0
        lens2.ZPosition = 520.0

        for z in BARE_SOURCE_DETECTOR_Z:
            detector.ZPosition = float(z)

            metrics, data = trace_and_measure()

            row = {
                "detector_z": z,
                "analysis_rays": analysis_rays,
                "wavelength_channel_count": WAVELENGTH_COUNT,
                "coverage": metrics["coverage"],
                "fixed_cv": metrics["fixed_cv"],
                "active_cv": metrics["active_cv"],
                "centroid_x": metrics["centroid_x"],
                "centroid_y": metrics["centroid_y"],
                "centroid_error": metrics["centroid_error"],
                "mean": metrics["mean"],
                "peak_smoothed": metrics["peak"],
            }
            rows.append(row)

            print(
                f"Z={z:6.2f} mm | "
                f"Coverage={metrics['coverage']:7.2f}% | "
                f"Fixed CV={metrics['fixed_cv']:8.2f}% | "
                f"Active CV={metrics['active_cv']:8.2f}% | "
                f"Smoothed peak={metrics['peak']:.8g}"
            )

            if (
                save_map
                and abs(z - CAPTURE_DETECTOR_Z) < 1e-9
            ):
                save_detector_map(
                    data,
                    OUTPUT_DIR / "bare_source_detector_map.png",
                    (
                        f"Bare source at detector Z={z:.1f} mm "
                        f"({analysis_rays} rays)"
                    ),
                )

    finally:
        lens1.ZPosition = original_l1_z
        lens2.ZPosition = original_l2_z
        detector.ZPosition = original_detector_z
        set_analysis_rays_and_seed(original_analysis_rays)

    with open(
        csv_path,
        "w",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)

    return rows


# ============================================================
# PARETO FRONT
# ============================================================

def pareto_front(records):
    """
    Objectives:
      maximize coverage
      minimize fixed-region CV
      minimize centroid error
    """
    feasible = [
        r for r in records
        if bool(r["feasible"])
    ]

    front = []

    for candidate in feasible:
        dominated = False

        for other in feasible:
            if other is candidate:
                continue

            no_worse = (
                other["coverage"] >= candidate["coverage"]
                and other["fixed_cv"] <= candidate["fixed_cv"]
                and other["centroid_error"]
                <= candidate["centroid_error"]
            )

            strictly_better = (
                other["coverage"] > candidate["coverage"]
                or other["fixed_cv"] < candidate["fixed_cv"]
                or other["centroid_error"]
                < candidate["centroid_error"]
            )

            if no_worse and strictly_better:
                dominated = True
                break

        if not dominated:
            front.append(candidate)

    front.sort(
        key=lambda r: (
            -r["coverage"],
            r["fixed_cv"],
            r["centroid_error"],
        )
    )

    return front


def save_pareto_plot(front):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not front:
        return

    x = [r["coverage"] for r in front]
    y = [r["fixed_cv"] for r in front]

    plt.figure(figsize=(6.5, 5.2))
    plt.scatter(x, y)
    plt.axvline(
        TARGET_COVERAGE,
        linestyle="--",
        label=f"Coverage floor={TARGET_COVERAGE:.1f}%",
    )
    plt.xlabel("Coverage [%]")
    plt.ylabel("Fixed-region CV [%]")
    plt.title("Feasible Pareto front")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "pareto_coverage_vs_fixed_cv.png",
        dpi=180,
    )
    plt.close()


# ============================================================
# VALIDATION
# ============================================================

def validate_top_candidates(records):
    feasible = [
        r for r in records
        if bool(r["feasible"])
    ]

    feasible.sort(
        key=lambda r: (
            r["fixed_cv"],
            r["centroid_error"],
        )
    )

    selected = feasible[:VALIDATION_COUNT]

    rows = []

    for idx, record in enumerate(selected, start=1):
        apply_design(record)

        metrics1, data1 = trace_and_measure()
        metrics2, data2 = trace_and_measure()

        reproducibility_delta = max(
            abs(
                metrics1[key] - metrics2[key]
            )
            for key in [
                "coverage",
                "fixed_cv",
                "active_cv",
                "centroid_error",
                "mean",
                "peak",
            ]
        )

        save_detector_map(
            data1,
            OUTPUT_DIR / f"validation_top{idx}.png",
            (
                f"Validation Top {idx}: "
                f"Coverage={metrics1['coverage']:.2f}%, "
                f"Fixed CV={metrics1['fixed_cv']:.2f}%"
            ),
        )

        row = {
            "validation_rank": idx,
            **{
                name: record[name]
                for name in VARIABLE_NAMES
            },
            "coverage": metrics1["coverage"],
            "fixed_cv": metrics1["fixed_cv"],
            "active_cv": metrics1["active_cv"],
            "centroid_x": metrics1["centroid_x"],
            "centroid_y": metrics1["centroid_y"],
            "centroid_error": metrics1["centroid_error"],
            "mean": metrics1["mean"],
            "peak": metrics1["peak"],
            "reproducibility_max_abs_delta": (
                reproducibility_delta
            ),
        }
        rows.append(row)

    if rows:
        with open(
            VALIDATION_CSV,
            "w",
            newline="",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=list(rows[0].keys()),
            )
            writer.writeheader()
            writer.writerows(rows)

    return selected, rows


# ============================================================
# MAIN
# ============================================================

def main():
    global source, detector, lens1, lens2

    print("\n" + "=" * 70)
    print("CORRECTED LED COLLIMATOR OPTIMIZATION")
    print("=" * 70)
    print(f"Baseline file: {ZOS_FILE}")
    print(f"Analysis rays: {ANALYSIS_RAYS}")
    print(f"Random seed:   {RANDOM_SEED}")
    print(
        f"Wavelength channels: {WAVELENGTH_COUNT}"
    )
    print(
        f"Source seed property supported: "
        f"{SOURCE_SEED_SUPPORTED}"
    )
    print(
        f"TARGET_COVERAGE floor: "
        f"{TARGET_COVERAGE:.2f}%"
    )
    print(
        "Ranking: feasible candidates -> fixed-region CV -> "
        "centroid error"
    )

    # --------------------------------------------------------
    # 1. Baseline single-lens map
    # --------------------------------------------------------
    print("\n[1/8] Baseline single-lens validation")
    configure_single_baseline()

    baseline_single_metrics, baseline_single_data = (
        run_and_save_map(
            "Baseline single-lens",
            "baseline_single_lens_map.png",
        )
    )

    # Historical 8.58% design map.
    configure_single_baseline()
    configure_standard_lens(
        lens1,
        single_baseline["lens_z"],
        HISTORICAL_8P58,
        material=single_baseline["material"],
    )
    historical_metrics, historical_data = (
        run_and_save_map(
            "Historical single-lens 8.58% design",
            "historical_8p58_single_lens_map.png",
        )
    )

    # --------------------------------------------------------
    # 2. Build two-lens architecture and baseline map
    # --------------------------------------------------------
    print("\n[2/8] Baseline two-lens architecture")
    make_two_lens_architecture()

    baseline_two_metrics, baseline_two_data = (
        run_and_save_map(
            "Baseline two-lens architecture",
            "baseline_two_lens_map.png",
        )
    )

    # --------------------------------------------------------
    # 3. Reconstruct old flawed candidate
    # --------------------------------------------------------
    print("\n[3/8] Reconstructing previous flawed candidate")
    apply_design(OLD_FLAWED_DESIGN)

    old_metrics, old_data = trace_and_measure()

    save_detector_map(
        old_data,
        OUTPUT_DIR / "old_flawed_4p02_map.png",
        (
            "Previous flawed objective candidate "
            f"(Coverage={old_metrics['coverage']:.2f}%, "
            f"Active CV={old_metrics['active_cv']:.2f}%)"
        ),
    )

    print(
        f"Old candidate: Coverage={old_metrics['coverage']:.3f}% | "
        f"Active pixels={old_metrics['active_pixels']} | "
        f"Active CV={old_metrics['active_cv']:.3f}% | "
        f"Fixed CV={old_metrics['fixed_cv']:.3f}%"
    )
    print(
        f"Historical 8.58%: Coverage={historical_metrics['coverage']:.3f}% | "
        f"Active pixels={historical_metrics['active_pixels']} | "
        f"Active CV={historical_metrics['active_cv']:.3f}% | "
        f"Fixed CV={historical_metrics['fixed_cv']:.3f}%"
    )

    # --------------------------------------------------------
    # 4. Capture diagnostic
    # --------------------------------------------------------
    print("\n[4/8] Capture diagnostic")
    capture_rows = run_capture_diagnostic()

    min_capture = min(
        r["capture_fraction_pct"]
        for r in capture_rows
    )
    mean_capture = float(np.mean([
        r["capture_fraction_pct"]
        for r in capture_rows
    ]))

    # Bare source comparison at the fixed search ray count.
    bare_rows = run_bare_source_sweep(
        ANALYSIS_RAYS,
        BARE_SOURCE_CSV,
        "fixed search sampling",
        save_map=True,
    )

    # Monte Carlo peak sanity check at 10x ray count.  This is not
    # used for the optimization decision; it is specifically intended
    # to test whether the smoothed peak behaves physically with Z.
    print("\n" + "=" * 70)
    print("10X BARE-SOURCE PEAK SANITY CHECK")
    print("=" * 70)
    print(
        f"Comparing {ANALYSIS_RAYS} rays against "
        f"{SANITY_ANALYSIS_RAYS} rays."
    )

    try:
        bare_sanity_rows = run_bare_source_sweep(
            SANITY_ANALYSIS_RAYS,
            BARE_SOURCE_SANITY_CSV,
            "10x Monte Carlo sanity check",
            save_map=False,
        )

        print("\n10x sanity-check peak trend:")
        for row in bare_sanity_rows:
            print(
                f"  Z={row['detector_z']:6.2f} mm | "
                f"mean={row['mean']:.8g} | "
                f"smoothed peak={row['peak_smoothed']:.8g}"
            )
    finally:
        # Restore the fixed search sampling before any architecture gate
        # or subsequent optimization work, even if the sanity sweep fails.
        set_analysis_rays_and_seed(ANALYSIS_RAYS)

    # Correlation between capture and coverage over the capture grid.
    capture_values = np.array([
        r["capture_fraction_pct"]
        for r in capture_rows
    ])
    coverage_values = np.array([
        r["coverage"]
        for r in capture_rows
    ])

    if (
        np.std(capture_values) > 0.0
        and np.std(coverage_values) > 0.0
    ):
        capture_coverage_corr = float(
            np.corrcoef(
                capture_values,
                coverage_values,
            )[0, 1]
        )
    else:
        capture_coverage_corr = np.nan

    print(
        f"Capture diagnostic: min={min_capture:.2f}%, "
        f"mean={mean_capture:.2f}%, "
        f"capture-vs-coverage correlation="
        f"{capture_coverage_corr:.3f}"
    )

    # --------------------------------------------------------
    # Architecture gate
    # --------------------------------------------------------
    if min_capture < CAPTURE_THRESHOLD:
        print("\n" + "=" * 70)
        print("PATH B — CAPTURE-LIMITED ARCHITECTURE")
        print("=" * 70)
        print(
            f"Minimum measured Lens-1 capture was "
            f"{min_capture:.2f}%, below the required "
            f"{CAPTURE_THRESHOLD:.2f}% threshold."
        )
        print(
            "The 21-variable search is intentionally NOT run."
        )
        print(
            "Interpretation: the present two-lens architecture "
            "cannot be fairly judged by shape optimization until "
            "the first-lens collection aperture is increased or "
            "the architecture is changed."
        )

        comparison_rows = [{
            "case": "baseline_single_lens",
            **baseline_single_metrics,
        }, {
            "case": "baseline_two_lens",
            **baseline_two_metrics,
        }, {
            "case": "historical_8p58_single_lens",
            **historical_metrics,
        }, {
            "case": "old_flawed_4p02_two_lens",
            **old_metrics,
        }]

        with open(
            COMPARISON_CSV,
            "w",
            newline="",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=list(comparison_rows[0].keys()),
            )
            writer.writeheader()
            writer.writerows(comparison_rows)

        write_methodology(
            status="STOPPED_AFTER_CAPTURE_DIAGNOSTIC",
            architecture_conclusion="capture-limited",
            min_capture=min_capture,
            mean_capture=mean_capture,
            capture_coverage_corr=capture_coverage_corr,
        )

        restore_original_single_lens()
        return

    # --------------------------------------------------------
    # 5. Coarse LHS search
    # --------------------------------------------------------
    print("\n[5/8] 150-point LHS coarse search")

    rng = np.random.default_rng(42)

    coarse_unit = lhs_samples(
        COARSE_RUNS,
        VARIABLE_NAMES,
        rng,
    )

    coarse_results = []

    for run in range(1, COARSE_RUNS + 1):
        x, reject_count, last_reason = (
            get_valid_lhs_design(
                "coarse",
                run,
                coarse_unit[run - 1],
                VARIABLE_NAMES,
                rng,
            )
        )

        result, _ = evaluate_design(
            x,
            "coarse",
            run,
            reject_count,
            last_reason,
        )
        coarse_results.append(result)

        print(
            f"[COARSE {run:03d}/{COARSE_RUNS}] "
            f"Coverage={result['coverage']:7.2f}% | "
            f"Fixed CV={result['fixed_cv']:8.2f}% | "
            f"Active CV={result['active_cv']:8.2f}% | "
            f"Centroid={result['centroid_error']:7.3f} | "
            f"{'FEASIBLE' if result['feasible'] else 'INFEASIBLE'}"
        )

    # --------------------------------------------------------
    # 6. Top-10 multi-start fine search
    # --------------------------------------------------------
    print("\n[6/8] Top-10 multi-start fine search")

    feasible_coarse = [
        r for r in coarse_results
        if r["feasible"]
    ]

    if feasible_coarse:
        feasible_coarse.sort(
            key=lambda r: (
                r["fixed_cv"],
                r["centroid_error"],
            )
        )
        fine_centers = feasible_coarse[
            :FINE_STARTS
        ]
    else:
        # If no coarse point reaches the floor, refine the best
        # coverage candidates only to determine whether the floor
        # was missed by coarse sampling. They remain infeasible
        # and cannot become a "best" unless they cross the floor.
        fine_centers = sorted(
            coarse_results,
            key=lambda r: (
                -r["coverage"],
                r["fixed_cv"],
            ),
        )[:FINE_STARTS]

    fine_results = []
    per_start = max(
        1,
        FINE_RUNS // max(1, len(fine_centers)),
    )

    fine_run = 0

    for start_idx, center in enumerate(
        fine_centers,
        start=1,
    ):
        local_names = VARIABLE_NAMES

        local_unit = lhs_samples(
            per_start,
            local_names,
            rng,
        )

        for local_idx in range(
            1,
            per_start + 1,
        ):
            fine_run += 1
            if fine_run > FINE_RUNS:
                break

            x = {}

            for value, name in zip(
                local_unit[local_idx - 1],
                local_names,
            ):
                low, high = VARIABLES[name]
                half = 0.15 * (high - low)

                local_low = max(
                    low,
                    center[name] - half,
                )
                local_high = min(
                    high,
                    center[name] + half,
                )

                x[name] = float(
                    local_low
                    + value
                    * (local_high - local_low)
                )

            ok, reason = valid_geometry(x)
            reject_count = 0
            last_reason = ""

            if not ok:
                reject_count += 1
                last_reason = reason
                log_rejection(
                    "fine",
                    fine_run,
                    reject_count,
                    x,
                    reason,
                )

                # Deterministic replacement.
                for _ in range(
                    MAX_LHS_REPLACEMENTS
                ):
                    replacement = lhs_samples(
                        1,
                        local_names,
                        rng,
                    )[0]

                    for value, name in zip(
                        replacement,
                        local_names,
                    ):
                        low, high = VARIABLES[name]
                        half = 0.15 * (
                            high - low
                        )
                        local_low = max(
                            low,
                            center[name] - half,
                        )
                        local_high = min(
                            high,
                            center[name] + half,
                        )
                        x[name] = float(
                            local_low
                            + value
                            * (
                                local_high
                                - local_low
                            )
                        )

                    ok, reason = valid_geometry(x)

                    if ok:
                        break

                    reject_count += 1
                    last_reason = reason
                    log_rejection(
                        "fine",
                        fine_run,
                        reject_count,
                        x,
                        reason,
                    )

                if not ok:
                    continue

            result, _ = evaluate_design(
                x,
                "fine",
                fine_run,
                reject_count,
                last_reason,
            )
            fine_results.append(result)

            print(
                f"[FINE {fine_run:03d}/{FINE_RUNS}] "
                f"Start={start_idx:02d} | "
                f"Coverage={result['coverage']:7.2f}% | "
                f"Fixed CV={result['fixed_cv']:8.2f}% | "
                f"Active CV={result['active_cv']:8.2f}% | "
                f"{'FEASIBLE' if result['feasible'] else 'INFEASIBLE'}"
            )

    all_results = coarse_results + fine_results

    # --------------------------------------------------------
    # Results CSV
    # --------------------------------------------------------
    fieldnames = (
        VARIABLE_NAMES
        + [
            "stage",
            "run",
            "coverage",
            "active_pixels",
            "active_cv",
            "fixed_cv",
            "centroid_x",
            "centroid_y",
            "centroid_error",
            "mean",
            "peak",
            "std_full",
            "capture_fraction_pct",
            "feasible",
            "objective_score",
            "ranking_rule",
            "rejection_count_before_accept",
            "last_rejection_reason",
        ]
    )

    with open(
        RESULTS_CSV,
        "w",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(all_results)

    # Geometry rejection CSV.
    if rejection_rows:
        rejection_fields = (
            [
                "stage",
                "run",
                "attempt",
                "reason",
            ]
            + VARIABLE_NAMES
        )

        with open(
            REJECTIONS_CSV,
            "w",
            newline="",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=rejection_fields,
            )
            writer.writeheader()
            writer.writerows(rejection_rows)

    # --------------------------------------------------------
    # Pareto front
    # --------------------------------------------------------
    front = pareto_front(all_results)

    if front:
        with open(
            PARETO_CSV,
            "w",
            newline="",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames,
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(front)

        save_pareto_plot(front)

    # --------------------------------------------------------
    # Corrected best candidate
    # --------------------------------------------------------
    feasible_all = [
        r for r in all_results
        if r["feasible"]
    ]

    if feasible_all:
        feasible_all.sort(
            key=lambda r: (
                r["fixed_cv"],
                r["centroid_error"],
            )
        )
        best = feasible_all[0]

        apply_design(best)
        best_metrics, best_data = trace_and_measure()

        save_detector_map(
            best_data,
            OUTPUT_DIR / "final_optimized_design_map.png",
            (
                "Corrected-objective final design: "
                f"Coverage={best_metrics['coverage']:.2f}%, "
                f"Fixed CV={best_metrics['fixed_cv']:.2f}%"
            ),
        )

        print("\n" + "=" * 70)
        print("CORRECTED SEARCH BEST FEASIBLE DESIGN")
        print("=" * 70)
        print(
            f"Coverage={best_metrics['coverage']:.4f}%"
        )
        print(
            f"Fixed-region CV={best_metrics['fixed_cv']:.4f}%"
        )
        print(
            f"Active-only CV={best_metrics['active_cv']:.4f}%"
        )
        print(
            f"Centroid=({best_metrics['centroid_x']:.6f}, "
            f"{best_metrics['centroid_y']:.6f}) mm"
        )

        # ----------------------------------------------------
        # Full-resolution validation of top 3
        # ----------------------------------------------------
        selected, validation_rows = (
            validate_top_candidates(all_results)
        )

        # ----------------------------------------------------
        # Before / after comparison
        # ----------------------------------------------------
        comparison = {
            "old_flawed_coverage": old_metrics["coverage"],
            "old_flawed_active_pixels": (
                old_metrics["active_pixels"]
            ),
            "old_flawed_active_cv": (
                old_metrics["active_cv"]
            ),
            "old_flawed_fixed_cv": (
                old_metrics["fixed_cv"]
            ),
            "corrected_best_coverage": (
                best_metrics["coverage"]
            ),
            "corrected_best_active_pixels": (
                best_metrics["active_pixels"]
            ),
            "corrected_best_active_cv": (
                best_metrics["active_cv"]
            ),
            "corrected_best_fixed_cv": (
                best_metrics["fixed_cv"]
            ),
            "corrected_best_centroid_error": (
                best_metrics["centroid_error"]
            ),
            "historical_8p58_coverage": (
                historical_metrics["coverage"]
            ),
            "historical_8p58_active_pixels": (
                historical_metrics["active_pixels"]
            ),
            "capture_min_pct": min_capture,
            "capture_mean_pct": mean_capture,
            "capture_coverage_correlation": (
                capture_coverage_corr
            ),
            "total_geometry_rejections": (
                global_rejection_count
            ),
        }

        with open(
            COMPARISON_CSV,
            "w",
            newline="",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=list(comparison.keys()),
            )
            writer.writeheader()
            writer.writerow(comparison)

        if (
            old_metrics["active_pixels"]
            < historical_metrics["active_pixels"]
        ):
            print(
                "\nShrink-to-uniform evidence: "
                "the flawed candidate illuminates fewer "
                "detector pixels than the earlier 8.58% design."
            )

        write_methodology(
            status="SEARCH_COMPLETED",
            architecture_conclusion=(
                "adequate within tested architecture"
                if best_metrics["coverage"]
                >= TARGET_COVERAGE
                else "shape-limited"
            ),
            min_capture=min_capture,
            mean_capture=mean_capture,
            capture_coverage_corr=capture_coverage_corr,
            best=best_metrics,
            pareto_count=len(front),
        )

    else:
        print("\n" + "=" * 70)
        print("NO FEASIBLE DESIGN REACHED THE COVERAGE FLOOR")
        print("=" * 70)
        print(
            f"No candidate reached "
            f"{TARGET_COVERAGE:.2f}% coverage."
        )
        print(
            "The corrected objective therefore produces no "
            "valid 'best' design. This is an architecture/parameter "
            "space finding, not a justification to relax the metric."
        )

        write_methodology(
            status="SEARCH_COMPLETED_NO_FEASIBLE_DESIGN",
            architecture_conclusion="shape-limited",
            min_capture=min_capture,
            mean_capture=mean_capture,
            capture_coverage_corr=capture_coverage_corr,
            pareto_count=len(front),
        )

    restore_original_single_lens()

    print("\nOutputs:")
    print(OUTPUT_DIR)
    print(
        "\nThe uploaded baseline .zos file was not saved."
    )


# ============================================================
# RESTORE ORIGINAL UPLOADED BASELINE IN MEMORY
# ============================================================

def restore_original_single_lens():
    global lens1, detector

    try:
        lens1 = get_object(LENS1_NUMBER)
        detector = get_object(DETECTOR_NUMBER)

        configure_standard_lens(
            lens1,
            single_baseline["lens_z"],
            {
                "R1": single_baseline["p1"],
                "C1": single_baseline["p2"],
                "Clear1": single_baseline["p3"],
                "Edge1": single_baseline["p4"],
                "Thickness": single_baseline["p5"],
                "R2": single_baseline["p6"],
                "C2": single_baseline["p7"],
                "Clear2": single_baseline["p8"],
                "Edge2": single_baseline["p9"],
            },
            material=single_baseline["material"],
        )

        detector.ZPosition = single_baseline["detector_z"]

        # If object #4 was created in memory, leave the unsaved
        # in-memory system alone; the original file remains intact.
        # A subsequent process always reloads the original file.

    except Exception:
        pass


# ============================================================
# METHODOLOGY REPORT
# ============================================================

def write_methodology(
    status,
    architecture_conclusion,
    min_capture,
    mean_capture,
    capture_coverage_corr,
    best=None,
    pareto_count=None,
):
    text = f"""
LED Collimator & Uniform Illumination System
Corrected optimization methodology

STATUS
{status}

MODEL
- Uploaded baseline file: {ZOS_FILE}
- Original baseline file is never saved by this script.
- Detector: {DETECTOR_SIZE:.1f} x {DETECTOR_SIZE:.1f} mm
- Detector sampling: {DETECTOR_PIXELS} x {DETECTOR_PIXELS}
- Coverage threshold: {COVERAGE_THRESHOLD_FRACTION:.2f} of peak
- Coverage floor: {TARGET_COVERAGE:.2f}%

WAVELENGTH AUDIT
- Configured wavelength channels: {WAVELENGTH_COUNT}
- The ZRD capture denominator is the actual number of primary ray
  records read from the ZRD, not NumberOfAnalysisRays.
- Captured ray identity includes wavelength index plus ray number so
  ray-number collisions between wavelength channels cannot undercount
  capture.

METRIC CORRECTION
- Active-only CV is retained as a diagnostic metric only.
- Corrected CV is calculated over the fixed full 100x100 detector grid.
- Peak irradiance used for the 10%-of-peak coverage threshold is the
  maximum of a 3x3 box-smoothed detector map, not the raw single-pixel
  maximum.
- Corrected ranking never uses active-only CV.
- Candidates below the coverage floor receive the same large fixed
  infeasibility penalty and are excluded from being called "best".
- Feasible candidates are ranked by fixed-region CV; centroid error
  is the secondary tiebreaker.

RAY-TRACE REPRODUCIBILITY
- Analysis rays: {ANALYSIS_RAYS}
- 10x bare-source sanity rays: {SANITY_ANALYSIS_RAYS}
- Source random seed: {RANDOM_SEED}
- Wavelength channel count: {WAVELENGTH_COUNT}
- Source seed property supported: {SOURCE_SEED_SUPPORTED}
- The NSC ray-trace random-seed API is attempted where supported.
- No claim of seed control is made when the installed API rejects it.

GEOMETRY VALIDATION
- L2 is required to be axially after L1 plus L1 thickness and safety gap.
- Detector is required to be after L2 plus L2 thickness and safety gap.
- Edge1 >= abs(Clear1) and Edge2 >= abs(Clear2) for both lenses.
- Conic sag domain is checked at clear and edge radii.
- Front/rear axial separation is sampled over radius to reject
  self-intersecting lens geometries.
- Every rejected sample is written to geometry_rejections.csv.

CAPTURE DIAGNOSTIC
- Minimum measured Lens-1 capture: {min_capture:.4f}%
- Mean measured Lens-1 capture: {mean_capture:.4f}%
- Capture-vs-coverage correlation: {capture_coverage_corr}
- ZRD capture denominator: actual primary ray records returned by the ZRD reader.
- Requested analysis-ray count is logged separately and is never used as the
  capture denominator.
- Wavelength channels audited from SystemData.Wavelengths: {WAVELENGTH_COUNT}
- Required architecture gate: > {CAPTURE_THRESHOLD:.1f}%
- Architecture conclusion: {architecture_conclusion}

BARE-SOURCE PEAK SANITY CHECK
- The normal bare-source sweep uses {ANALYSIS_RAYS} rays.
- A separate 10x sweep uses {SANITY_ANALYSIS_RAYS} rays.
- The 10x sweep is diagnostic only and is not mixed into optimization
  ranking or feasibility decisions.

SEARCH
- Coarse: {COARSE_RUNS} LHS proposals.
- Fine: {FINE_RUNS} multi-start evaluations distributed over the
  top {FINE_STARTS} coarse candidates.
- Search evaluations use full 100x100 detector data.

VALIDATION
- Top feasible candidates are re-traced at full 100x100 resolution.
- The same fixed ray count and seed are used for reproducibility.

INTERPRETATION
Simulation outputs are irradiance distributions produced by the
specified OpticStudio NSC model. Engineering interpretation is kept
separate from simulation output. In particular, a low coverage result
cannot be interpreted as an efficient uniform illuminator merely
because active-only CV is low.

LITERATURE GROUNDING
1. Zeng, Li & Ge, "Design of LED collimator for uniform illumination
   using two freeform lenses", Optica Applicata 48(3), 413-420 (2018).
   DOI: 10.5277/oa180307
2. Aslanov, Doskolovich & Moiseev, "Thin LED collimator with free-form
   lens array for illumination applications", Applied Optics 51,
   7200-7205 (2012). DOI: 10.1364/AO.51.007200
3. Zheng, Hao & Xu, "Freeform surface lens for LED uniform illumination",
   Applied Optics 48, 6627-6634 (2009). DOI: 10.1364/AO.48.006627

RESULT FILES
- optimization_results.csv
- geometry_rejections.csv
- capture_diagnostic.csv
- bare_source_sweep.csv
- bare_source_sweep_sanity_10x.csv
- pareto_front.csv
- validation_top3.csv
- before_after_comparison.csv
- detector-map PNG files
"""
    if best is not None:
        text += f"""
CORRECTED BEST FEASIBLE DESIGN
- Coverage: {best['coverage']:.6f}%
- Fixed CV: {best['fixed_cv']:.6f}%
- Active-only CV: {best['active_cv']:.6f}%
- Centroid error: {best['centroid_error']:.6f} mm
- Mean: {best['mean']:.8e}
- Peak: {best['peak']:.8e}
- Pareto-front size: {pareto_count}
"""

    METHODOLOGY_TXT.write_text(
        textwrap.dedent(text).strip() + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("\nFATAL ERROR")
        traceback.print_exc()
        restore_original_single_lens()
        raise
