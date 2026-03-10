#!/usr/bin/env python3
"""
hibeam_vertex_finder_control_full.py

Standalone HIBEAM-style vertex finder + rich control plots.

- Reads a "truth_entry_exit" style CSV (entry/exit points in TPC gas).
- Can optionally read a steering card (hibeam_vertex_new.card) to override defaults.
- Produces:
    * vertices CSV  (one row per fitted vertex candidate)
    * residuals CSV (one row per trackvertex association)
    * diagnostics PDF with many control plots matching the user's spec.

Assumptions about input CSV columns (can be overridden via CLI):
    event            : event id
    track            : track id
    gas_entry_x_cm   : entry x [cm]
    gas_entry_y_cm   : entry y [cm]
    gas_entry_z_cm   : entry z [cm]
    gas_exit_x_cm    : exit x [cm]
    gas_exit_y_cm    : exit y [cm]
    gas_exit_z_cm    : exit z [cm]
    gas_exited       : 1 if track exits gas, 0 otherwise
    vtx_x_cm         : true vertex x [cm]
    vtx_y_cm         : true vertex y [cm]
    vtx_z_cm         : true vertex z [cm]
    is_signal        : 1 for signal tracks, 0 for background (if missing, all=1)

The vertex fit:
    - Minimises sum_i ||(I - u_i u_i^T)(v - p_i)||^2 with v_z constrained to z_constraint.
    - Chi2 is built from per-track effective sigmas (sigma_x,y,z) which include:
        * intrinsic TPC resolution (tpc_sigma_x,y,z)
        * optional multiple-scattering based broadening from a beampipe.

Cuts (all steerable via CLI or steering card):
    - chi2/ndf <= chi2ndf_max       (if chi2ndf_max >= 0, else no cut)
    - max DCA to fitted vertex <= max_dca   (if max_dca >= 0, else no cut)
    - vertex radius r_vtx <= max_r_vertex   (if max_r_vertex >= 0, else no cut)

Usage example:
    python hibeam_vertex_finder_control_full.py truth_entry_exit.csv \
        --card hibeam_vertex_new.card \
        --output-vertices vertices_entry_exit.csv \
        --output-residuals vertex_residuals_entry_exit.csv \
        --output-pdf vertex_diagnostics.pdf
"""

import argparse
import itertools
import math
from typing import List, Tuple, Dict, Any, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import uproot 
import awkward as ak


DEBUG_THETA90 = False

# ----------------------------------------------------------------------
# ROOT helper
# ----------------------------------------------------------------------

def load_from_root(root_path: str, tree_name: str = "digitizedHits") -> pd.DataFrame:
    """
    Load per-track entry/exit info from ROOT TTree with vector branches.
    Returns a flat pandas DataFrame with one row per (event, track).
    Assumes positions are in cm.
    """
    t = uproot.open(root_path)[tree_name]

    branches = [
        "TPC_firstPosX", "TPC_firstPosY", "TPC_firstPosZ",
        "TPC_lastPosX",  "TPC_lastPosY",  "TPC_lastPosZ",
        "TPC_PathLength",
        "TPC_pdg",
        "PrimaryPosX", "PrimaryPosY", "PrimaryPosZ",
    ]

    # arr = t.arrays(branches, how="ak")
    arr = t.arrays(branches, library="ak")

    # track index 0..n-1 per event (used as track_id if no explicit ID exists)
    track_idx = ak.local_index(arr["TPC_firstPosX"])

    # event_id = entry number in the tree (0..Nevents-1)
    # Create an event index array with length = number of events
    nevents = len(arr["TPC_firstPosX"])
    event_id = ak.Array(np.arange(nevents, dtype=np.int64))

    # broadcast event_id to per-track shape
    event_id_bt = ak.broadcast_arrays(event_id, arr["TPC_firstPosX"])[0]

    # Truth vertex: take first primary vertex per event (PrimaryPosX/Y/Z are vectors)
    # If empty, set NaN.
    def first_or_nan(x):
        return ak.where(ak.num(x) > 0, x[:, 0], np.nan)

    truth_x = first_or_nan(arr["PrimaryPosX"])
    truth_y = first_or_nan(arr["PrimaryPosY"])
    truth_z = first_or_nan(arr["PrimaryPosZ"])

    # broadcast truth vertex to per-track shape
    truth_x_bt = ak.broadcast_arrays(truth_x, arr["TPC_firstPosX"])[0]
    truth_y_bt = ak.broadcast_arrays(truth_y, arr["TPC_firstPosX"])[0]
    truth_z_bt = ak.broadcast_arrays(truth_z, arr["TPC_firstPosX"])[0]

    # Entry/exit points
    ex = arr["TPC_firstPosX"]; ey = arr["TPC_firstPosY"]; ez = arr["TPC_firstPosZ"]
    lx = arr["TPC_lastPosX"];  ly = arr["TPC_lastPosY"];  lz = arr["TPC_lastPosZ"]

    path = arr["TPC_PathLength"]

    # Define "gas_exited_flag" robustly:
    # valid if pathlength > 0 and coordinates finite
    # finite_entry = ak.isfinite(ex) & ak.isfinite(ey) & ak.isfinite(ez)
    # finite_exit  = ak.isfinite(lx) & ak.isfinite(ly) & ak.isfinite(lz)
    finite_entry = np.isfinite(ex) & np.isfinite(ey) & np.isfinite(ez)
    finite_exit  = np.isfinite(lx) & np.isfinite(ly) & np.isfinite(lz)
    gas_exited_flag = (path > 0) & finite_entry & finite_exit

    # is_signal: default all 1 (keep old behavior)
    is_signal = ak.ones_like(ex, dtype=np.int32)

    # Flatten to 1D columns
    df = pd.DataFrame({
        "event_id": ak.to_numpy(ak.flatten(event_id_bt)),
        "track_id": ak.to_numpy(ak.flatten(track_idx)).astype(np.int64),

        "entry_x": ak.to_numpy(ak.flatten(ex)),
        "entry_y": ak.to_numpy(ak.flatten(ey)),
        "entry_z": ak.to_numpy(ak.flatten(ez)),

        "exit_x": ak.to_numpy(ak.flatten(lx)),
        "exit_y": ak.to_numpy(ak.flatten(ly)),
        "exit_z": ak.to_numpy(ak.flatten(lz)),

        "gas_exited_flag": ak.to_numpy(ak.flatten(gas_exited_flag)).astype(np.int32),

        "truth_vtx_x": ak.to_numpy(ak.flatten(truth_x_bt)),
        "truth_vtx_y": ak.to_numpy(ak.flatten(truth_y_bt)),
        "truth_vtx_z": ak.to_numpy(ak.flatten(truth_z_bt)),

        "is_signal": ak.to_numpy(ak.flatten(is_signal)).astype(np.int32),

        # optional extras (handy later)
        "TPC_PathLength": ak.to_numpy(ak.flatten(path)),
        "TPC_pdg": ak.to_numpy(ak.flatten(arr["TPC_pdg"])).astype(np.int32),
    })

    return df


def write_root_outputs(vertices_rows, residuals_rows, output_root: str) -> None:
    import uproot
    import awkward as ak
    import numpy as np

    def rows_to_columns(rows):
        if not rows:
            return {}
        keys = list(rows[0].keys())
        cols = {k: [] for k in keys}
        for r in rows:
            for k in keys:
                cols[k].append(r.get(k))
        return cols

    def normalize_columns(cols, tree_name):
        """
        Convert python lists to numpy/awkward arrays and determine branch dtypes.
        Returns (branch_types, branch_arrays).
        """
        branch_types = {}
        arrays = {}

        for k, v in cols.items():
            if len(v) == 0:
                continue

            # Jagged arrays (e.g. track_ids = list[int] per row)
            if isinstance(v[0], (list, tuple, np.ndarray)):
                arrays[k] = ak.Array(v)
                # Declare as variable-length integer array
                branch_types[k] = "var * int32"
                continue

            # Booleans -> uint8 (ROOT-friendly)
            if isinstance(v[0], (bool, np.bool_)):
                arr = np.asarray(v, dtype=np.uint8)
                arrays[k] = arr
                branch_types[k] = np.uint8
                continue

            # Integers
            if isinstance(v[0], (int, np.integer)):
                arr = np.asarray(v, dtype=np.int64)
                arrays[k] = arr
                branch_types[k] = np.int64
                continue

            # Floats (including NaN)
            if isinstance(v[0], (float, np.floating)) or v[0] is None:
                arr = np.asarray(v, dtype=np.float64)
                arrays[k] = arr
                branch_types[k] = np.float64
                continue

            # Fallback: try float, then int
            try:
                arr = np.asarray(v, dtype=np.float64)
                arrays[k] = arr
                branch_types[k] = np.float64
            except Exception:
                arr = np.asarray(v, dtype=np.int64)
                arrays[k] = arr
                branch_types[k] = np.int64

        return branch_types, arrays

    vcols = rows_to_columns(vertices_rows)
    rcols = rows_to_columns(residuals_rows)

    vtypes, varrays = normalize_columns(vcols, "vertices")
    rtypes, rarrays = normalize_columns(rcols, "residuals")

    with uproot.recreate(output_root) as f:
        if varrays:
            f.mktree("vertices", vtypes)
            f["vertices"].extend(varrays)
        if rarrays:
            f.mktree("residuals", rtypes)
            f["residuals"].extend(rarrays)

    print(f"Wrote ROOT output to: {output_root}")

# ----------------------------------------------------------------------
# Geometry / linear algebra helpers
# ----------------------------------------------------------------------

def build_track_from_entry_exit(entry: np.ndarray, exit_: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (point p, unit direction u) for a straight track from entry to exit."""
    diff = exit_ - entry
    norm = float(np.linalg.norm(diff))
    if norm <= 0.0 or (not np.isfinite(norm)):
        u = np.array([np.nan, np.nan, np.nan], dtype=float)
        p = entry.astype(float)
    else:
        u = diff / norm
        p = entry.astype(float)
    return p, u


def fit_vertex(points: List[np.ndarray], dirs: List[np.ndarray], z_constraint: float = 0.0) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Fit a vertex v = (x_v, y_v, z_constraint) by minimising sum_i ||(I - u_i u_i^T)(v - p_i)||^2.

    Returns:
        vtx: best-fit vertex [x, y, z]
        info: dict with 'A' matrix, 'cond' condition number.
    """
    if len(points) != len(dirs):
        raise ValueError("points and dirs must have same length")

    I = np.eye(3)
    A = np.zeros((3, 3), dtype=float)
    b = np.zeros(3, dtype=float)

    n_used = 0
    for p, u in zip(points, dirs):
        if not np.all(np.isfinite(u)):
            continue
        P = I - np.outer(u, u)
        A += P
        b += P @ p
        n_used += 1

    if n_used == 0 or not np.any(A):
        return np.full(3, np.nan, dtype=float), {"A": A, "cond": float("inf")}

    # Constrain z = z_constraint : solve reduced 2x2 system for x,y
    z_c = float(z_constraint)

    Axx, Axy, Axz = A[0, 0], A[0, 1], A[0, 2]
    Ayx, Ayy, Ayz = A[1, 0], A[1, 1], A[1, 2]
    bx, by, bz = b[0], b[1], b[2]

    # Effective RHS: b' = b - A * v_fixed_z where v_fixed_z = (0, 0, z_c)
    bx_eff = bx - Axz * z_c
    by_eff = by - Ayz * z_c

    A2 = np.array([[Axx, Axy],
                   [Ayx, Ayy]], dtype=float)
    b2 = np.array([bx_eff, by_eff], dtype=float)

    try:
        xy = np.linalg.solve(A2, b2)
        x, y = float(xy[0]), float(xy[1])
    except np.linalg.LinAlgError:
        x = y = float("nan")

    vtx = np.array([x, y, z_c], dtype=float)
    try:
        cond = float(np.linalg.cond(A2))
    except np.linalg.LinAlgError:
        cond = float("inf")

    return vtx, {"A": A, "cond": cond}


# ----------------------------------------------------------------------
# Multiple scattering estimate
# ----------------------------------------------------------------------

def estimate_ms_sigma_for_track(
    entry_point: np.ndarray,
    direction: np.ndarray,
    ms_density: float,
    ms_radlen_gcm2: float,
    ms_thickness_cm: float,
    ms_lever_arm_cm: float,
    ms_p_MeV: float,
    ms_min_cos: float,
) -> float:
    """
    Estimate an RMS positional uncertainty at the vertex from multiple scattering
    in a cylindrical beampipe, using the Highland formula.

    Returns sigma_ms in cm, to be added in quadrature to TPC sigmas.
    """
    if ms_density <= 0 or ms_radlen_gcm2 <= 0 or ms_thickness_cm <= 0 or ms_p_MeV <= 0:
        return 0.0

    # approximate incident angle with respect to radial normal at entry
    ex, ey = entry_point[0], entry_point[1]
    r_xy = math.hypot(ex, ey)
    if r_xy > 0:
        n_hat = np.array([ex / r_xy, ey / r_xy, 0.0], dtype=float)
        cos_inc = abs(float(np.dot(direction, n_hat)))
    else:
        cos_inc = 1.0

    cos_eff = max(cos_inc, ms_min_cos)

    # material thickness in g/cm^2
    t_gcm2 = ms_density * ms_thickness_cm / cos_eff
    # x/X0
    x_over_X0 = t_gcm2 / ms_radlen_gcm2

    if x_over_X0 <= 0.0:
        return 0.0

    # Highland formula (radians)
    theta0 = (13.6 / ms_p_MeV) * math.sqrt(x_over_X0) * (1.0 + 0.038 * math.log(x_over_X0))
    sigma_ms = theta0 * ms_lever_arm_cm  # cm

    return float(abs(sigma_ms))


# ----------------------------------------------------------------------
# Residuals / chi2
# ----------------------------------------------------------------------

def compute_residuals(vtx: np.ndarray, tracks: List[Dict[str, Any]]) -> Tuple[List[Dict[str, float]], float, int]:
    """
    For a fitted vertex, compute per-track residuals and chi2.

    Each track dict must contain:
      "p"        : point on track (3-vector)
      "u"        : unit direction (3-vector)
      "sigma_x"  : effective uncertainty in x  [cm] (for chi2)
      "sigma_y"  : effective uncertainty in y  [cm]
      "sigma_z"  : effective uncertainty in z  [cm]

    The chi2 is:
      chi2 = sum_i ((dx_i / sigma_x_i)^2 + (dy_i / sigma_y_i)^2 + (dz_i / sigma_z_i)^2)
    """
    residuals: List[Dict[str, float]] = []
    chi2 = 0.0
    N = len(tracks)
    if N == 0:
        return residuals, float("nan"), 0

    for trk in tracks:
        p = trk["p"]
        u = trk["u"]
        if not np.all(np.isfinite(u)):
            continue

        sx = float(trk.get("sigma_x", 1.0)) or 1.0
        sy = float(trk.get("sigma_y", 1.0)) or 1.0
        sz = float(trk.get("sigma_z", 1.0)) or 1.0

        t = float(np.dot(u, (vtx - p)))
        r_closest = p + t * u
        delta = r_closest - vtx

        # Full 3D DCA and transverse (xy) DCA
        dca_3d = float(np.linalg.norm(delta))
        dca_xy = float(math.hypot(delta[0], delta[1]))

        # chi2 += (delta[0] / sx) ** 2 + (delta[1] / sy) ** 2 + (delta[2] / sz) ** 2
        chi2 += (delta[0] / sx) ** 2 + (delta[1] / sy) ** 2

        residuals.append({
            "track_id": int(trk.get("track_id", -1)), #added this as it was used later but never included -Lucas 
            "delta_x": float(delta[0]),
            "delta_y": float(delta[1]),
            "delta_z": float(delta[2]),
            # By convention in this script, 'dca' is the transverse DCA in the xy plane.
            # The full 3D DCA is stored as 'dca_3d'.
            "dca": dca_xy,
            "dca_3d": dca_3d,
            "sigma_x": sx,
            "sigma_y": sy,
            "sigma_z": sz,
            "p_x": float(p[0]),
            "p_y": float(p[1]),
            "p_z": float(p[2]),
            "u_x": float(u[0]),
            "u_y": float(u[1]),
            "u_z": float(u[2]),
            "is_signal": int(trk.get("is_signal", 1)),
        })

    # ndf = 3 * N - 2  # x,y,z for each track minus 2 free vertex coords (x,y)
    ndf = 2 * N - 2
    return residuals, chi2, ndf


def add_truth_residuals(
    residuals: List[Dict[str, float]],
    tracks: List[Dict[str, Any]],
    truth_vtx: Optional[np.ndarray],
) -> None:
    """
    If a true vertex is known, compute per-track DCA and residual components to the true vertex.
    """
    if truth_vtx is None:
        return
    for res, trk in zip(residuals, tracks):
        p = trk["p"]
        u = trk["u"]
        if not np.all(np.isfinite(u)):
            res["dca_truth"] = float("nan")
            res["delta_truth_x"] = float("nan")
            res["delta_truth_y"] = float("nan")
            res["delta_truth_z"] = float("nan")
            continue
        t_true = float(np.dot(u, (truth_vtx - p)))
        r_closest_true = p + t_true * u
        delta_true = r_closest_true - truth_vtx

        # Full 3D DCA and transverse (xy) DCA to the true vertex
        dca_true_3d = float(np.linalg.norm(delta_true))
        dca_true_xy = float(math.hypot(delta_true[0], delta_true[1]))

        # Also compute the "projected" transverse miss in the foil plane (same z as the true vertex)
        # by intersecting the reconstructed track with the plane z = truth_vtx[2].
        # This is the quantity that becomes poorly constrained as � � 90�, because it involves 1/u_z.
        dca_plane_xy = float("nan")
        plane_dx = float("nan")
        plane_dy = float("nan")
        uz = float(u[2])
        if abs(uz) > 1e-6:
            t_plane = float((truth_vtx[2] - p[2]) / uz)
            r_plane = p + t_plane * u
            plane_dx = float(r_plane[0] - truth_vtx[0])
            plane_dy = float(r_plane[1] - truth_vtx[1])
            dca_plane_xy = float(math.hypot(plane_dx, plane_dy))

        # Store transverse DCA as 'dca_truth' and 3D as 'dca_truth_3d'
        # and the plane-projected transverse miss as 'dca_plane_truth'
        res["dca_truth"] = dca_true_xy
        res["dca_truth_3d"] = dca_true_3d
        res["dca_plane_truth"] = dca_plane_xy
        res["delta_truth_x"] = float(delta_true[0])
        res["delta_truth_y"] = float(delta_true[1])
        res["delta_truth_z"] = float(delta_true[2])


# ----------------------------------------------------------------------
# Steering card
# ----------------------------------------------------------------------

def load_steering_card(path: str) -> Dict[str, str]:
    """Read key = value card, stripping inline comments starting with '#'."""
    params: Dict[str, str] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "#" in line:
                line = line.split("#", 1)[0].rstrip()
            if not line:
                continue
            if "=" not in line:
                continue
            key, val = line.split("=", 1)
            params[key.strip()] = val.strip()
    return params


def apply_card_overrides(args, card_params: Dict[str, str]) -> None:
    """Override argparse defaults with values from card."""
    float_keys = {
        "tpc_sigma_x", "tpc_sigma_y", "tpc_sigma_z",
        "chi2ndf_max", "max_dca", "max_r_vertex",
        "z_constraint",
        "ms_density", "ms_radlen_gcm2", "ms_thickness_cm",
        "ms_lever_arm_cm", "ms_p_MeV", "ms_min_cos",
    }
    int_keys = {"min_comb", "max_comb", "seed", "max_events"}
    bool_keys = {"use_ms_sigma", "use_tpc_smearing"}

    for key, val in card_params.items():
        if not hasattr(args, key):
            continue

        if key in bool_keys:
            low = val.lower()
            if low in ("1", "true", "yes", "on"):
                new_value = True
            elif low in ("0", "false", "no", "off"):
                new_value = False
            else:
                print(f"[steering card] Warning: could not parse bool for '{key}', keeping {getattr(args, key)!r}")
                continue
        elif key in int_keys:
            try:
                new_value = int(val)
            except ValueError:
                print(f"[steering card] Warning: could not parse int for '{key}', keeping {getattr(args, key)!r}")
                continue
        elif key in float_keys:
            try:
                new_value = float(val)
            except ValueError:
                print(f"[steering card] Warning: could not parse float for '{key}', keeping {getattr(args, key)!r}")
                continue
        else:
            new_value = val

        print(f"[steering card] Overriding {key} = {getattr(args, key)!r} -> {new_value!r}")
        setattr(args, key, new_value)


# ----------------------------------------------------------------------
# Per-event processing
# ----------------------------------------------------------------------

def process_event(
    df_event: pd.DataFrame,
    comb_sizes: List[int],
    tpc_sigma_x: float,
    tpc_sigma_y: float,
    tpc_sigma_z: float,
    chi2ndf_max: float,
    max_dca_cut: float,
    max_r_vertex: float,
    z_constraint: float,
    use_tpc_smearing: bool,
    use_ms_sigma: bool,
    ms_density: float,
    ms_radlen_gcm2: float,
    ms_thickness_cm: float,
    ms_lever_arm_cm: float,
    ms_p_MeV: float,
    ms_min_cos: float,
) -> Dict[str, Any]:
    """
    Process a single event:
      - build tracks,
      - scan all N-track combinations (N in comb_sizes),
      - fit vertices, compute residuals and chi2,
      - decide acceptance,
      - enforce "no track in more than one accepted vertex" by choosing
        a non-overlapping subset of accepted vertices.
    """
    # Build track_meta with truth/gas-exit info
    track_meta: Dict[int, Dict[str, Any]] = {}
    for _, row in df_event.iterrows():
        tid = int(row["track_id"])
        gas_exited = int(row["gas_exited_flag"]) == 1
        is_sig = int(row.get("is_signal", 1))
        track_meta[tid] = {
            "gas_exited": gas_exited,
            "invalid_exit": (not gas_exited),
            "used_in_any_comb": False,
            "used_in_any_accepted_raw": False,   # before overlap cleanup
            "used_in_any_accepted": False,       # after overlap cleanup
            "in_combo_fail_chi2": False,
            "in_combo_fail_dca": False,
            "in_combo_fail_r": False,
            "is_signal": is_sig,
        }

    df_valid = df_event[df_event["gas_exited_flag"] == 1]

    diag_by_N: Dict[int, Dict[str, bool]] = {
        n: {
            "any_combo_considered": False,
            "any_accepted": False,
            "any_chi2_fail": False,
            "any_dca_fail": False,
            "any_r_fail": False,
        } for n in comb_sizes
    }

    if df_valid.empty:
        return {
            "vertices": [],
            "found_by_N": {n: False for n in comb_sizes},
            "track_meta": track_meta,
            "n_valid_tracks": 0,
            "diag_by_N": diag_by_N,
        }

    # Build fit tracks
    tracks: List[Dict[str, Any]] = []
    for _, row in df_valid.iterrows():
        entry = np.array([row["entry_x"], row["entry_y"], row["entry_z"]], dtype=float)
        exit_ = np.array([row["exit_x"], row["exit_y"], row["exit_z"]], dtype=float)

        if use_tpc_smearing:
            entry_sm = entry + np.random.normal(0.0, [tpc_sigma_x, tpc_sigma_y, tpc_sigma_z])
            exit_sm = exit_ + np.random.normal(0.0, [tpc_sigma_x, tpc_sigma_y, tpc_sigma_z])
        else:
            entry_sm = entry
            exit_sm = exit_

        p, u = build_track_from_entry_exit(entry_sm, exit_sm)

        sigma_tpc_x = max(tpc_sigma_x, 1e-6)
        sigma_tpc_y = max(tpc_sigma_y, 1e-6)
        sigma_tpc_z = max(tpc_sigma_z, 1e-6)

        sigma_ms = 0.0
        if use_ms_sigma:
            sigma_ms = estimate_ms_sigma_for_track(
                entry_point=entry_sm,
                direction=u,
                ms_density=ms_density,
                ms_radlen_gcm2=ms_radlen_gcm2,
                ms_thickness_cm=ms_thickness_cm,
                ms_lever_arm_cm=ms_lever_arm_cm,
                ms_p_MeV=ms_p_MeV,
                ms_min_cos=ms_min_cos,
            )

        eff_sigma_x = math.sqrt(sigma_tpc_x ** 2 + sigma_ms ** 2)
        eff_sigma_y = math.sqrt(sigma_tpc_y ** 2 + sigma_ms ** 2)
        eff_sigma_z = math.sqrt(sigma_tpc_z ** 2 + sigma_ms ** 2)

        tracks.append({
            "track_id": int(row["track_id"]),
            "p": p,
            "u": u,
            "sigma_x": eff_sigma_x,
            "sigma_y": eff_sigma_y,
            "sigma_z": eff_sigma_z,
            "is_signal": bool(row.get("is_signal", 1)),
        })

    n_valid_tracks = len(tracks)
    results_vertices: List[Dict[str, Any]] = []
    found_by_N = {n: False for n in comb_sizes}

    # True vertex (assume common for the event)
    truth_x = df_event["truth_vtx_x"].iloc[0]
    truth_y = df_event["truth_vtx_y"].iloc[0]
    truth_z = df_event["truth_vtx_z"].iloc[0]
    if np.all(np.isfinite([truth_x, truth_y, truth_z])):
        truth_vtx = np.array([truth_x, truth_y, truth_z], dtype=float)
    else:
        truth_vtx = None

    # Loop over multiplicities and combinations
    for n in comb_sizes:
        if len(tracks) < n:
            continue

        for combo in itertools.combinations(tracks, n):
            combo_tracks = list(combo)
            track_ids = [trk["track_id"] for trk in combo_tracks]
            points = [trk["p"] for trk in combo_tracks]
            dirs = [trk["u"] for trk in combo_tracks]

            vtx, info = fit_vertex(points, dirs, z_constraint=z_constraint)

            diag_by_N[n]["any_combo_considered"] = True

            if not np.all(np.isfinite(vtx)):
                continue

            r_vtx = math.hypot(vtx[0], vtx[1])

            residuals, chi2, ndf = compute_residuals(vtx, combo_tracks)
            if ndf <= 0 or (not math.isfinite(chi2)):
                continue

            add_truth_residuals(residuals, combo_tracks, truth_vtx)

            chi2ndf = chi2 / ndf
            if chi2ndf_max >= 0.0:
                chi2_ok = (chi2ndf <= chi2ndf_max)
            else:
                chi2_ok = True

            max_dca_val = max(r["dca"] for r in residuals) if residuals else float("nan")
            avg_dca_val = sum(r["dca"] for r in residuals) / len(residuals) if residuals else float("nan")

            if max_dca_cut >= 0.0:
                dca_ok = (max_dca_val <= max_dca_cut)
            else:
                dca_ok = True

            if max_r_vertex >= 0.0:
                r_ok = (r_vtx <= max_r_vertex)
            else:
                r_ok = True

            accepted = chi2_ok and dca_ok and r_ok

            if accepted:
                found_by_N[n] = True
                diag_by_N[n]["any_accepted"] = True
            else:
                if not chi2_ok:
                    diag_by_N[n]["any_chi2_fail"] = True
                if not dca_ok:
                    diag_by_N[n]["any_dca_fail"] = True
                if not r_ok:
                    diag_by_N[n]["any_r_fail"] = True

            for tid in track_ids:
                meta = track_meta.get(tid)
                if meta is None:
                    continue
                meta["used_in_any_comb"] = True
                if accepted:
                    meta["used_in_any_accepted_raw"] = True
                else:
                    if not chi2_ok:
                        meta["in_combo_fail_chi2"] = True
                    if not dca_ok:
                        meta["in_combo_fail_dca"] = True
                    if not r_ok:
                        meta["in_combo_fail_r"] = True

            # Truth-space residuals of the vertex itself
            if truth_vtx is not None:
                dv = vtx - truth_vtx
                vtx_dx_true = float(dv[0])
                vtx_dy_true = float(dv[1])
                vtx_dz_true = float(dv[2])
                vtx_d3_true = float(np.linalg.norm(dv))
            else:
                vtx_dx_true = float("nan")
                vtx_dy_true = float("nan")
                vtx_dz_true = float("nan")
                vtx_d3_true = float("nan")

            results_vertices.append({
                "n_tracks": n,
                "track_ids": track_ids,
                "vtx_x": float(vtx[0]),
                "vtx_y": float(vtx[1]),
                "vtx_z": float(vtx[2]),
                "vtx_r": float(r_vtx),
                "vtx_dx_true": vtx_dx_true,
                "vtx_dy_true": vtx_dy_true,
                "vtx_dz_true": vtx_dz_true,
                "vtx_d3_true": vtx_d3_true,
                "chi2": float(chi2),
                "ndf": int(ndf),
                "chi2ndf": float(chi2ndf),
                "cond": float(info.get("cond", float("inf"))),
                "max_dca": float(max_dca_val),
                "avg_dca": float(avg_dca_val),
                "accepted_raw": bool(accepted),
                "accepted": False,  # will be set after overlap resolution
                "chi2ndf_max_used": float(chi2ndf_max),
                "max_dca_cut_used": float(max_dca_cut),
                "max_r_vertex_cut_used": float(max_r_vertex),
                "fail_chi2": bool(not chi2_ok),
                "fail_dca": bool(not dca_ok),
                "fail_r": bool(not r_ok),
                "residuals": residuals,
            })

    # Enforce "no track in more than one accepted vertex":
    #   sort raw-accepted vertices by (n_tracks desc, chi2ndf asc)
    #   greedily accept if all its tracks are unused so far.
    used_tracks = set()
    # mark accepted flag false initially
    for v in results_vertices:
        v["accepted"] = False

    # Candidate list: only those that passed cuts
    sorted_vertices = sorted(
        [v for v in results_vertices if v["accepted_raw"]],
        key=lambda v: (-v["n_tracks"], v["chi2ndf"])
    )

    for v in sorted_vertices:
        if any(tid in used_tracks for tid in v["track_ids"]):
            continue
        v["accepted"] = True
        for tid in v["track_ids"]:
            used_tracks.add(tid)

    # Update track_meta used_in_any_accepted (exclusive)
    for meta in track_meta.values():
        meta["used_in_any_accepted"] = False
    for v in results_vertices:
        if v["accepted"]:
            for tid in v["track_ids"]:
                meta = track_meta.get(tid)
                if meta is not None:
                    meta["used_in_any_accepted"] = True

    return {
        "vertices": results_vertices,
        "found_by_N": found_by_N,
        "track_meta": track_meta,
        "n_valid_tracks": n_valid_tracks,
        "diag_by_N": diag_by_N,
    }


# ----------------------------------------------------------------------
# Control plots
# ----------------------------------------------------------------------

THETA_RANGES_DEG = [
    (40.0, 75.0),
    (75.0, 100.0),
    (100.0, 125.0),
    (125.0, 140.0),
    (87.0, 93.0),  # special around 90 deg
]


def compute_theta_deg(u_x: float, u_y: float, u_z: float) -> float:
    """Polar angle theta from track direction (degrees)."""
    u = np.array([u_x, u_y, u_z], dtype=float)
    norm = float(np.linalg.norm(u))
    if norm == 0.0 or (not np.isfinite(norm)):
        return float("nan")
    uz = u[2] / norm
    uz = max(-1.0, min(1.0, uz))
    return math.degrees(math.acos(uz))


def make_control_plots(
    vertices_csv: str,
    residuals_csv: str,
    output_pdf: str,
) -> None:
    """Create diagnostics PDF based on vertices and residuals CSVs."""
    try:
        df_vtx = pd.read_csv(vertices_csv)
    except Exception as exc:
        print(f"[control plots] Could not read vertices CSV {vertices_csv!r}: {exc}")
        return

    try:
        df_res = pd.read_csv(residuals_csv)
    except Exception as exc:
        print(f"[control plots] Could not read residuals CSV {residuals_csv!r}: {exc}")
        return

    # Basic sanity
    if "accepted" not in df_vtx.columns:
        print("[control plots] 'accepted' column missing in vertices CSV  cannot proceed.")
        return

    # Merge residuals with vertex info: to know which residual belongs to accepted vertex,
    # and access vtx_x, vtx_y, truth vertex etc.
    # Assume (event_id, vertex_index) in residuals, and we have the same in vertices.
    # If vertex_index not present, create a simple running index per event.
    if "vertex_index" not in df_vtx.columns:
        # create it: per event, 0..n_vtx-1
        df_vtx = df_vtx.copy()
        df_vtx["vertex_index"] = df_vtx.groupby("event_id").cumcount()
    if "vertex_index" not in df_res.columns:
        df_res = df_res.copy()
        df_res["vertex_index"] = df_res.groupby(["event_id", "n_tracks"]).cumcount()

    # Event-level fractions
    with PdfPages(output_pdf) as pdf:
        # 1) Fractions of events with accepted vertices by multiplicity
        if "event_id" in df_vtx.columns:
            total_events = df_vtx["event_id"].nunique()
            df_acc = df_vtx[df_vtx["accepted"] == True]

            # any accepted vertex
            events_with_any = df_acc["event_id"].nunique() if total_events > 0 else 0

            # per multiplicity
            mults = sorted(df_vtx["n_tracks"].unique())
            frac_any = events_with_any / total_events if total_events > 0 else 0.0

            fig, ax = plt.subplots(figsize=(6, 4))
            labels = ["any"]
            fracs = [frac_any]

            # 2,3,4,>4 as requested
            for N in [2, 3, 4]:
                evN = df_acc.loc[df_acc["n_tracks"] == N, "event_id"].nunique()
                fracN = evN / total_events if total_events > 0 else 0.0
                labels.append(f"N={N}")
                fracs.append(fracN)

            # N>4
            ev_gt4 = df_acc.loc[df_acc["n_tracks"] > 4, "event_id"].nunique()
            frac_gt4 = ev_gt4 / total_events if total_events > 0 else 0.0
            labels.append("N>4")
            fracs.append(frac_gt4)

            ax.bar(range(len(labels)), fracs, align="center")
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels)
            ax.set_ylim(0.0, 1.05)
            ax.set_ylabel("Fraction of events")
            ax.set_title("Fractions of events with accepted vertices")
            for i, f in enumerate(fracs):
                ax.text(i, f + 0.01, f"{f:.3f}", ha="center", va="bottom", fontsize=8)
            ax.grid(True, axis="y", alpha=0.3)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        # 2) Chi2/ndf distributions (broad and zoomed [0,4])
        if "chi2ndf" in df_vtx.columns:
            chi_all = df_vtx["chi2ndf"].replace([np.inf, -np.inf], np.nan).dropna()
            chi_acc = df_vtx.loc[df_vtx["accepted"] == True, "chi2ndf"]
            chi_acc = chi_acc.replace([np.inf, -np.inf], np.nan).dropna()

            if len(chi_all) > 0:
                fig, ax = plt.subplots(figsize=(6, 4))
                ax.hist(chi_all, bins=50, histtype="step", label="All candidates", density=True)
                if len(chi_acc) > 0:
                    ax.hist(chi_acc, bins=50, histtype="step", linestyle="--", label="Accepted", density=True)
                ax.set_xlabel("chi2 / ndf")
                ax.set_ylabel("Normalised counts")
                ax.set_title("Vertex chi2/ndf (full range)")
                ax.legend()
                ax.grid(True, alpha=0.3)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

                # Zoomed 0-4
                fig, ax = plt.subplots(figsize=(6, 4))
                ax.hist(chi_all[(chi_all >= 0.0) & (chi_all <= 4.0)], bins=40, histtype="step",
                        label="All candidates", density=True)
                ax.hist(chi_acc[(chi_acc >= 0.0) & (chi_acc <= 4.0)], bins=40, histtype="step",
                        linestyle="--", label="Accepted", density=True)
                ax.set_xlim(0.0, 4.0)
                ax.set_xlabel("chi2 / ndf")
                ax.set_ylabel("Normalised counts")
                ax.set_title("Vertex chi2/ndf (0-4)")
                ax.legend()
                ax.grid(True, alpha=0.3)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

        # 3) Vertex multiplicity distributions (all vs accepted)
        if "n_tracks" in df_vtx.columns:
            fig, ax = plt.subplots(figsize=(6, 4))
            all_mult = df_vtx["n_tracks"]
            acc_mult = df_vtx.loc[df_vtx["accepted"] == True, "n_tracks"]
            if len(all_mult) > 0:
                bins = range(int(all_mult.min()), int(all_mult.max()) + 2)
                ax.hist(all_mult, bins=bins, histtype="step", label="All candidates")
                if len(acc_mult) > 0:
                    ax.hist(acc_mult, bins=bins, histtype="step", linestyle="--", label="Accepted")
                ax.set_xlabel("Vertex multiplicity N (tracks)")
                ax.set_ylabel("Counts")
                ax.set_title("Vertex multiplicity")
                ax.legend()
                ax.grid(True, alpha=0.3)
                fig.tight_layout()
                pdf.savefig(fig)
            plt.close(fig)

        # 4) Distance to truth of vertex
        if "vtx_d3_true" in df_vtx.columns:
            fig, ax = plt.subplots(figsize=(6, 4))
            d3_all = df_vtx["vtx_d3_true"].replace([np.inf, -np.inf], np.nan).dropna()
            d3_acc = df_vtx.loc[df_vtx["accepted"] == True, "vtx_d3_true"]
            d3_acc = d3_acc.replace([np.inf, -np.inf], np.nan).dropna()
            if len(d3_all) > 0:
                ax.hist(d3_all, bins=50, histtype="step", label="All candidates", density=True)
            if len(d3_acc) > 0:
                ax.hist(d3_acc, bins=50, histtype="step", linestyle="--", label="Accepted", density=True)
            ax.set_xlabel("|v_reco - v_truth| [cm]")
            ax.set_ylabel("Normalised counts")
            ax.set_title("Vertex distance to truth")
            ax.legend()
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        # 5) Global signal / background fractions in accepted vertices
        if {"event_id", "vertex_index", "is_signal"}.issubset(df_res.columns) and {"accepted"}.issubset(df_vtx.columns):
            df_merge = df_res.merge(df_vtx[["event_id", "vertex_index", "accepted"]],
                                    on=["event_id", "vertex_index"], how="left")
            df_acc_trk = df_merge[df_merge["accepted"] == True]
            tot_sig = int((df_acc_trk["is_signal"] == 1).sum())
            tot_bg = int((df_acc_trk["is_signal"] == 0).sum())
            tot = tot_sig + tot_bg
            if tot > 0:
                frac_sig = tot_sig / tot
                frac_bg = tot_bg / tot
                fig, ax = plt.subplots(figsize=(5, 4))
                ax.bar([0, 1], [frac_sig, frac_bg])
                ax.set_xticks([0, 1])
                ax.set_xticklabels(["Signal", "Background"])
                ax.set_ylim(0.0, 1.05)
                ax.set_ylabel("Fraction of tracks in\naccepted vertices")
                ax.set_title("Signal vs background content")
                for x, y, label in zip([0, 1], [frac_sig, frac_bg], ["sig", "bg"]):
                    ax.text(x, y + 0.02, f"{y:.2f}", ha="center", va="bottom", fontsize=8)
                ax.grid(True, axis="y", alpha=0.3)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

        # For the more detailed track-parameter plots we need u_x,y,z, sigma_x,y,z, delta_x,y,z, dca_truth, etc.
        needed_cols = {"delta_x", "delta_y", "delta_z", "dca", "sigma_x", "sigma_y", "sigma_z",
                       "u_x", "u_y", "u_z", "dca_truth", "delta_truth_x", "delta_truth_y", "delta_truth_z"}
        if not needed_cols.issubset(df_res.columns):
            print("[control plots] Residuals CSV lacks some detailed columns; skipping detailed track plots.")
            print("  Missing:", needed_cols - set(df_res.columns))
            return

        # Add theta column
        df_res = df_res.copy()
        df_res["theta_deg"] = df_res.apply(
            lambda r: compute_theta_deg(r["u_x"], r["u_y"], r["u_z"]), axis=1
        )

        # Optional verbose debugging of sigma_x/sigma_y vs theta around 90 degrees.
        if DEBUG_THETA90:
            print("\n[debug theta90] === Detailed sigma_x/sigma_y vs theta diagnostics ===")
            print(f"[debug theta90] Number of residual rows: {len(df_res)}")
            print(f"[debug theta90] Residual columns: {list(df_res.columns)}")

            theta_all = df_res["theta_deg"].astype(float)
            sigx_all = df_res["sigma_x"].astype(float)
            sigy_all = df_res["sigma_y"].astype(float)

            print("[debug theta90] Global theta_deg stats:")
            print(theta_all.describe())
            print("[debug theta90] Global sigma_x stats:")
            print(sigx_all.describe())
            print("[debug theta90] Global sigma_y stats:")
            print(sigy_all.describe())

            mask_90 = (theta_all >= 87.0) & (theta_all <= 93.0)
            sub_90 = df_res[mask_90].copy()
            print(f"[debug theta90] Rows in theta slice [87,93] deg: {len(sub_90)}")

            if len(sub_90) > 0:
                print("[debug theta90] theta_deg stats in [87,93]:")
                print(sub_90["theta_deg"].describe())
                print("[debug theta90] sigma_x stats in [87,93]:")
                print(sub_90["sigma_x"].describe())
                print("[debug theta90] sigma_y stats in [87,93]:")
                print(sub_90["sigma_y"].describe())

                # Rows closest to 90 degrees
                sub_90["abs_theta_minus_90"] = (sub_90["theta_deg"] - 90.0).abs()
                sub_sorted = sub_90.sort_values("abs_theta_minus_90")
                print("\n[debug theta90] 10 rows closest to theta = 90 deg:")
                print(sub_sorted[["theta_deg", "sigma_x", "sigma_y"]].head(10).to_string(index=False))

                # Unique sigma values
                import numpy as _np
                uniq_sigx = _np.sort(sub_90["sigma_x"].dropna().unique())
                uniq_sigy = _np.sort(sub_90["sigma_y"].dropna().unique())

                print("\n[debug theta90] Unique sigma_x values in [87,93] (up to 20):")
                for v in uniq_sigx[:20]:
                    print(f"  {v:.10f}")
                if len(uniq_sigx) > 20:
                    print(f"  ... ({len(uniq_sigx) - 20} more)")

                print("\n[debug theta90] Unique sigma_y values in [87,93] (up to 20):")
                for v in uniq_sigy[:20]:
                    print(f"  {v:.10f}")
                if len(uniq_sigy) > 20:
                    print(f"  ... ({len(uniq_sigy) - 20} more)")

                # Correlations
                if len(sub_90) > 1:
                    try:
                        corr_x = _np.corrcoef(sub_90["theta_deg"], sub_90["sigma_x"])[0, 1]
                        corr_y = _np.corrcoef(sub_90["theta_deg"], sub_90["sigma_y"])[0, 1]
                        print(f"\n[debug theta90] Corr(theta, sigma_x) in [87,93] = {corr_x:.6g}")
                        print(f"[debug theta90] Corr(theta, sigma_y) in [87,93] = {corr_y:.6g}")
                    except Exception as _exc:
                        print(f"[debug theta90] Warning: failed to compute correlations: {_exc}")

                # Extremes
                print("\n[debug theta90] Largest 10 sigma_x values in [87,93]:")
                print(sub_90.sort_values("sigma_x", ascending=False)[["theta_deg", "sigma_x", "sigma_y"]]
                          .head(10).to_string(index=False))
                print("\n[debug theta90] Largest 10 sigma_y values in [87,93]:")
                print(sub_90.sort_values("sigma_y", ascending=False)[["theta_deg", "sigma_x", "sigma_y"]]
                          .head(10).to_string(index=False))
            else:
                print("[debug theta90] WARNING: No rows found in theta slice [87,93] deg.")

        # (2) DCA and angles: before vs after fit
        # Here 'dca' and 'dca_truth' are the transverse (x-y) DCA values.
        # "before" = DCA to true vertex, "after" = DCA to fitted vertex.
        # Fractional residual = (before - after)/before.
        df_res["dca_residual"] = df_res["dca_truth"] - df_res["dca"]
        df_res["dca_frac_residual"] = df_res.apply(
            lambda r: (r["dca_truth"] - r["dca"]) / r["dca_truth"] if r["dca_truth"] != 0 else np.nan,
            axis=1
        )

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        ax1, ax2 = axes
        ax1.hist(df_res["dca_truth"].dropna(), bins=60, histtype="step", label="before (to true vertex)", density=True)
        ax1.hist(df_res["dca"].dropna(), bins=60, histtype="step", linestyle="--", label="after (to fitted vertex)", density=True)
        ax1.set_xlabel("DCA [cm]")
        ax1.set_ylabel("Normalised counts")
        ax1.set_title("DCA before/after fit")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        ax2.hist(df_res["dca_residual"].dropna(), bins=60, histtype="step", density=True)
        ax2.set_xlabel("DCA_before - DCA_after [cm]")
        ax2.set_ylabel("Normalised counts")
        ax2.set_title("DCA residual")
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        
        # 6) Uncertainties sigma_x,y,z and fractional uncertainties, theta-binned
        for theta_range in THETA_RANGES_DEG:
            lo, hi = theta_range
            mask = (df_res["theta_deg"] >= lo) & (df_res["theta_deg"] < hi)
            sub = df_res[mask]
            if len(sub) == 0:
                continue
            label = f"{lo:.0f}-{hi:.0f} theta"
            fig, axes = plt.subplots(1, 3, figsize=(12, 4))
            for ax, col in zip(axes, ["sigma_x", "sigma_y", "sigma_z"]):
                ax.hist(sub[col].dropna(), bins=50, histtype="step", density=True)
                ax.set_xlabel(f"{col} [cm]")
                ax.set_ylabel("Normalised counts")
                ax.set_title(f"{col} in theta [{label}]")
                ax.grid(True, alpha=0.3)
            fig.suptitle(f"Track parameter uncertainties in theta range {label}", y=1.03)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        # 6b) DCA residual (error on DCA) in transverse plane, theta-binned
        if "dca_residual" in df_res.columns:
            for theta_range in THETA_RANGES_DEG:
                lo, hi = theta_range
                mask = (df_res["theta_deg"] >= lo) & (df_res["theta_deg"] < hi)
                sub = df_res[mask]
                if len(sub) == 0:
                    continue
                label = f"{lo:.0f}-{hi:.0f} theta"
                fig, ax = plt.subplots(figsize=(6, 4))
                ax.hist(sub["dca_residual"].dropna(), bins=60, histtype="step", density=True)
                ax.set_xlabel("DELTA DCA = DCA_truth - DCA_fit [cm]")
                ax.set_ylabel("Normalised counts")
                ax.set_title(f"DELTA DCA in theta in[{label}]")
                ax.grid(True, alpha=0.3)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

        # 6c) Transverse DCA to true vertex and plane-projected miss, theta-binned, and RMS vs theta
        if "dca_truth" in df_res.columns:
            theta_midpoints = []
            dca_truth_rms = []
            dca_plane_rms = []
            have_plane = "dca_plane_truth" in df_res.columns

            for theta_range in THETA_RANGES_DEG:
                lo, hi = theta_range
                mask = (df_res["theta_deg"] >= lo) & (df_res["theta_deg"] < hi)
                sub = df_res[mask]
                if len(sub) == 0:
                    continue
                label = f"{lo:.0f}-{hi:.0f} theta"

                vals = sub["dca_truth"].dropna().to_numpy()
                vals_plane = sub["dca_plane_truth"].dropna().to_numpy() if have_plane else np.array([])
                if vals.size == 0 and vals_plane.size == 0:
                    continue

                # Histogram(s) in this theta slice:
                #  - solid line: minimal transverse DCA to true vertex
                #  - dashed line: plane-projected transverse miss at z = z_true
                fig, ax = plt.subplots(figsize=(6, 4))
                if vals.size > 0:
                    ax.hist(vals, bins=60, histtype="step", density=True, label="minimal DCA_T(true)")
                if vals_plane.size > 0:
                    ax.hist(vals_plane, bins=60, histtype="step", density=True,
                            linestyle="--", label="plane-projected DCA_T(true)")
                ax.set_xlabel("DCA_T to true vertex [cm]")
                ax.set_ylabel("Normalised counts")
                ax.set_title(f"DCA_T(true) in theta in [{label}]")
                ax.grid(True, alpha=0.3)
                if (vals.size > 0) or (vals_plane.size > 0):
                    ax.legend()
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

                theta_midpoints.append(0.5 * (lo + hi))
                if vals.size > 0:
                    dca_truth_rms.append(float(np.sqrt(np.mean(vals ** 2))))
                else:
                    dca_truth_rms.append(float("nan"))
                if vals_plane.size > 0:
                    dca_plane_rms.append(float(np.sqrt(np.mean(vals_plane ** 2))))
                else:
                    dca_plane_rms.append(float("nan"))

            # Summary plot: RMS transverse DCA to true vertex vs theta
            if len(theta_midpoints) > 0:
                fig, ax = plt.subplots(figsize=(6, 4))
                ax.plot(theta_midpoints, dca_truth_rms, marker="o", label="minimal DCA_T(true)")
                if any(np.isfinite(dca_plane_rms)):
                    ax.plot(theta_midpoints, dca_plane_rms, marker="s", label="plane-projected DCA_T(true)")
                ax.set_xlabel("Theta [deg]")
                ax.set_ylabel("RMS(DCA_T to true vertex) [cm]")
                ax.set_title("RMS transverse DCA to true vertex vs �")
                ax.grid(True, alpha=0.3)
                ax.legend()
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

# 7) Track residuals x_track - x_vertex and x_track - x_true
        # For straight lines, r_closest to fitted vertex is vtx + delta; r_closest_true is truth_vtx + delta_truth.
        if {"vtx_x", "vtx_y", "vtx_z"}.issubset(df_vtx.columns):
            df_merge = df_res.merge(df_vtx[["event_id", "vertex_index", "vtx_x", "vtx_y", "vtx_z", "accepted"]],
                                    on=["event_id", "vertex_index"], how="left")
            df_acc = df_merge[df_merge["accepted"] == True].copy()

            # Reconstruct r_closest and r_closest_true
            df_acc["x_track_fit"] = df_acc["vtx_x"] + df_acc["delta_x"]
            df_acc["y_track_fit"] = df_acc["vtx_y"] + df_acc["delta_y"]
            df_acc["z_track_fit"] = df_acc["vtx_z"] + df_acc["delta_z"]

            df_acc["x_track_true"] = df_acc["vtx_x"] + df_acc["delta_truth_x"]
            df_acc["y_track_true"] = df_acc["vtx_y"] + df_acc["delta_truth_y"]
            df_acc["z_track_true"] = df_acc["vtx_z"] + df_acc["delta_truth_z"]

            # Residuals (x_track_fit - x_vertex) is just delta_x; similarly vs true is delta_truth_x
            for theta_range in THETA_RANGES_DEG + [("all", "all")]:
                if theta_range == ("all", "all"):
                    sub = df_acc
                    theta_label = "all theta"
                else:
                    lo, hi = theta_range
                    mask = (df_acc["theta_deg"] >= lo) & (df_acc["theta_deg"] < hi)
                    sub = df_acc[mask]
                    theta_label = f"{lo:.0f}-{hi:.0f}�"
                if len(sub) == 0:
                    continue
                fig, axes = plt.subplots(1, 3, figsize=(12, 4))
                for ax, comp, label_axis in zip(
                    axes,
                    ["x", "y", "z"],
                    ["Delta x [cm]", "Delta y [cm]", "Delta z [cm]"],
                ):
                    dx_fit = sub[f"delta_{comp}"]
                    dx_true = sub[f"delta_truth_{comp}"]
                    ax.hist(dx_fit.dropna(), bins=60, histtype="step", label="track - fitted vtx", density=True)
                    ax.hist(dx_true.dropna(), bins=60, histtype="step", linestyle="--",
                            label="track - true vtx", density=True)
                    ax.set_xlabel(label_axis)
                    ax.set_ylabel("Normalised counts")
                    ax.set_title(f"{comp}-residuals, {theta_label}")
                    ax.grid(True, alpha=0.3)
                    ax.legend()
                fig.suptitle(f"Track residuals to fitted and true vertex ({theta_label})", y=1.03)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

    print(f"[control plots] Wrote diagnostics PDF to: {output_pdf}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="HIBEAM vertex finder using entry/exit points with steering card and rich control plots."
    )
    # parser.add_argument("input_csv", help="Input CSV with entry/exit track parameters.")
    parser.add_argument("input_root", help="Input ROOT file containing the digitizedHits TTree.")
    parser.add_argument("--tree", dest="tree_name", default="digitizedHits",help="TTree name (default: digitizedHits).")
    parser.add_argument("--output-root", dest="output_root", default="vertex_products.root", help="Output ROOT file with TTrees (vertices, residuals).")

    parser.add_argument("--output-vertices", dest="output_vertices", default="vertices_entry_exit.csv")
    parser.add_argument("--output-residuals", dest="output_residuals", default="vertex_residuals_entry_exit.csv")
    parser.add_argument("--output-pdf", dest="output_pdf", default="vertex_diagnostics.pdf")

    parser.add_argument("--event-col", default="event", help="Name of event ID column.")
    parser.add_argument("--track-col", default="track", help="Name of track ID column.")

    parser.add_argument("--entry-x-col", dest="entry_x_col", default="gas_entry_x_cm")
    parser.add_argument("--entry-y-col", dest="entry_y_col", default="gas_entry_y_cm")
    parser.add_argument("--entry-z-col", dest="entry_z_col", default="gas_entry_z_cm")

    parser.add_argument("--exit-x-col", dest="exit_x_col", default="gas_exit_x_cm")
    parser.add_argument("--exit-y-col", dest="exit_y_col", default="gas_exit_y_cm")
    parser.add_argument("--exit-z-col", dest="exit_z_col", default="gas_exit_z_cm")

    parser.add_argument("--gas-exited-col", dest="gas_exited_col", default="gas_exited")

    parser.add_argument("--truth-x-col", dest="truth_x_col", default="vtx_x_cm")
    parser.add_argument("--truth-y-col", dest="truth_y_col", default="vtx_y_cm")
    parser.add_argument("--truth-z-col", dest="truth_z_col", default="vtx_z_cm")

    # TPC resolution
    # parser.add_argument("--tpc-sigma-x", dest="tpc_sigma_x", type=float, default=0.05)
    # parser.add_argument("--tpc-sigma-y", dest="tpc_sigma_y", type=float, default=0.05)
    # parser.add_argument("--tpc-sigma-z", dest="tpc_sigma_z", type=float, default=0.10)

    # parser.add_argument("--tpc-sigma-x", dest="tpc_sigma_x", type=float, default=0.78)
    # parser.add_argument("--tpc-sigma-y", dest="tpc_sigma_y", type=float, default=1.06)
    # parser.add_argument("--tpc-sigma-z", dest="tpc_sigma_z", type=float, default=0.10)

    parser.add_argument("--tpc-sigma-x", dest="tpc_sigma_x", type=float, default=0.42)
    parser.add_argument("--tpc-sigma-y", dest="tpc_sigma_y", type=float, default=0.58)
    parser.add_argument("--tpc-sigma-z", dest="tpc_sigma_z", type=float, default=0.10)

    # Quality cuts
    # parser.add_argument("--chi2ndf-max", dest="chi2ndf_max", type=float, default=5.0,
    #                     help="Max chi2/ndf. Set <0 to disable chi2 cut.")
    parser.add_argument("--chi2ndf-max", dest="chi2ndf_max", type=float, default=15.0,
                        help="Max chi2/ndf. Set <0 to disable chi2 cut.")
    # parser.add_argument("--max-dca", dest="max_dca", type=float, default=-1.0,
    #                     help="Max DCA to fitted vertex [cm]. Set <0 to disable.")
    parser.add_argument("--max-dca", dest="max_dca", type=float, default=3,
                        help="Max DCA to fitted vertex [cm]. Set <0 to disable.")
    parser.add_argument("--max-r-vertex", dest="max_r_vertex", type=float, default=20,
                        help="Max vertex radius in foil plane [cm]. Set <0 to disable.")

    parser.add_argument("--min-comb", dest="min_comb", type=int, default=2)
    parser.add_argument("--max-comb", dest="max_comb", type=int, default=4)

    parser.add_argument("--z-constraint", dest="z_constraint", type=float, default=0.0)

    parser.add_argument("--use-tpc-smearing", dest="use_tpc_smearing", action="store_true")
    # parser.set_defaults(use_tpc_smearing=True)
    parser.set_defaults(use_tpc_smearing=False)

    parser.add_argument("--use-ms-sigma", dest="use_ms_sigma", action="store_true")
    parser.set_defaults(use_ms_sigma=False)

    # MS parameters
    parser.add_argument("--ms-density", dest="ms_density", type=float, default=2.7)
    parser.add_argument("--ms-radlen-gcm2", dest="ms_radlen_gcm2", type=float, default=24.0)
    parser.add_argument("--ms-thickness-cm", dest="ms_thickness_cm", type=float, default=0.5)
    parser.add_argument("--ms-lever-arm-cm", dest="ms_lever_arm_cm", type=float, default=20.0)
    parser.add_argument("--ms-p-MeV", dest="ms_p_MeV", type=float, default=200.0)
    parser.add_argument("--ms-min-cos", dest="ms_min_cos", type=float, default=0.05)

    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--card", type=str, default=None,
                        help="Optional steering card to override defaults.")
    parser.add_argument("--max-events", dest="max_events", type=int, default=0,
                        help="If >0, stop after processing this many events.")

    parser.add_argument("--debug-theta90", dest="debug_theta90", action="store_true",
                        help="Print detailed debug info for sigma_x/sigma_y near thetaH90� in control plots.")

    args = parser.parse_args()

    global DEBUG_THETA90
    DEBUG_THETA90 = bool(getattr(args, "debug_theta90", False))
    if DEBUG_THETA90:
        print("[debug theta90] Enabled: will print sigma_x/sigma_y vs theta diagnostics in make_control_plots().")

    # Apply steering card overrides
    card_params: Dict[str, str] = {}
    if args.card is not None:
        card_params = load_steering_card(args.card)
        apply_card_overrides(args, card_params)

    # Steering card can also override output filenames
    if "output_vertices" in card_params:
        args.output_vertices = card_params["output_vertices"]
    if "output_residuals" in card_params:
        args.output_residuals = card_params["output_residuals"]

    # # Read input CSV
    # df = pd.read_csv(args.input_csv)

    # rename_map = {
    #     args.event_col: "event_id",
    #     args.track_col: "track_id",
    #     args.entry_x_col: "entry_x",
    #     args.entry_y_col: "entry_y",
    #     args.entry_z_col: "entry_z",
    #     args.exit_x_col: "exit_x",
    #     args.exit_y_col: "exit_y",
    #     args.exit_z_col: "exit_z",
    #     args.gas_exited_col: "gas_exited_flag",
    #     args.truth_x_col: "truth_vtx_x",
    #     args.truth_y_col: "truth_vtx_y",
    #     args.truth_z_col: "truth_vtx_z",
    # }
    # df = df.rename(columns=rename_map)

    # Read input ROOT (digitizedHits)
    df = load_from_root(args.input_root, tree_name=args.tree_name)

    # Only tracks that exited gas (keep same behavior as before)
    df = df[df["gas_exited_flag"] == 1].copy()

    required_cols = [
        "event_id", "track_id",
        "entry_x", "entry_y", "entry_z",
        "exit_x", "exit_y", "exit_z",
        "gas_exited_flag",
        "truth_vtx_x", "truth_vtx_y", "truth_vtx_z",
        "is_signal",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns after ROOT load: {missing}")

    if "is_signal" not in df.columns:
        df["is_signal"] = 1

    required_cols = [
        "event_id", "track_id",
        "entry_x", "entry_y", "entry_z",
        "exit_x", "exit_y", "exit_z",
        "gas_exited_flag",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns in input: {missing}")

    # Only tracks that exited gas
    df = df[df["gas_exited_flag"] == 1]

    if args.seed is not None:
        np.random.seed(args.seed)

    comb_sizes = list(range(args.min_comb, args.max_comb + 1))

    print("[config] comb_sizes =", comb_sizes)
    print("[config] chi2ndf_max =", args.chi2ndf_max,
          "max_dca =", args.max_dca,
          "max_r_vertex =", args.max_r_vertex)
    print(f"[chi2 model] tpc_sigma = ({args.tpc_sigma_x}, {args.tpc_sigma_y}, {args.tpc_sigma_z}) [cm]")
    print(f"[chi2 model] use_tpc_smearing = {args.use_tpc_smearing}, use_ms_sigma = {args.use_ms_sigma}")

    if args.use_ms_sigma:
        print(
            "[chi2 model] MS params: density={:.3f} g/cm^3, X0={:.1f} g/cm^2, "
            "thickness={:.3f} cm, lever_arm={:.1f} cm, p={:.1f} MeV/c, min_cos={:.3f}".format(
                args.ms_density, args.ms_radlen_gcm2,
                args.ms_thickness_cm, args.ms_lever_arm_cm,
                args.ms_p_MeV, args.ms_min_cos
            )
        )

    # Prepare writers
    vtx_header_written = False
    res_header_written = False
    vtx_writer = None
    res_writer = None

    total_events = 0
    found_counts = {n: 0 for n in comb_sizes}

    # Event-level diagnostics per N
    event_diag: Dict[int, Dict[str, int]] = {
        n: {
            "too_few_tracks": 0,
            "geN_no_accepted": 0,
            "fail_chi2": 0,
            "fail_dca": 0,
            "fail_r": 0,
            "other_fail": 0,
        } for n in comb_sizes
    }

    # Track-level rejection indices
    total_tracks = 0
    idx_p_no_exit = 0
    idx_q_no_comb = 0
    idx_r_chi2 = 0
    idx_s_dca = 0
    idx_t_r = 0

    all_vertex_rows = []
    all_residual_rows = []

    with open(args.output_vertices, "w", newline="") as f_vtx, open(args.output_residuals, "w", newline="") as f_res:
        import csv
        # Group by event
        for event_id, df_event in df.groupby("event_id"):
            total_events += 1
            if args.max_events > 0 and total_events > args.max_events:
                break

            ev_result = process_event(
                df_event=df_event,
                comb_sizes=comb_sizes,
                tpc_sigma_x=args.tpc_sigma_x,
                tpc_sigma_y=args.tpc_sigma_y,
                tpc_sigma_z=args.tpc_sigma_z,
                chi2ndf_max=args.chi2ndf_max,
                max_dca_cut=args.max_dca,
                max_r_vertex=args.max_r_vertex,
                z_constraint=args.z_constraint,
                use_tpc_smearing=args.use_tpc_smearing,
                use_ms_sigma=args.use_ms_sigma,
                ms_density=args.ms_density,
                ms_radlen_gcm2=args.ms_radlen_gcm2,
                ms_thickness_cm=args.ms_thickness_cm,
                ms_lever_arm_cm=args.ms_lever_arm_cm,
                ms_p_MeV=args.ms_p_MeV,
                ms_min_cos=args.ms_min_cos,
            )

            n_valid_tracks = ev_result.get("n_valid_tracks", 0)
            diag_by_N = ev_result.get("diag_by_N", {
                n: {
                    "any_combo_considered": False,
                    "any_accepted": False,
                    "any_chi2_fail": False,
                    "any_dca_fail": False,
                    "any_r_fail": False,
                } for n in comb_sizes
            })

            for n in comb_sizes:
                if n_valid_tracks < n:
                    event_diag[n]["too_few_tracks"] += 1
                else:
                    if not ev_result["found_by_N"].get(n, False):
                        event_diag[n]["geN_no_accepted"] += 1
                        dN = diag_by_N.get(n, {})
                        if dN.get("any_chi2_fail", False):
                            event_diag[n]["fail_chi2"] += 1
                        if dN.get("any_dca_fail", False):
                            event_diag[n]["fail_dca"] += 1
                        if dN.get("any_r_fail", False):
                            event_diag[n]["fail_r"] += 1
                        if dN.get("any_combo_considered", False) and not (
                            dN.get("any_chi2_fail", False)
                            or dN.get("any_dca_fail", False)
                            or dN.get("any_r_fail", False)
                        ):
                            event_diag[n]["other_fail"] += 1

            for n in comb_sizes:
                if ev_result["found_by_N"].get(n, False):
                    found_counts[n] += 1

            # Track-level indices
            for _, meta in ev_result["track_meta"].items():
                total_tracks += 1
                used_acc = meta["used_in_any_accepted"]
                used_comb = meta["used_in_any_comb"]
                invalid = meta["invalid_exit"]

                if invalid and not used_acc:
                    idx_p_no_exit += 1
                if (not invalid) and (not used_comb) and (not used_acc):
                    idx_q_no_comb += 1

                if used_comb and (not used_acc):
                    if meta["in_combo_fail_chi2"]:
                        idx_r_chi2 += 1
                    if meta["in_combo_fail_dca"]:
                        idx_s_dca += 1
                    if meta["in_combo_fail_r"]:
                        idx_t_r += 1

            # Write out vertices and residuals
            vertices = ev_result["vertices"]
            # assign vertex_index per event
            for vtx_index, v in enumerate(vertices):
                vtx_row = {
                    "event_id": event_id,
                    "vertex_index": vtx_index,
                    "n_tracks": v["n_tracks"],
                    "vtx_x": v["vtx_x"],
                    "vtx_y": v["vtx_y"],
                    "vtx_z": v["vtx_z"],
                    "vtx_r": v["vtx_r"],
                    "vtx_dx_true": v["vtx_dx_true"],
                    "vtx_dy_true": v["vtx_dy_true"],
                    "vtx_dz_true": v["vtx_dz_true"],
                    "vtx_d3_true": v["vtx_d3_true"],
                    "chi2": v["chi2"],
                    "ndf": v["ndf"],
                    "chi2ndf": v["chi2ndf"],
                    "cond": v["cond"],
                    "max_dca": v["max_dca"],
                    "avg_dca": v["avg_dca"],
                    "accepted_raw": v["accepted_raw"],
                    "accepted": v["accepted"],
                    "chi2ndf_max_used": v["chi2ndf_max_used"],
                    "max_dca_cut_used": v["max_dca_cut_used"],
                    "max_r_vertex_cut_used": v["max_r_vertex_cut_used"],
                    "fail_chi2": v["fail_chi2"],
                    "fail_dca": v["fail_dca"],
                    "fail_r": v["fail_r"],
                }
                if not vtx_header_written:
                    import csv
                    vtx_writer = csv.DictWriter(f_vtx, fieldnames=list(vtx_row.keys()))
                    vtx_writer.writeheader()
                    vtx_header_written = True
                vtx_writer.writerow(vtx_row)
                all_vertex_rows.append(vtx_row)

                for res in v["residuals"]:
                    res_row = {
                        "event_id": event_id,
                        "vertex_index": vtx_index,
                        "n_tracks": v["n_tracks"],
                        "track_id": res.get("track_id", None),
                        "is_signal": res.get("is_signal", 1),
                        "delta_x": res["delta_x"],
                        "delta_y": res["delta_y"],
                        "delta_z": res["delta_z"],
                        "dca": res["dca"],
                        "dca_truth": res.get("dca_truth", float("nan")),
                        "dca_truth_3d": res.get("dca_truth_3d", float("nan")),
                        "dca_plane_truth": res.get("dca_plane_truth", float("nan")),
                        "delta_truth_x": res.get("delta_truth_x", float("nan")),
                        "delta_truth_y": res.get("delta_truth_y", float("nan")),
                        "delta_truth_z": res.get("delta_truth_z", float("nan")),
                        "sigma_x": res.get("sigma_x", float("nan")),
                        "sigma_y": res.get("sigma_y", float("nan")),
                        "sigma_z": res.get("sigma_z", float("nan")),
                        "p_x": res.get("p_x", float("nan")),
                        "p_y": res.get("p_y", float("nan")),
                        "p_z": res.get("p_z", float("nan")),
                        "u_x": res.get("u_x", float("nan")),
                        "u_y": res.get("u_y", float("nan")),
                        "u_z": res.get("u_z", float("nan")),
                    }
                    if not res_header_written:
                        import csv
                        res_writer = csv.DictWriter(f_res, fieldnames=list(res_row.keys()))
                        res_writer.writeheader()
                        res_header_written = True
                    res_writer.writerow(res_row)
                    all_residual_rows.append(res_row)

    print(f"Wrote vertex-level information to: {args.output_vertices}")
    print(f"Wrote per-track residuals to: {args.output_residuals}")

    print("\n=== Event-level vertex finding summary ===")
    print(f"Total events (with >=1 exiting track): {total_events}")
    if total_events > 0:
        for n in comb_sizes:
            frac = found_counts[n] / total_events
            print(f"Events with at least one accepted {n}-track vertex: "
                  f"{found_counts[n]} / {total_events}  (fraction = {frac:.3f})")

    print("\n=== Track-level rejection indices (p, q, r, s, t) ===")
    print(f"Total tracks seen (after gas_exited_flag == 1): {total_tracks}")
    print("p: invalid exit (gas_exit ~ (0,0,0); track never considered in vertex combos)")
    print(f"   count p = {idx_p_no_exit}")
    print("q: valid exit but never appears in any N>=min_comb combination")
    print(f"   count q = {idx_q_no_comb}")
    print("r: appears only in combinations whose vertices fail the chi2 cut")
    print(f"   count r = {idx_r_chi2}")
    print("s: appears only in combinations whose vertices fail the max_dca cut")
    print(f"   count s = {idx_s_dca}")
    print("t: appears only in combinations whose vertices fail the max_r_vertex cut")
    print(f"   count t = {idx_t_r}")
    print("(Note: a single track can contribute to multiple indices r/s/t.)")

    # Detailed N-wise diagnostics
    for n in comb_sizes:
        diag = event_diag[n]
        print(f"\nN = {n}")
        print(f"  Events with <{n} valid tracks: {diag['too_few_tracks']}")
        print(f"  Events with >={n} tracks but no accepted {n}-track vertex: {diag['geN_no_accepted']}")
        if diag["geN_no_accepted"] > 0:
            print("    Among those:")
            print(f"      with at least one {n}-combo failing chi2 cut     : {diag['fail_chi2']}")
            print(f"      with at least one {n}-combo failing max_dca cut : {diag['fail_dca']}")
            print(f"      with at least one {n}-combo failing r_vertex cut: {diag['fail_r']}")
            print(f"      'other' (combos considered but no explicit cut flag): {diag['other_fail']}")

    try:
        write_root_outputs(all_vertex_rows, all_residual_rows, args.output_root)
    except Exception as exc:
        print(f"[write_root_outputs] Warning: failed to write to root: {exc}")


    # Finally, control plots
    try:
        make_control_plots(
            vertices_csv=args.output_vertices,
            residuals_csv=args.output_residuals,
            output_pdf=args.output_pdf,
        )
    except Exception as exc:
        print(f"[control plots] Warning: failed to produce PDF: {exc}")


if __name__ == "__main__":
    main()

