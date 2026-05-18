#!/usr/bin/env python3
"""
hibeam_vertex_finder_streaming.py

HIBEAM-style vertex finder + rich control plots.
Streams ROOT input event-by-event and writes output incrementally —
no large in-memory accumulation.

Key changes vs. original:
  - load_from_root() replaced by iter_events_from_root() which yields
    one event's DataFrame at a time via uproot.iterate (step=1).
  - CSV output removed; ROOT output written with incremental uproot extend().
  - PDF control plots read directly from the (now fully-written) output ROOT
    trees rather than from CSV files.
  - all_vertex_rows / all_residual_rows lists eliminated.
"""

import argparse
import itertools
import math
from typing import List, Tuple, Dict, Any, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import uproot
import awkward as ak


DEBUG_THETA90 = False

# -----------------------------------------------------------------------
# Streaming ROOT reader — yields one event's DataFrame at a time
# -----------------------------------------------------------------------

_ROOT_BRANCHES = [
    "TPC_firstPosX", "TPC_firstPosY", "TPC_firstPosZ",
    "TPC_lastPosX",  "TPC_lastPosY",  "TPC_lastPosZ",
    # "TPC_PathLength",
    "TPC_pdg",
    "PrimaryPosX", "PrimaryPosY", "PrimaryPosZ",
]


def iter_events_from_root(root_path: str, tree_name: str = "digitizedHits"):
    """
    Yield one pandas DataFrame per event (TTree entry).
    Uses uproot.iterate with step=1 so only one entry lives in memory at a time.
    """
    with uproot.open(root_path) as f:
        tree = f[tree_name]
        for batch_idx, arrays in enumerate(
            tree.iterate(_ROOT_BRANCHES, library="ak", step_size=1)
        ):
            # 'arrays' contains exactly one TTree entry (one event).
            event_id = batch_idx  # entry number == event id

            ex = arrays["TPC_firstPosX"][0]
            ey = arrays["TPC_firstPosY"][0]
            ez = arrays["TPC_firstPosZ"][0]
            lx = arrays["TPC_lastPosX"][0]
            ly = arrays["TPC_lastPosY"][0]
            lz = arrays["TPC_lastPosZ"][0]
            # path = arrays["TPC_PathLength"][0]
            path = np.sqrt((lx - ex)**2 + (ly - ey)**2 + (lz - ez)**2)
            pdg  = arrays["TPC_pdg"][0]

            pvx = arrays["PrimaryPosX"][0]
            pvy = arrays["PrimaryPosY"][0]
            pvz = arrays["PrimaryPosZ"][0]

            n_tracks = len(ex)
            if n_tracks == 0:
                continue

            # Truth vertex — first primary vertex or NaN
            truth_x = float(pvx[0]) if len(pvx) > 0 else float("nan")
            truth_y = float(pvy[0]) if len(pvy) > 0 else float("nan")
            truth_z = float(pvz[0]) if len(pvz) > 0 else float("nan")

            ex_np = ak.to_numpy(ex).astype(float)
            ey_np = ak.to_numpy(ey).astype(float)
            ez_np = ak.to_numpy(ez).astype(float)
            lx_np = ak.to_numpy(lx).astype(float)
            ly_np = ak.to_numpy(ly).astype(float)
            lz_np = ak.to_numpy(lz).astype(float)
            path_np = ak.to_numpy(path).astype(float)

            finite_entry = np.isfinite(ex_np) & np.isfinite(ey_np) & np.isfinite(ez_np)
            finite_exit  = np.isfinite(lx_np) & np.isfinite(ly_np) & np.isfinite(lz_np)
            gas_exited   = ((path_np > 0) & finite_entry & finite_exit).astype(np.int32)

            df = pd.DataFrame({
                "event_id":       np.full(n_tracks, event_id, dtype=np.int64),
                "track_id":       np.arange(n_tracks, dtype=np.int64),
                "entry_x":        ex_np,
                "entry_y":        ey_np,
                "entry_z":        ez_np,
                "exit_x":         lx_np,
                "exit_y":         ly_np,
                "exit_z":         lz_np,
                "gas_exited_flag": gas_exited,
                "truth_vtx_x":    np.full(n_tracks, truth_x),
                "truth_vtx_y":    np.full(n_tracks, truth_y),
                "truth_vtx_z":    np.full(n_tracks, truth_z),
                "is_signal":      np.ones(n_tracks, dtype=np.int32),
                "TPC_PathLength": path_np,
                "TPC_pdg":        ak.to_numpy(pdg).astype(np.int32),
            })

            yield event_id, df


# -----------------------------------------------------------------------
# Incremental ROOT writer
# -----------------------------------------------------------------------

class RootWriter:
    """
    Writes vertex and residual rows incrementally to a ROOT file.
    Call .flush(vtx_rows, res_rows) after each event, then .close().

    We use uproot.recreate + mktree on first flush, then extend() on all
    subsequent flushes so memory is bounded by one event's worth of rows.
    """

    # Branch type maps (Python/numpy type -> uproot dtype string)
    _VTX_DTYPES: Dict[str, Any] = {
        "event_id": np.int64, "vertex_index": np.int64, "n_tracks": np.int64,
        "vtx_x": np.float64, "vtx_y": np.float64, "vtx_z": np.float64,
        "vtx_r": np.float64,
        "vtx_dx_true": np.float64, "vtx_dy_true": np.float64,
        "vtx_dz_true": np.float64, "vtx_d3_true": np.float64,
        "chi2": np.float64, "ndf": np.int64, "chi2ndf": np.float64,
        "cond": np.float64, "max_dca": np.float64, "avg_dca": np.float64,
        "accepted_raw": np.int32, "accepted": np.int32,
        "chi2ndf_max_used": np.float64, "max_dca_cut_used": np.float64,
        "max_r_vertex_cut_used": np.float64,
        "fail_chi2": np.int32, "fail_dca": np.int32, "fail_r": np.int32,
    }

    _RES_DTYPES: Dict[str, Any] = {
        "event_id": np.int64, "vertex_index": np.int64, "n_tracks": np.int64,
        "track_id": np.int64, "is_signal": np.int32,
        "delta_x": np.float64, "delta_y": np.float64, "delta_z": np.float64,
        "dca": np.float64,
        "dca_truth": np.float64, "dca_truth_3d": np.float64,
        "dca_plane_truth": np.float64,
        "delta_truth_x": np.float64, "delta_truth_y": np.float64,
        "delta_truth_z": np.float64,
        "sigma_x": np.float64, "sigma_y": np.float64, "sigma_z": np.float64,
        "p_x": np.float64, "p_y": np.float64, "p_z": np.float64,
        "u_x": np.float64, "u_y": np.float64, "u_z": np.float64,
    }

    def __init__(self, path: str):
        self._path = path
        self._file = uproot.recreate(path)
        self._vtx_created = False
        self._res_created = False

    def _rows_to_arrays(
        self, rows: List[Dict], dtype_map: Dict[str, Any]
    ) -> Dict[str, np.ndarray]:
        """Convert list-of-dicts to column arrays with the declared dtypes."""
        out = {}
        for col, dtype in dtype_map.items():
            vals = [r.get(col, 0 if np.issubdtype(dtype, np.integer) else float("nan"))
                    for r in rows]
            out[col] = np.asarray(vals, dtype=dtype)
        return out

    def flush(
        self,
        vtx_rows: List[Dict],
        res_rows: List[Dict],
    ) -> None:
        if vtx_rows:
            arrays = self._rows_to_arrays(vtx_rows, self._VTX_DTYPES)
            if not self._vtx_created:
                self._file.mktree("vertices", {k: v.dtype for k, v in arrays.items()})
                self._vtx_created = True
            self._file["vertices"].extend(arrays)

        if res_rows:
            arrays = self._rows_to_arrays(res_rows, self._RES_DTYPES)
            if not self._res_created:
                self._file.mktree("residuals", {k: v.dtype for k, v in arrays.items()})
                self._res_created = True
            self._file["residuals"].extend(arrays)

    def close(self) -> None:
        self._file.close()
        print(f"[RootWriter] Closed {self._path}")


# -----------------------------------------------------------------------
# Geometry helpers
# -----------------------------------------------------------------------

def build_track_from_entry_exit(
    entry: np.ndarray, exit_: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    diff = exit_ - entry
    norm = float(np.linalg.norm(diff))
    if norm <= 0.0 or not np.isfinite(norm):
        return entry.astype(float), np.full(3, np.nan)
    return entry.astype(float), diff / norm


# def fit_vertex(
#     points: List[np.ndarray],
#     dirs: List[np.ndarray],
#     z_constraint: float = 0.0,
# ) -> Tuple[np.ndarray, Dict[str, Any]]:
#     """
#     Fit vertex v = (x, y, z_constraint) minimising
#         sum_i || (I - u_i u_i^T)(v - p_i) ||^2
#     with z fixed to z_constraint.
#     """
#     I3 = np.eye(3)
#     A = np.zeros((3, 3))
#     b = np.zeros(3)

#     for p, u in zip(points, dirs):
#         if not np.all(np.isfinite(u)):
#             continue
#         P = I3 - np.outer(u, u)
#         A += P
#         b += P @ p

#     z_c = float(z_constraint)
#     A2 = A[:2, :2]
#     b2 = b[:2] - A[:2, 2] * z_c

#     try:
#         xy = np.linalg.solve(A2, b2)
#         x, y = float(xy[0]), float(xy[1])
#     except np.linalg.LinAlgError:
#         x = y = float("nan")

#     vtx = np.array([x, y, z_c])
#     try:
#         cond = float(np.linalg.cond(A2))
#     except np.linalg.LinAlgError:
#         cond = float("inf")

#     return vtx, {"A": A, "cond": cond}

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

# def fit_vertex(points, dirs, z_constraint=0.0):
#     I3 = np.eye(3)
#     A = np.zeros((3, 3))
#     b = np.zeros(3)

#     for p, u in zip(points, dirs):
#         if not np.all(np.isfinite(u)):
#             continue

#         u = u / np.linalg.norm(u)

#         P = I3 - np.outer(u, u)

#         Cp = np.diag([0.05**2, 0.05**2, 0.10**2])

#         Wi = P @ np.linalg.pinv(Cp) @ P

#         A += Wi
#         b += Wi @ p

#     z = float(z_constraint)

#     A2 = A[:2, :2]
#     b2 = b[:2] - A[:2, 2] * z

#     try:
#         xy = np.linalg.solve(A2, b2)
#         vtx = np.array([xy[0], xy[1], z])
#     except np.linalg.LinAlgError:
#         vtx = np.array([np.nan, np.nan, z])

#     cond = np.linalg.cond(A2) if np.all(np.isfinite(A2)) else np.inf

#     return vtx, {"A": A, "cond": cond}

# -----------------------------------------------------------------------
# Residuals / chi2
# -----------------------------------------------------------------------

# def compute_residuals(
#     vtx: np.ndarray,
#     tracks: List[Dict],
# ) -> Tuple[List[Dict], float, int]:
#     residuals = []
#     chi2 = 0.0
#     for trk in tracks:
#         p, u = trk["p"], trk["u"]
#         if not np.all(np.isfinite(u)):
#             continue
#         sx = float(trk.get("sigma_x", 1.0)) or 1.0
#         sy = float(trk.get("sigma_y", 1.0)) or 1.0
#         sz = float(trk.get("sigma_z", 1.0)) or 1.0
#         t = float(np.dot(u, vtx - p))
#         delta = (p + t * u) - vtx
#         dca_xy = float(math.hypot(delta[0], delta[1]))
#         dca_3d = float(np.linalg.norm(delta))
#         chi2 += (delta[0]/sx)**2 + (delta[1]/sy)**2 # + (delta[2]/sz)**2
#         # chi2 += (d1/sigma_perp)**2 + (d2/sigma_perp)**2
#         residuals.append({
#             "track_id": int(trk.get("track_id", -1)),
#             "delta_x": float(delta[0]), "delta_y": float(delta[1]), "delta_z": float(delta[2]),
#             "dca": dca_xy, "dca_3d": dca_3d,
#             "sigma_x": sx, "sigma_y": sy, "sigma_z": sz,
#             "p_x": float(p[0]), "p_y": float(p[1]), "p_z": float(p[2]),
#             "u_x": float(u[0]), "u_y": float(u[1]), "u_z": float(u[2]),
#             "is_signal": int(trk.get("is_signal", 1)),
#         })
#     # ndf = 3 * len(tracks) - 2
#     ndf = 2 * len(tracks) - 2
#     # ndf = 2 * (len(tracks) - 2)

#     return residuals, chi2, ndf

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

        # Full 3D DCA and transverse (x-y) DCA
        dca_3d = float(np.linalg.norm(delta))
        dca_xy = float(math.hypot(delta[0], delta[1]))

        # chi2 += (delta[0] / sx) ** 2 + (delta[1] / sy) ** 2 + (delta[2] / sz) ** 2
        chi2 += (delta[0] / sx) ** 2 + (delta[1] / sy) ** 2

        residuals.append({
            "track_id": int(trk.get("track_id", -1)), #added this as it was used later but never included -Lucas 
            "delta_x": float(delta[0]),
            "delta_y": float(delta[1]),
            "delta_z": float(delta[2]),
            # By convention in this script, 'dca' is the transverse DCA in the x-y plane.
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

    ndf = 2 * N - 2
    return residuals, chi2, ndf


def add_truth_residuals(
    residuals: List[Dict],
    tracks: List[Dict],
    truth_vtx: Optional[np.ndarray],
) -> None:
    if truth_vtx is None:
        return
    for res, trk in zip(residuals, tracks):
        p, u = trk["p"], trk["u"]
        if not np.all(np.isfinite(u)):
            for k in ("dca_truth", "dca_truth_3d", "dca_plane_truth",
                      "delta_truth_x", "delta_truth_y", "delta_truth_z"):
                res[k] = float("nan")
            continue
        t_true = float(np.dot(u, truth_vtx - p))
        delta_true = (p + t_true * u) - truth_vtx
        res["dca_truth"]    = float(math.hypot(delta_true[0], delta_true[1]))
        res["dca_truth_3d"] = float(np.linalg.norm(delta_true))
        res["delta_truth_x"] = float(delta_true[0])
        res["delta_truth_y"] = float(delta_true[1])
        res["delta_truth_z"] = float(delta_true[2])
        uz = float(u[2])
        if abs(uz) > 1e-6:
            t_pl = (truth_vtx[2] - p[2]) / uz
            r_pl = p + t_pl * u
            res["dca_plane_truth"] = float(math.hypot(r_pl[0]-truth_vtx[0], r_pl[1]-truth_vtx[1]))
        else:
            res["dca_plane_truth"] = float("nan")


# -----------------------------------------------------------------------
# Per-event processing  (unchanged logic, no accumulation)
# -----------------------------------------------------------------------

def process_event(
    df_event: pd.DataFrame,
    comb_sizes: List[int],
    tpc_sigma_x: float, tpc_sigma_y: float, tpc_sigma_z: float,
    chi2ndf_max: float, max_dca_cut: float, max_r_vertex: float,
    z_constraint: float,
    use_tpc_smearing: bool
) -> Dict[str, Any]:

    track_meta: Dict[int, Dict] = {}
    for _, row in df_event.iterrows():
        tid = int(row["track_id"])
        track_meta[tid] = {
            "gas_exited": int(row["gas_exited_flag"]) == 1,
            "invalid_exit": int(row["gas_exited_flag"]) != 1,
            "used_in_any_comb": False,
            "used_in_any_accepted_raw": False,
            "used_in_any_accepted": False,
            "in_combo_fail_chi2": False,
            "in_combo_fail_dca": False,
            "in_combo_fail_r": False,
            "is_signal": int(row.get("is_signal", 1)),
        }

    df_valid = df_event[df_event["gas_exited_flag"] == 1]
    diag_by_N = {n: dict(any_combo_considered=False, any_accepted=False,
                         any_chi2_fail=False, any_dca_fail=False, any_r_fail=False)
                 for n in comb_sizes}

    if df_valid.empty:
        return {
            "vertices": [], "found_by_N": {n: False for n in comb_sizes},
            "track_meta": track_meta, "n_valid_tracks": 0, "diag_by_N": diag_by_N,
        }

    # Build fit tracks
    tracks: List[Dict] = []
    for _, row in df_valid.iterrows():
        pdg = int(row["TPC_pdg"])
        if abs(pdg) == 11:
            continue
        entry = np.array([row["entry_x"], row["entry_y"], row["entry_z"]], dtype=float)
        exit_ = np.array([row["exit_x"],  row["exit_y"],  row["exit_z"]],  dtype=float)
        if use_tpc_smearing:
            entry = entry + np.random.normal(0.0, [tpc_sigma_x, tpc_sigma_y, tpc_sigma_z])
            exit_ = exit_  + np.random.normal(0.0, [tpc_sigma_x, tpc_sigma_y, tpc_sigma_z])
        p, u = build_track_from_entry_exit(entry, exit_)
        sx = max(tpc_sigma_x, 1e-6); sy = max(tpc_sigma_y, 1e-6); sz = max(tpc_sigma_z, 1e-6)
        L = float(np.linalg.norm(exit_ - entry)) # track length 
        tracks.append({"track_id": int(row["track_id"]), "p": p, "u": u,
                        "sigma_x": sx, "sigma_y": sy, "sigma_z": sz,
                        "track_length": L,
                        "is_signal": bool(row.get("is_signal", 1))})

    truth_x, truth_y, truth_z = (df_event["truth_vtx_x"].iloc[0],
                                  df_event["truth_vtx_y"].iloc[0],
                                  df_event["truth_vtx_z"].iloc[0])
    truth_vtx = (np.array([truth_x, truth_y, truth_z])
                 if np.all(np.isfinite([truth_x, truth_y, truth_z])) else None)

    results_vertices: List[Dict] = []
    found_by_N = {n: False for n in comb_sizes}

    for n in comb_sizes:
        if len(tracks) < n:
            continue
        for combo in itertools.combinations(tracks, n):
            combo_tracks = list(combo)
            track_ids = [t["track_id"] for t in combo_tracks]
            vtx, info = fit_vertex([t["p"] for t in combo_tracks],
                                   [t["u"] for t in combo_tracks],
                                   z_constraint=z_constraint)
            diag_by_N[n]["any_combo_considered"] = True
            if not np.all(np.isfinite(vtx)):
                continue
            r_vtx = math.hypot(vtx[0], vtx[1])
            residuals, chi2, ndf = compute_residuals(vtx, combo_tracks)
            if ndf <= 0 or not math.isfinite(chi2):
                continue
            add_truth_residuals(residuals, combo_tracks, truth_vtx)
            chi2ndf = chi2 / ndf
            max_dca_val = max(r["dca"] for r in residuals) if residuals else float("nan")
            avg_dca_val = sum(r["dca"] for r in residuals) / len(residuals) if residuals else float("nan")
            chi2_ok = (chi2ndf <= chi2ndf_max) if chi2ndf_max >= 0 else True
            dca_ok  = (max_dca_val <= max_dca_cut) if max_dca_cut >= 0 else True
            r_ok    = (r_vtx <= max_r_vertex) if max_r_vertex >= 0 else True
            accepted = chi2_ok and dca_ok and r_ok

            if accepted:
                found_by_N[n] = True; diag_by_N[n]["any_accepted"] = True
            else:
                if not chi2_ok: diag_by_N[n]["any_chi2_fail"] = True
                if not dca_ok:  diag_by_N[n]["any_dca_fail"]  = True
                if not r_ok:    diag_by_N[n]["any_r_fail"]    = True

            for tid in track_ids:
                m = track_meta.get(tid)
                if m is None: continue
                m["used_in_any_comb"] = True
                if accepted:
                    m["used_in_any_accepted_raw"] = True
                else:
                    if not chi2_ok: m["in_combo_fail_chi2"] = True
                    if not dca_ok:  m["in_combo_fail_dca"]  = True
                    if not r_ok:    m["in_combo_fail_r"]    = True

            dv = vtx - truth_vtx if truth_vtx is not None else np.full(3, np.nan)
            results_vertices.append({
                "n_tracks": n, "track_ids": track_ids,
                "vtx_x": float(vtx[0]), "vtx_y": float(vtx[1]), "vtx_z": float(vtx[2]),
                "vtx_r": float(r_vtx),
                "vtx_dx_true": float(dv[0]), "vtx_dy_true": float(dv[1]),
                "vtx_dz_true": float(dv[2]), "vtx_d3_true": float(np.linalg.norm(dv)),
                "chi2": float(chi2), "ndf": int(ndf), "chi2ndf": float(chi2ndf),
                "cond": float(info.get("cond", float("inf"))),
                "max_dca": float(max_dca_val), "avg_dca": float(avg_dca_val),
                "accepted_raw": bool(accepted), "accepted": False,
                "chi2ndf_max_used": float(chi2ndf_max),
                "max_dca_cut_used": float(max_dca_cut),
                "max_r_vertex_cut_used": float(max_r_vertex),
                "fail_chi2": bool(not chi2_ok), "fail_dca": bool(not dca_ok),
                "fail_r": bool(not r_ok),
                "residuals": residuals,
            })

    # Overlap resolution: greedy, best chi2ndf first, no shared tracks
    used_tracks: set = set()
    for v in results_vertices:
        v["accepted"] = False
    for v in sorted([v for v in results_vertices if v["accepted_raw"]],
                    key=lambda v: (-v["n_tracks"], v["chi2ndf"])):
        if any(tid in used_tracks for tid in v["track_ids"]):
            continue
        v["accepted"] = True
        used_tracks.update(v["track_ids"])

    for m in track_meta.values():
        m["used_in_any_accepted"] = False
    for v in results_vertices:
        if v["accepted"]:
            for tid in v["track_ids"]:
                if tid in track_meta:
                    track_meta[tid]["used_in_any_accepted"] = True

    return {
        "vertices": results_vertices, "found_by_N": found_by_N,
        "track_meta": track_meta, "n_valid_tracks": len(tracks),
        "diag_by_N": diag_by_N,
    }


# -----------------------------------------------------------------------
# Control plots — read from the output ROOT file
# -----------------------------------------------------------------------

THETA_RANGES_DEG = [
    (40.0, 75.0), (75.0, 100.0), (100.0, 125.0), (125.0, 140.0), (87.0, 93.0),
]


def _compute_theta_deg(ux, uy, uz) -> np.ndarray:
    """Vectorised theta calculation from direction arrays."""
    norm = np.sqrt(ux**2 + uy**2 + uz**2)
    safe = norm > 0
    cos_theta = np.where(safe, uz / np.where(safe, norm, 1.0), 0.0)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))


def make_control_plots(output_root: str, output_pdf: str) -> None:
    """
    Read vertices and residuals trees from the output ROOT file and
    produce the diagnostics PDF.  Uses chunked reads so the plots
    can also be generated without loading the full file into memory.
    """
    print(f"[control plots] Reading from {output_root} ...")

    # ---- load vertices tree (usually much smaller than residuals) ----
    try:
        with uproot.open(output_root) as f:
            vtx_tree = f["vertices"]
            df_vtx = vtx_tree.arrays(library="pd")
    except Exception as exc:
        print(f"[control plots] Cannot read 'vertices' tree: {exc}")
        return

    # ---- load residuals tree in chunks and concatenate ----
    # For very large files you could build the plots incrementally,
    # but for typical HIBEAM sizes a full read is fine.  If you hit
    # memory limits here, reduce step_size below.
    try:
        with uproot.open(output_root) as f:
            res_chunks = []
            for chunk in f["residuals"].iterate(library="pd", step_size=50_000):
                res_chunks.append(chunk)
            df_res = pd.concat(res_chunks, ignore_index=True) if res_chunks else pd.DataFrame()
    except Exception as exc:
        print(f"[control plots] Cannot read 'residuals' tree: {exc}")
        return

    if df_vtx.empty:
        print("[control plots] vertices tree is empty — nothing to plot.")
        return

    # Ensure vertex_index column exists
    if "vertex_index" not in df_vtx.columns:
        df_vtx["vertex_index"] = df_vtx.groupby("event_id").cumcount()
    if not df_res.empty and "vertex_index" not in df_res.columns:
        df_res["vertex_index"] = df_res.groupby(["event_id", "n_tracks"]).cumcount()

    with PdfPages(output_pdf) as pdf:

        # 1) Event-level acceptance fractions
        total_events = df_vtx["event_id"].nunique()
        df_acc = df_vtx[df_vtx["accepted"] == 1]

        labels = ["any"]
        fracs  = [df_acc["event_id"].nunique() / total_events if total_events else 0.0]
        for N in [2, 3, 4]:
            ev_N = df_acc.loc[df_acc["n_tracks"] == N, "event_id"].nunique()
            fracs.append(ev_N / total_events if total_events else 0.0)
            labels.append(f"N={N}")
        ev_gt4 = df_acc.loc[df_acc["n_tracks"] > 4, "event_id"].nunique()
        fracs.append(ev_gt4 / total_events if total_events else 0.0)
        labels.append("N>4")

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(range(len(labels)), fracs)
        ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels)
        ax.set_ylim(0, 1.05); ax.set_ylabel("Fraction of events")
        ax.set_title("Fractions of events with accepted vertices")
        for i, f in enumerate(fracs):
            ax.text(i, f + 0.01, f"{f:.3f}", ha="center", fontsize=8)
        ax.grid(True, axis="y", alpha=0.3); fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # 2) chi2/ndf distributions
        if "chi2ndf" in df_vtx.columns:
            chi_all = df_vtx["chi2ndf"].replace([np.inf, -np.inf], np.nan).dropna()
            chi_acc = df_acc["chi2ndf"].replace([np.inf, -np.inf], np.nan).dropna()
            if len(chi_all) > 0:
                for xlim, title_sfx in [(None, "full range"), ((0, 4), "0–4")]:
                    fig, ax = plt.subplots(figsize=(6, 4))
                    data_all = chi_all[(chi_all>=xlim[0])&(chi_all<=xlim[1])] if xlim else chi_all
                    data_acc = chi_acc[(chi_acc>=xlim[0])&(chi_acc<=xlim[1])] if xlim else chi_acc
                    ax.hist(data_all, bins=50, histtype="step", label="All", density=True)
                    if len(data_acc):
                        ax.hist(data_acc, bins=50, histtype="step", ls="--", label="Accepted", density=True)
                    if xlim: ax.set_xlim(*xlim)
                    ax.set_xlabel("chi2/ndf"); ax.set_ylabel("Norm. counts")
                    ax.set_title(f"Vertex chi2/ndf ({title_sfx})")
                    ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
                    pdf.savefig(fig); plt.close(fig)

        # 3) Vertex multiplicity
        if "n_tracks" in df_vtx.columns:
            bins = range(int(df_vtx["n_tracks"].min()), int(df_vtx["n_tracks"].max()) + 2)
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.hist(df_vtx["n_tracks"], bins=bins, histtype="step", label="All")
            if len(df_acc):
                ax.hist(df_acc["n_tracks"], bins=bins, histtype="step", ls="--", label="Accepted")
            ax.set_xlabel("N tracks"); ax.set_ylabel("Counts")
            ax.set_title("Vertex multiplicity")
            ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
            pdf.savefig(fig); plt.close(fig)

        # 4) Distance to truth
        if "vtx_d3_true" in df_vtx.columns:
            d3_all = df_vtx["vtx_d3_true"].replace([np.inf,-np.inf],np.nan).dropna()
            d3_acc = df_acc["vtx_d3_true"].replace([np.inf,-np.inf],np.nan).dropna()
            fig, ax = plt.subplots(figsize=(6, 4))
            if len(d3_all): ax.hist(d3_all, bins=50, histtype="step", label="All", density=True)
            if len(d3_acc): ax.hist(d3_acc, bins=50, histtype="step", ls="--", label="Accepted", density=True)
            ax.set_xlabel("|v_reco - v_truth| [cm]"); ax.set_ylabel("Norm. counts")
            ax.set_title("Vertex distance to truth")
            ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
            pdf.savefig(fig); plt.close(fig)

        if df_res.empty:
            print("[control plots] residuals tree is empty — skipping track-level plots.")
            pdf.savefig(plt.figure()); plt.close()
            print(f"[control plots] Wrote {output_pdf}")
            return

        # 5) Signal/background content in accepted vertices
        if "is_signal" in df_res.columns:
            df_res_merge = df_res.merge(
                df_vtx[["event_id", "vertex_index", "accepted"]],
                on=["event_id", "vertex_index"], how="left"
            )
            df_acc_trk = df_res_merge[df_res_merge["accepted"] == 1]
            tot_sig = int((df_acc_trk["is_signal"] == 1).sum())
            tot_bg  = int((df_acc_trk["is_signal"] == 0).sum())
            tot = tot_sig + tot_bg
            if tot > 0:
                fig, ax = plt.subplots(figsize=(5, 4))
                ax.bar([0, 1], [tot_sig/tot, tot_bg/tot])
                ax.set_xticks([0,1]); ax.set_xticklabels(["Signal","Background"])
                ax.set_ylim(0, 1.05); ax.set_ylabel("Fraction of tracks in accepted vtx")
                ax.set_title("Signal vs background content")
                for x, y in zip([0,1],[tot_sig/tot, tot_bg/tot]):
                    ax.text(x, y+0.02, f"{y:.2f}", ha="center", fontsize=8)
                ax.grid(True, axis="y", alpha=0.3); fig.tight_layout()
                pdf.savefig(fig); plt.close(fig)

        needed = {"delta_x","delta_y","delta_z","dca","sigma_x","sigma_y","sigma_z",
                  "u_x","u_y","u_z","dca_truth","delta_truth_x","delta_truth_y","delta_truth_z"}
        if not needed.issubset(df_res.columns):
            print("[control plots] Missing columns for detailed track plots:",
                  needed - set(df_res.columns))
            return

        # Compute theta once
        df_res = df_res.copy()
        df_res["theta_deg"] = _compute_theta_deg(
            df_res["u_x"].to_numpy(), df_res["u_y"].to_numpy(), df_res["u_z"].to_numpy()
        )
        df_res["dca_residual"] = df_res["dca_truth"] - df_res["dca"]

        # 6) DCA before/after fit
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        ax1.hist(df_res["dca_truth"].dropna(), bins=60, histtype="step", label="before (true)", density=True)
        ax1.hist(df_res["dca"].dropna(),        bins=60, histtype="step", ls="--", label="after (fit)", density=True)
        ax1.set_xlabel("DCA [cm]"); ax1.set_ylabel("Norm. counts")
        ax1.set_title("DCA before/after fit"); ax1.legend(); ax1.grid(True, alpha=0.3)
        ax2.hist(df_res["dca_residual"].dropna(), bins=60, histtype="step", density=True)
        ax2.set_xlabel("DCA_before - DCA_after [cm]"); ax2.set_ylabel("Norm. counts")
        ax2.set_title("DCA residual"); ax2.grid(True, alpha=0.3)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # 7) Sigma x,y,z per theta bin
        for lo, hi in THETA_RANGES_DEG:
            sub = df_res[(df_res["theta_deg"] >= lo) & (df_res["theta_deg"] < hi)]
            if sub.empty: continue
            fig, axes = plt.subplots(1, 3, figsize=(12, 4))
            for ax, col in zip(axes, ["sigma_x","sigma_y","sigma_z"]):
                ax.hist(sub[col].dropna(), bins=50, histtype="step", density=True)
                ax.set_xlabel(f"{col} [cm]"); ax.set_ylabel("Norm. counts")
                ax.set_title(f"{col} theta [{lo:.0f}–{hi:.0f}]")
                ax.grid(True, alpha=0.3)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # 8) Delta DCA per theta bin
        for lo, hi in THETA_RANGES_DEG:
            sub = df_res[(df_res["theta_deg"] >= lo) & (df_res["theta_deg"] < hi)]
            if sub.empty: continue
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.hist(sub["dca_residual"].dropna(), bins=60, histtype="step", density=True)
            ax.set_xlabel("DELTA DCA [cm]"); ax.set_ylabel("Norm. counts")
            ax.set_title(f"DELTA DCA theta [{lo:.0f}–{hi:.0f}]")
            ax.grid(True, alpha=0.3); fig.tight_layout()
            pdf.savefig(fig); plt.close(fig)

        # 9) DCA to true vertex + plane-projected, and RMS vs theta
        have_plane = "dca_plane_truth" in df_res.columns
        theta_mids, rms_dca, rms_plane = [], [], []
        for lo, hi in THETA_RANGES_DEG:
            sub = df_res[(df_res["theta_deg"] >= lo) & (df_res["theta_deg"] < hi)]
            if sub.empty: continue
            vals       = sub["dca_truth"].dropna().to_numpy()
            vals_plane = sub["dca_plane_truth"].dropna().to_numpy() if have_plane else np.array([])
            if vals.size == 0 and vals_plane.size == 0: continue
            fig, ax = plt.subplots(figsize=(6, 4))
            if vals.size:
                ax.hist(vals,       bins=60, histtype="step", density=True, label="minimal DCA_T(true)")
            if vals_plane.size:
                ax.hist(vals_plane, bins=60, histtype="step", density=True, ls="--",
                        label="plane-projected DCA_T(true)")
            ax.set_xlabel("DCA_T to true vtx [cm]"); ax.set_ylabel("Norm. counts")
            ax.set_title(f"DCA_T(true) theta [{lo:.0f}–{hi:.0f}]")
            ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
            pdf.savefig(fig); plt.close(fig)
            theta_mids.append(0.5*(lo+hi))
            rms_dca.append(float(np.sqrt(np.mean(vals**2)))    if vals.size   else float("nan"))
            rms_plane.append(float(np.sqrt(np.mean(vals_plane**2))) if vals_plane.size else float("nan"))

        if theta_mids:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(theta_mids, rms_dca,   marker="o", label="minimal DCA_T(true)")
            if any(np.isfinite(rms_plane)):
                ax.plot(theta_mids, rms_plane, marker="s", label="plane-projected DCA_T(true)")
            ax.set_xlabel("Theta [deg]"); ax.set_ylabel("RMS(DCA_T to true vtx) [cm]")
            ax.set_title("RMS transverse DCA to true vtx vs theta")
            ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
            pdf.savefig(fig); plt.close(fig)

        # 10) Track residuals to fitted and true vertex, theta-binned
        if {"vtx_x","vtx_y","vtx_z"}.issubset(df_vtx.columns):
            df_acc_merged = df_res.merge(
                df_vtx.loc[df_vtx["accepted"]==1,
                           ["event_id","vertex_index","vtx_x","vtx_y","vtx_z"]],
                on=["event_id","vertex_index"], how="inner"
            )
            for theta_spec in THETA_RANGES_DEG + [("all","all")]:
                if theta_spec == ("all","all"):
                    sub = df_acc_merged; lbl = "all theta"
                else:
                    lo, hi = theta_spec
                    sub = df_acc_merged[(df_acc_merged["theta_deg"]>=lo) & (df_acc_merged["theta_deg"]<hi)]
                    lbl = f"{lo:.0f}–{hi:.0f} deg"
                if sub.empty: continue
                fig, axes = plt.subplots(1, 3, figsize=(12, 4))
                for ax, comp in zip(axes, ["x","y","z"]):
                    ax.hist(sub[f"delta_{comp}"].dropna(),       bins=60, histtype="step",
                            label="track - fitted vtx", density=True)
                    ax.hist(sub[f"delta_truth_{comp}"].dropna(), bins=60, histtype="step",
                            ls="--", label="track - true vtx", density=True)
                    ax.set_xlabel(f"Delta {comp} [cm]"); ax.set_ylabel("Norm. counts")
                    ax.set_title(f"{comp}-residuals, {lbl}")
                    ax.legend(); ax.grid(True, alpha=0.3)
                fig.suptitle(f"Track residuals ({lbl})", y=1.03)
                fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    print(f"[control plots] Wrote {output_pdf}")


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="HIBEAM streaming vertex finder (ROOT I/O, no OOM)."
    )
    parser.add_argument("input_root",
                        help="Input ROOT file with digitizedHits TTree.")
    parser.add_argument("--tree", dest="tree_name", default="digitizedHits")
    parser.add_argument("--output-root", dest="output_root",
                        default="vertex_products.root")
    parser.add_argument("--output-pdf", dest="output_pdf",
                        default="vertex_diagnostics.pdf")

    # TPC resolution
    parser.add_argument("--tpc-sigma-x", dest="tpc_sigma_x", type=float, default=0.05)
    parser.add_argument("--tpc-sigma-y", dest="tpc_sigma_y", type=float, default=0.05)
    parser.add_argument("--tpc-sigma-z", dest="tpc_sigma_z", type=float, default=0.10)

    # Quality cuts
    parser.add_argument("--chi2ndf-max",  dest="chi2ndf_max",  type=float, default=15.0)
    parser.add_argument("--max-dca",      dest="max_dca",       type=float, default=3.0)
    parser.add_argument("--max-r-vertex", dest="max_r_vertex",  type=float, default=20.0)

    parser.add_argument("--min-comb", dest="min_comb", type=int, default=2)
    parser.add_argument("--max-comb", dest="max_comb", type=int, default=5)
    parser.add_argument("--z-constraint", dest="z_constraint", type=float, default=0.0)

    parser.add_argument("--use-tpc-smearing", dest="use_tpc_smearing",
                        action="store_true", default=False)


    parser.add_argument("--seed",       type=int, default=None)
    parser.add_argument("--max-events", dest="max_events", type=int, default=0)
    parser.add_argument("--debug-theta90", dest="debug_theta90",
                        action="store_true", default=False)

    args = parser.parse_args()

    global DEBUG_THETA90
    DEBUG_THETA90 = args.debug_theta90

    # if args.card:
    #     apply_card_overrides(args, load_steering_card(args.card))

    if args.seed is not None:
        np.random.seed(args.seed)

    comb_sizes = list(range(args.min_comb, args.max_comb + 1))
    print("[config] comb_sizes =", comb_sizes)
    print("[config] chi2ndf_max =", args.chi2ndf_max,
          "  max_dca =", args.max_dca,
          "  max_r_vertex =", args.max_r_vertex)
    print(f"[config] tpc_sigma = ({args.tpc_sigma_x}, {args.tpc_sigma_y}, {args.tpc_sigma_z}) cm")

    # Counters
    total_events = 0
    found_counts = {n: 0 for n in comb_sizes}
    total_tracks = idx_p = idx_q = idx_r = idx_s = idx_t = 0
    event_diag: Dict[int, Dict[str, int]] = {
        n: dict(too_few_tracks=0, geN_no_accepted=0,
                fail_chi2=0, fail_dca=0, fail_r=0, other_fail=0)
        for n in comb_sizes
    }

    writer = RootWriter(args.output_root)

    # ---- Main event loop: one event in memory at a time ----
    for event_id, df_event in iter_events_from_root(args.input_root, args.tree_name):
        total_events += 1
        if args.max_events > 0 and total_events > args.max_events:
            break

        ev_result = process_event(
            df_event=df_event,
            comb_sizes=comb_sizes,
            tpc_sigma_x=args.tpc_sigma_x, tpc_sigma_y=args.tpc_sigma_y,
            tpc_sigma_z=args.tpc_sigma_z,
            chi2ndf_max=args.chi2ndf_max,  max_dca_cut=args.max_dca,
            max_r_vertex=args.max_r_vertex, z_constraint=args.z_constraint,
            use_tpc_smearing=args.use_tpc_smearing
            )

        # ---- Build per-event row lists and flush immediately ----
        vtx_rows: List[Dict] = []
        res_rows: List[Dict] = []

        for vtx_index, v in enumerate(ev_result["vertices"]):
            vtx_rows.append({
                "event_id": event_id, "vertex_index": vtx_index,
                "n_tracks": v["n_tracks"],
                "vtx_x": v["vtx_x"], "vtx_y": v["vtx_y"], "vtx_z": v["vtx_z"],
                "vtx_r": v["vtx_r"],
                "vtx_dx_true": v["vtx_dx_true"], "vtx_dy_true": v["vtx_dy_true"],
                "vtx_dz_true": v["vtx_dz_true"], "vtx_d3_true": v["vtx_d3_true"],
                "chi2": v["chi2"], "ndf": v["ndf"], "chi2ndf": v["chi2ndf"],
                "cond": v["cond"],
                "max_dca": v["max_dca"], "avg_dca": v["avg_dca"],
                "accepted_raw": int(v["accepted_raw"]), "accepted": int(v["accepted"]),
                "chi2ndf_max_used": v["chi2ndf_max_used"],
                "max_dca_cut_used": v["max_dca_cut_used"],
                "max_r_vertex_cut_used": v["max_r_vertex_cut_used"],
                "fail_chi2": int(v["fail_chi2"]), "fail_dca": int(v["fail_dca"]),
                "fail_r": int(v["fail_r"]),
            })
            for res in v["residuals"]:
                res_rows.append({
                    "event_id": event_id, "vertex_index": vtx_index,
                    "n_tracks": v["n_tracks"],
                    "track_id":  res.get("track_id", -1),
                    "is_signal": res.get("is_signal", 1),
                    "delta_x": res["delta_x"], "delta_y": res["delta_y"],
                    "delta_z": res["delta_z"],
                    "dca":          res["dca"],
                    "dca_truth":    res.get("dca_truth",      float("nan")),
                    "dca_truth_3d": res.get("dca_truth_3d",   float("nan")),
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
                })

        # Write this event's rows immediately — no accumulation
        writer.flush(vtx_rows, res_rows)

        # ---- Update counters (no row storage) ----
        n_valid = ev_result["n_valid_tracks"]
        diag_by_N = ev_result["diag_by_N"]
        for n in comb_sizes:
            if n_valid < n:
                event_diag[n]["too_few_tracks"] += 1
            else:
                if not ev_result["found_by_N"].get(n, False):
                    event_diag[n]["geN_no_accepted"] += 1
                    dN = diag_by_N.get(n, {})
                    if dN.get("any_chi2_fail"): event_diag[n]["fail_chi2"] += 1
                    if dN.get("any_dca_fail"):  event_diag[n]["fail_dca"]  += 1
                    if dN.get("any_r_fail"):    event_diag[n]["fail_r"]    += 1
                    if dN.get("any_combo_considered") and not (
                        dN.get("any_chi2_fail") or dN.get("any_dca_fail") or dN.get("any_r_fail")
                    ):
                        event_diag[n]["other_fail"] += 1
            if ev_result["found_by_N"].get(n, False):
                found_counts[n] += 1

        for _, meta in ev_result["track_meta"].items():
            total_tracks += 1
            if meta["invalid_exit"] and not meta["used_in_any_accepted"]:
                idx_p += 1
            if (not meta["invalid_exit"]) and (not meta["used_in_any_comb"]) \
                    and (not meta["used_in_any_accepted"]):
                idx_q += 1
            if meta["used_in_any_comb"] and not meta["used_in_any_accepted"]:
                if meta["in_combo_fail_chi2"]: idx_r += 1
                if meta["in_combo_fail_dca"]:  idx_s += 1
                if meta["in_combo_fail_r"]:    idx_t += 1

        if total_events % 500 == 0:
            print(f"  ... processed {total_events} events")

    writer.close()

    # ---- Summary ----
    print(f"\n=== Summary: {total_events} events processed ===")
    for n in comb_sizes:
        frac = found_counts[n] / total_events if total_events else 0.0
        print(f"  Accepted {n}-track vertices: {found_counts[n]}/{total_events} ({frac:.3f})")
    print(f"\nTrack indices: total={total_tracks}, p={idx_p}, q={idx_q}, "
          f"r={idx_r}, s={idx_s}, t={idx_t}")
    for n in comb_sizes:
        d = event_diag[n]
        print(f"\nN={n}: too_few={d['too_few_tracks']}, no_accepted={d['geN_no_accepted']} "
              f"[chi2={d['fail_chi2']}, dca={d['fail_dca']}, r={d['fail_r']}, "
              f"other={d['other_fail']}]")

    # ---- Control plots ----
    try:
        make_control_plots(args.output_root, args.output_pdf)
    except Exception as exc:
        print(f"[control plots] Warning: {exc}")
        import traceback; traceback.print_exc()


if __name__ == "__main__":
    main()

# #!/usr/bin/env python3
# """
# hibeam_vertex_finder_streaming.py  —  Route-1 rewrite

# Mathematical basis
# ------------------
# The chi2 is formulated entirely in **hit space**, not DCA space, so that
# the classical theorem chi2 ~ chi2(ndf) applies exactly.

# Model
# -----
# Each track i is represented by a line through the vertex v = (x, y, 0):

#     q_ij(s) = v + s * u_i          j = entry (j=0), exit (j=1)

# The observed hit positions p_ij are drawn from

#     p_ij ~ N( q_ij(s_ij*), Sigma )     Sigma = diag(sx^2, sy^2, sz^2)

# where s_ij* is the true parameter value.  The chi2 is therefore

#     chi2 = sum_{i,j} (p_ij - v - s_ij u_i)^T Sigma^{-1} (p_ij - v - s_ij u_i)

# Free parameters
# ---------------
#     v   : (x, y)                                  — 2 parameters
#     u_i : unit direction (2 angles per track)     — 2N parameters
#     s_i0, s_i1 : scalars for entry/exit           — 2N parameters
#     ----------------------------------------------------------------
#     Total k = 2 + 4N

# Observations
# ------------
#     2 points x 3 components x N tracks = 6N scalar observations

# ndf = 6N - (2 + 4N) = 2N - 2

# This is identical to the old ndf formula but now **derived from first
# principles**, so chi2/ndf ~ 1 at the minimum when the model is correct.

# Fitting strategy
# ----------------
# The problem is nonlinear (u_i must be a unit vector) but separable.
# For fixed v the optimal (u_i, s_i0, s_i1) have closed-form solutions,
# and for fixed (u_i) the optimal v is a linear system.  We therefore
# use coordinate descent (alternating minimisation):

#     1. Initialise u_i from raw entry/exit points.
#     2. Solve for v = (x, y) via the standard linear vertex fit.
#     3. For each track refit u_i as the principal axis of the two
#        centred points (p_i0 - v, p_i1 - v).  This is the eigenvector
#        of the 2-point scatter matrix corresponding to the largest
#        eigenvalue.  Recompute s_ij as projections onto u_i.
#     4. Repeat steps 2-3 until |delta_v| < tol (typically < 5 iter).
#     5. Evaluate chi2 in hit space with known Sigma.

# All other logic (quality cuts, ROOT I/O, control plots, steering card,
# MS sigma, …) is unchanged from the original.
# """

# import argparse
# import itertools
# import math
# from typing import List, Tuple, Dict, Any, Optional

# import numpy as np
# import pandas as pd
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages

# import uproot
# import awkward as ak


# DEBUG_THETA90 = False

# # -----------------------------------------------------------------------
# # Streaming ROOT reader
# # -----------------------------------------------------------------------

# _ROOT_BRANCHES = [
#     "TPC_firstPosX", "TPC_firstPosY", "TPC_firstPosZ",
#     "TPC_lastPosX",  "TPC_lastPosY",  "TPC_lastPosZ",
#     "TPC_pdg",
#     "PrimaryPosX", "PrimaryPosY", "PrimaryPosZ",
# ]


# def iter_events_from_root(root_path: str, tree_name: str = "digitizedHits"):
#     """Yield (event_id, DataFrame) one TTree entry at a time."""
#     with uproot.open(root_path) as f:
#         tree = f[tree_name]
#         for batch_idx, arrays in enumerate(
#             tree.iterate(_ROOT_BRANCHES, library="ak", step_size=1)
#         ):
#             event_id = batch_idx

#             ex = arrays["TPC_firstPosX"][0]
#             ey = arrays["TPC_firstPosY"][0]
#             ez = arrays["TPC_firstPosZ"][0]
#             lx = arrays["TPC_lastPosX"][0]
#             ly = arrays["TPC_lastPosY"][0]
#             lz = arrays["TPC_lastPosZ"][0]
#             path = np.sqrt((lx - ex)**2 + (ly - ey)**2 + (lz - ez)**2)
#             pdg  = arrays["TPC_pdg"][0]

#             pvx = arrays["PrimaryPosX"][0]
#             pvy = arrays["PrimaryPosY"][0]
#             pvz = arrays["PrimaryPosZ"][0]

#             n_tracks = len(ex)
#             if n_tracks == 0:
#                 continue

#             truth_x = float(pvx[0]) if len(pvx) > 0 else float("nan")
#             truth_y = float(pvy[0]) if len(pvy) > 0 else float("nan")
#             truth_z = float(pvz[0]) if len(pvz) > 0 else float("nan")

#             ex_np = ak.to_numpy(ex).astype(float)
#             ey_np = ak.to_numpy(ey).astype(float)
#             ez_np = ak.to_numpy(ez).astype(float)
#             lx_np = ak.to_numpy(lx).astype(float)
#             ly_np = ak.to_numpy(ly).astype(float)
#             lz_np = ak.to_numpy(lz).astype(float)
#             path_np = ak.to_numpy(path).astype(float)

#             finite_entry = np.isfinite(ex_np) & np.isfinite(ey_np) & np.isfinite(ez_np)
#             finite_exit  = np.isfinite(lx_np) & np.isfinite(ly_np) & np.isfinite(lz_np)
#             gas_exited   = ((path_np > 0) & finite_entry & finite_exit).astype(np.int32)

#             df = pd.DataFrame({
#                 "event_id":        np.full(n_tracks, event_id, dtype=np.int64),
#                 "track_id":        np.arange(n_tracks, dtype=np.int64),
#                 "entry_x":         ex_np, "entry_y": ey_np, "entry_z": ez_np,
#                 "exit_x":          lx_np, "exit_y":  ly_np, "exit_z":  lz_np,
#                 "gas_exited_flag": gas_exited,
#                 "truth_vtx_x":     np.full(n_tracks, truth_x),
#                 "truth_vtx_y":     np.full(n_tracks, truth_y),
#                 "truth_vtx_z":     np.full(n_tracks, truth_z),
#                 "is_signal":       np.ones(n_tracks, dtype=np.int32),
#                 "TPC_PathLength":  path_np,
#                 "TPC_pdg":         ak.to_numpy(pdg).astype(np.int32),
#             })

#             yield event_id, df


# # -----------------------------------------------------------------------
# # Incremental ROOT writer  (unchanged)
# # -----------------------------------------------------------------------

# class RootWriter:
#     _VTX_DTYPES: Dict[str, Any] = {
#         "event_id": np.int64, "vertex_index": np.int64, "n_tracks": np.int64,
#         "vtx_x": np.float64, "vtx_y": np.float64, "vtx_z": np.float64,
#         "vtx_r": np.float64,
#         "vtx_dx_true": np.float64, "vtx_dy_true": np.float64,
#         "vtx_dz_true": np.float64, "vtx_d3_true": np.float64,
#         "chi2": np.float64, "ndf": np.int64, "chi2ndf": np.float64,
#         "cond": np.float64, "max_dca": np.float64, "avg_dca": np.float64,
#         "accepted_raw": np.int32, "accepted": np.int32,
#         "chi2ndf_max_used": np.float64, "max_dca_cut_used": np.float64,
#         "max_r_vertex_cut_used": np.float64,
#         "fail_chi2": np.int32, "fail_dca": np.int32, "fail_r": np.int32,
#     }

#     _RES_DTYPES: Dict[str, Any] = {
#         "event_id": np.int64, "vertex_index": np.int64, "n_tracks": np.int64,
#         "track_id": np.int64, "is_signal": np.int32,
#         "delta_x": np.float64, "delta_y": np.float64, "delta_z": np.float64,
#         "dca": np.float64,
#         "dca_truth": np.float64, "dca_truth_3d": np.float64,
#         "dca_plane_truth": np.float64,
#         "delta_truth_x": np.float64, "delta_truth_y": np.float64,
#         "delta_truth_z": np.float64,
#         "sigma_x": np.float64, "sigma_y": np.float64, "sigma_z": np.float64,
#         "p_x": np.float64, "p_y": np.float64, "p_z": np.float64,
#         "u_x": np.float64, "u_y": np.float64, "u_z": np.float64,
#     }

#     def __init__(self, path: str):
#         self._path = path
#         self._file = uproot.recreate(path)
#         self._vtx_created = False
#         self._res_created = False

#     def _rows_to_arrays(self, rows, dtype_map):
#         out = {}
#         for col, dtype in dtype_map.items():
#             vals = [r.get(col, 0 if np.issubdtype(dtype, np.integer) else float("nan"))
#                     for r in rows]
#             out[col] = np.asarray(vals, dtype=dtype)
#         return out

#     def flush(self, vtx_rows, res_rows):
#         if vtx_rows:
#             arrays = self._rows_to_arrays(vtx_rows, self._VTX_DTYPES)
#             if not self._vtx_created:
#                 self._file.mktree("vertices", {k: v.dtype for k, v in arrays.items()})
#                 self._vtx_created = True
#             self._file["vertices"].extend(arrays)
#         if res_rows:
#             arrays = self._rows_to_arrays(res_rows, self._RES_DTYPES)
#             if not self._res_created:
#                 self._file.mktree("residuals", {k: v.dtype for k, v in arrays.items()})
#                 self._res_created = True
#             self._file["residuals"].extend(arrays)

#     def close(self):
#         self._file.close()
#         print(f"[RootWriter] Closed {self._path}")


# # -----------------------------------------------------------------------
# # Multiple scattering  (unchanged)
# # -----------------------------------------------------------------------

# def estimate_ms_sigma_for_track(
#     entry_point, direction,
#     ms_density, ms_radlen_gcm2, ms_thickness_cm,
#     ms_lever_arm_cm, ms_p_MeV, ms_min_cos,
# ) -> float:
#     if ms_density <= 0 or ms_radlen_gcm2 <= 0 or ms_thickness_cm <= 0 or ms_p_MeV <= 0:
#         return 0.0
#     ex, ey = entry_point[0], entry_point[1]
#     r_xy = math.hypot(ex, ey)
#     cos_inc = abs(float(np.dot(direction, np.array([ex, ey, 0.0]) / r_xy))) if r_xy > 0 else 1.0
#     cos_eff = max(cos_inc, ms_min_cos)
#     x_over_X0 = ms_density * ms_thickness_cm / cos_eff / ms_radlen_gcm2
#     if x_over_X0 <= 0.0:
#         return 0.0
#     theta0 = (13.6 / ms_p_MeV) * math.sqrt(x_over_X0) * (1.0 + 0.038 * math.log(x_over_X0))
#     return abs(theta0 * ms_lever_arm_cm)


# # -----------------------------------------------------------------------
# # Route-1 core: joint fit of vertex + track directions in hit space
# # -----------------------------------------------------------------------

# def _refit_direction(p_entry: np.ndarray, p_exit: np.ndarray,
#                      vtx: np.ndarray) -> np.ndarray:
#     """
#     For fixed vertex v, find the unit direction u_i that minimises
#     the sum of squared distances from p_entry and p_exit to the line
#     through v along u_i (unweighted, i.e. sigma_x = sigma_y = sigma_z).

#     This is the eigenvector of the largest eigenvalue of the scatter
#     matrix  M = (a a^T + b b^T)  where a = p_entry - v, b = p_exit - v.

#     For the weighted case (sigma_x != sigma_y != sigma_z) replace M
#     with  Sigma^{-1/2} M Sigma^{-1/2}  — see weighted variant below.
#     """
#     a = p_entry - vtx
#     b = p_exit  - vtx
#     M = np.outer(a, a) + np.outer(b, b)
#     # eigenvalues in ascending order
#     eigvals, eigvecs = np.linalg.eigh(M)
#     u = eigvecs[:, -1]          # largest eigenvalue → best-fit direction
#     return u / np.linalg.norm(u)


# def _refit_direction_weighted(p_entry: np.ndarray, p_exit: np.ndarray,
#                                vtx: np.ndarray,
#                                sigma: np.ndarray) -> np.ndarray:
#     """
#     Weighted version: minimise
#         (a - s*u)^T Sigma^{-1} (a - s*u) + (b - t*u)^T Sigma^{-1} (b - t*u)
#     over unit vector u (and scalars s, t).

#     Substituting the optimal s = u^T Sigma^{-1} a / (u^T Sigma^{-1} u)
#     and similarly for t, and maximising over u, leads to the
#     generalised eigenproblem  M_w u = lambda u  where

#         M_w = Sigma^{-1/2} (a a^T + b b^T) Sigma^{-1/2}

#     and u is expressed in the whitened basis.  We work in whitened
#     coordinates throughout.
#     """
#     inv_sigma = 1.0 / sigma                      # shape (3,)
#     sqrt_inv  = np.sqrt(inv_sigma)               # Sigma^{-1/2} (diagonal)

#     a_w = sqrt_inv * (p_entry - vtx)
#     b_w = sqrt_inv * (p_exit  - vtx)

#     M_w = np.outer(a_w, a_w) + np.outer(b_w, b_w)
#     eigvals, eigvecs = np.linalg.eigh(M_w)
#     u_w = eigvecs[:, -1]                         # in whitened space

#     # back-transform to original space
#     u_orig = sqrt_inv * u_w                      # Sigma^{-1/2} u_w
#     norm   = np.linalg.norm(u_orig)
#     if norm < 1e-12:
#         # fallback: unweighted direction
#         return _refit_direction(p_entry, p_exit, vtx)
#     return u_orig / norm


# def _fit_vertex_from_directions(
#     points: List[np.ndarray],
#     dirs:   List[np.ndarray],
#     sigma:  np.ndarray,
#     z_constraint: float = 0.0,
# ) -> Tuple[np.ndarray, float]:
#     """
#     Standard linear vertex fit with z fixed to z_constraint.
#     Weights each track by its contribution in the whitened metric.

#     Returns (vtx, condition_number).
#     """
#     inv_sigma2 = 1.0 / sigma**2
#     I3 = np.eye(3)
#     A  = np.zeros((3, 3))
#     b  = np.zeros(3)

#     for p, u in zip(points, dirs):
#         if not np.all(np.isfinite(u)):
#             continue
#         # Weighted projection: P_i = Sigma^{-1} (I - u u^T)
#         P = (I3 - np.outer(u, u)) * inv_sigma2[np.newaxis, :]
#         A += P
#         b += P @ p

#     z_c  = float(z_constraint)
#     A2   = A[:2, :2]
#     b2   = b[:2] - A[:2, 2] * z_c

#     try:
#         xy   = np.linalg.solve(A2, b2)
#         cond = float(np.linalg.cond(A2))
#         x, y = float(xy[0]), float(xy[1])
#     except np.linalg.LinAlgError:
#         x = y = float("nan")
#         cond = float("inf")

#     return np.array([x, y, z_c]), cond


# def fit_vertex_route1(
#     tracks:        List[Dict],
#     sigma:         np.ndarray,
#     z_constraint:  float = 0.0,
#     max_iter:      int   = 10,
#     tol:           float = 1e-6,
# ) -> Tuple[np.ndarray, np.ndarray, float, float, int]:
#     """
#     Joint fit of vertex position and track directions in hit space.

#     Parameters
#     ----------
#     tracks       : list of dicts with keys 'p_entry', 'p_exit'
#                    (raw smeared hit positions), plus metadata.
#     sigma        : 1-D array [sx, sy, sz] — known hit resolutions.
#     z_constraint : fixed z value of the vertex plane.
#     max_iter     : maximum alternating-minimisation iterations.
#     tol          : convergence criterion on |delta_v| [cm].

#     Returns
#     -------
#     vtx      : fitted vertex (x, y, z_constraint)
#     dirs     : array of fitted unit directions, shape (N, 3)
#     chi2     : chi2 evaluated in hit space  (scalar)
#     cond     : condition number of the 2x2 vertex fit matrix
#     ndf      : 2N - 2   (derived from first principles, see module docstring)
#     """
#     N = len(tracks)
#     if N == 0:
#         return np.full(3, np.nan), np.full((0, 3), np.nan), float("nan"), float("inf"), 0

#     p_entries = np.array([t["p_entry"] for t in tracks])   # (N, 3)
#     p_exits   = np.array([t["p_exit"]  for t in tracks])   # (N, 3)

#     # Step 1: initialise directions from raw entry/exit vectors
#     raw_diff = p_exits - p_entries                          # (N, 3)
#     norms    = np.linalg.norm(raw_diff, axis=1, keepdims=True)
#     norms    = np.where(norms > 1e-12, norms, 1.0)
#     dirs     = raw_diff / norms                             # (N, 3)

#     # Seed vertex: midpoint of DCA-based estimate using initial directions
#     vtx, cond = _fit_vertex_from_directions(
#         list(p_entries), list(dirs), sigma, z_constraint
#     )

#     # Alternating minimisation
#     for _ in range(max_iter):
#         # Step 2: refit directions given current vertex
#         new_dirs = np.zeros_like(dirs)
#         for i, trk in enumerate(tracks):
#             new_dirs[i] = _refit_direction_weighted(
#                 p_entries[i], p_exits[i], vtx, sigma
#             )
#         dirs = new_dirs

#         # Step 3: refit vertex given new directions
#         vtx_new, cond = _fit_vertex_from_directions(
#             list(p_entries), list(dirs), sigma, z_constraint
#         )
#         if np.linalg.norm(vtx_new - vtx) < tol:
#             vtx = vtx_new
#             break
#         vtx = vtx_new

#     # ----------------------------------------------------------------
#     # Evaluate chi2 in hit space
#     # ----------------------------------------------------------------
#     # For each track i and each endpoint j in {entry, exit}:
#     #   residual_ij = p_ij - (vtx + s_ij * u_i)
#     # where the optimal  s_ij = u_i . (p_ij - vtx)  (projection).
#     #
#     # Because  residual_ij = (I - u_i u_i^T)(p_ij - vtx)
#     # this is just the component of (p_ij - vtx) perpendicular to u_i.
#     #
#     # chi2 = sum_{i,j} residual_ij^T Sigma^{-1} residual_ij
#     # ----------------------------------------------------------------
#     inv_sigma2 = 1.0 / sigma**2
#     chi2 = 0.0
#     for i, u in enumerate(dirs):
#         if not np.all(np.isfinite(u)):
#             continue
#         for p in (p_entries[i], p_exits[i]):
#             d    = p - vtx
#             proj = np.dot(u, d) * u          # component along track
#             res  = d - proj                  # component perpendicular to track
#             chi2 += float(np.dot(res**2, inv_sigma2))

#     # ndf from first principles:  6N observations - (2 + 4N) parameters = 2N - 2
#     ndf = 2 * N - 2

#     return vtx, dirs, chi2, cond, ndf


# # -----------------------------------------------------------------------
# # Residuals and truth matching  (adapted for Route-1 fitted directions)
# # -----------------------------------------------------------------------

# def compute_residuals_route1(
#     vtx:    np.ndarray,
#     tracks: List[Dict],
#     dirs:   np.ndarray,
#     sigma:  np.ndarray,
# ) -> Tuple[List[Dict], float, int]:
#     """
#     Compute per-track hit-space residuals and the total chi2.

#     The chi2 here is identical to what fit_vertex_route1 returns —
#     it is recalculated here so that per-track contributions are available
#     for diagnostics.

#     Also computes the DCA (distance of closest approach in the xy plane)
#     for use in the max_dca quality cut, using the *fitted* direction.
#     """
#     residuals = []
#     chi2      = 0.0
#     inv_sigma2 = 1.0 / sigma**2

#     for i, trk in enumerate(tracks):
#         u = dirs[i]
#         if not np.all(np.isfinite(u)):
#             continue

#         p_entry = trk["p_entry"]
#         p_exit  = trk["p_exit"]

#         trk_chi2 = 0.0
#         for p in (p_entry, p_exit):
#             d    = p - vtx
#             proj = np.dot(u, d) * u
#             res  = d - proj
#             trk_chi2 += float(np.dot(res**2, inv_sigma2))
#         chi2 += trk_chi2

#         # DCA in the xy plane: distance from vtx to line through p_entry
#         # along fitted direction u (for the max_dca quality cut)
#         d_entry = p_entry - vtx
#         dca_vec = d_entry - np.dot(u, d_entry) * u
#         dca_xy  = float(math.hypot(dca_vec[0], dca_vec[1]))
#         dca_3d  = float(np.linalg.norm(dca_vec))

#         residuals.append({
#             "track_id":  int(trk.get("track_id", -1)),
#             "is_signal": int(trk.get("is_signal", 1)),
#             # residual of entry point perpendicular to fitted direction
#             "delta_x": float(dca_vec[0]),
#             "delta_y": float(dca_vec[1]),
#             "delta_z": float(dca_vec[2]),
#             "dca":     dca_xy,
#             "dca_3d":  dca_3d,
#             "sigma_x": float(sigma[0]),
#             "sigma_y": float(sigma[1]),
#             "sigma_z": float(sigma[2]),
#             # store fitted track reference point and direction for truth matching
#             "p_x": float(p_entry[0]),
#             "p_y": float(p_entry[1]),
#             "p_z": float(p_entry[2]),
#             "u_x": float(u[0]),
#             "u_y": float(u[1]),
#             "u_z": float(u[2]),
#         })

#     ndf = 2 * len(tracks) - 2
#     return residuals, chi2, ndf


# def add_truth_residuals(
#     residuals:  List[Dict],
#     tracks:     List[Dict],
#     dirs:       np.ndarray,
#     truth_vtx:  Optional[np.ndarray],
# ) -> None:
#     """Add truth-matching columns to residual dicts (in-place)."""
#     if truth_vtx is None:
#         nan = float("nan")
#         for res in residuals:
#             for k in ("dca_truth", "dca_truth_3d", "dca_plane_truth",
#                       "delta_truth_x", "delta_truth_y", "delta_truth_z"):
#                 res[k] = nan
#         return

#     for res, trk, u in zip(residuals, tracks, dirs):
#         p = trk["p_entry"]
#         if not np.all(np.isfinite(u)):
#             for k in ("dca_truth", "dca_truth_3d", "dca_plane_truth",
#                       "delta_truth_x", "delta_truth_y", "delta_truth_z"):
#                 res[k] = float("nan")
#             continue
#         t_true       = float(np.dot(u, truth_vtx - p))
#         delta_true   = (p + t_true * u) - truth_vtx
#         res["dca_truth"]       = float(math.hypot(delta_true[0], delta_true[1]))
#         res["dca_truth_3d"]    = float(np.linalg.norm(delta_true))
#         res["delta_truth_x"]   = float(delta_true[0])
#         res["delta_truth_y"]   = float(delta_true[1])
#         res["delta_truth_z"]   = float(delta_true[2])
#         uz = float(u[2])
#         if abs(uz) > 1e-6:
#             t_pl = (truth_vtx[2] - p[2]) / uz
#             r_pl = p + t_pl * u
#             res["dca_plane_truth"] = float(
#                 math.hypot(r_pl[0] - truth_vtx[0], r_pl[1] - truth_vtx[1])
#             )
#         else:
#             res["dca_plane_truth"] = float("nan")


# # -----------------------------------------------------------------------
# # Steering card  (unchanged)
# # -----------------------------------------------------------------------

# def load_steering_card(path: str) -> Dict[str, str]:
#     params: Dict[str, str] = {}
#     with open(path) as f:
#         for line in f:
#             line = line.strip()
#             if not line or line.startswith("#"):
#                 continue
#             if "#" in line:
#                 line = line.split("#", 1)[0].rstrip()
#             if "=" not in line:
#                 continue
#             k, v = line.split("=", 1)
#             params[k.strip()] = v.strip()
#     return params


# def apply_card_overrides(args, card_params: Dict[str, str]) -> None:
#     float_keys = {
#         "tpc_sigma_x", "tpc_sigma_y", "tpc_sigma_z",
#         "chi2ndf_max", "max_dca", "max_r_vertex", "z_constraint",
#         "ms_density", "ms_radlen_gcm2", "ms_thickness_cm",
#         "ms_lever_arm_cm", "ms_p_MeV", "ms_min_cos",
#     }
#     int_keys  = {"min_comb", "max_comb", "seed", "max_events"}
#     bool_keys = {"use_ms_sigma", "use_tpc_smearing"}

#     for key, val in card_params.items():
#         if not hasattr(args, key):
#             continue
#         try:
#             if key in bool_keys:
#                 new_value = val.lower() in ("1", "true", "yes", "on")
#             elif key in int_keys:
#                 new_value = int(val)
#             elif key in float_keys:
#                 new_value = float(val)
#             else:
#                 new_value = val
#         except ValueError:
#             print(f"[card] Could not parse {key}={val!r}, keeping {getattr(args,key)!r}")
#             continue
#         print(f"[card] {key}: {getattr(args,key)!r} -> {new_value!r}")
#         setattr(args, key, new_value)


# # -----------------------------------------------------------------------
# # Per-event processing
# # -----------------------------------------------------------------------

# def process_event(
#     df_event:        pd.DataFrame,
#     comb_sizes:      List[int],
#     tpc_sigma_x:     float, tpc_sigma_y: float, tpc_sigma_z: float,
#     chi2ndf_max:     float, max_dca_cut: float, max_r_vertex: float,
#     z_constraint:    float,
#     use_tpc_smearing: bool, use_ms_sigma: bool,
#     ms_density:      float, ms_radlen_gcm2: float, ms_thickness_cm: float,
#     ms_lever_arm_cm: float, ms_p_MeV: float, ms_min_cos: float,
# ) -> Dict[str, Any]:

#     track_meta: Dict[int, Dict] = {}
#     for _, row in df_event.iterrows():
#         tid = int(row["track_id"])
#         track_meta[tid] = {
#             "gas_exited":               int(row["gas_exited_flag"]) == 1,
#             "invalid_exit":             int(row["gas_exited_flag"]) != 1,
#             "used_in_any_comb":         False,
#             "used_in_any_accepted_raw": False,
#             "used_in_any_accepted":     False,
#             "in_combo_fail_chi2":       False,
#             "in_combo_fail_dca":        False,
#             "in_combo_fail_r":          False,
#             "is_signal":                int(row.get("is_signal", 1)),
#         }

#     df_valid = df_event[df_event["gas_exited_flag"] == 1]
#     diag_by_N = {
#         n: dict(any_combo_considered=False, any_accepted=False,
#                 any_chi2_fail=False, any_dca_fail=False, any_r_fail=False)
#         for n in comb_sizes
#     }

#     if df_valid.empty:
#         return {
#             "vertices": [], "found_by_N": {n: False for n in comb_sizes},
#             "track_meta": track_meta, "n_valid_tracks": 0, "diag_by_N": diag_by_N,
#         }

#     # Build track objects — store raw hit positions (p_entry, p_exit)
#     # plus metadata.  Smearing is applied here if requested, exactly
#     # as in the original code.
#     sigma_base = np.array([tpc_sigma_x, tpc_sigma_y, tpc_sigma_z])

#     tracks: List[Dict] = []
#     for _, row in df_valid.iterrows():
#         pdg = int(row["TPC_pdg"])
#         if abs(pdg) == 11:
#             continue

#         entry = np.array([row["entry_x"], row["entry_y"], row["entry_z"]], dtype=float)
#         exit_ = np.array([row["exit_x"],  row["exit_y"],  row["exit_z"]],  dtype=float)

#         # Optional additional smearing (use_tpc_smearing adds noise on
#         # top of the G4 smearing — keep disabled if G4 already smears)
#         if use_tpc_smearing:
#             entry = entry + np.random.normal(0.0, sigma_base)
#             exit_ = exit_  + np.random.normal(0.0, sigma_base)

#         # Effective sigma: base resolution + MS contribution
#         sigma = sigma_base.copy()
#         if use_ms_sigma:
#             u_init = exit_ - entry
#             n      = np.linalg.norm(u_init)
#             if n > 1e-12:
#                 u_init /= n
#             sms = estimate_ms_sigma_for_track(
#                 entry, u_init,
#                 ms_density, ms_radlen_gcm2, ms_thickness_cm,
#                 ms_lever_arm_cm, ms_p_MeV, ms_min_cos,
#             )
#             sigma = np.sqrt(sigma**2 + sms**2)

#         tracks.append({
#             "track_id": int(row["track_id"]),
#             "p_entry":  entry,
#             "p_exit":   exit_,
#             "sigma":    sigma,          # per-track effective resolution
#             "is_signal": bool(row.get("is_signal", 1)),
#         })

#     truth_x = df_event["truth_vtx_x"].iloc[0]
#     truth_y = df_event["truth_vtx_y"].iloc[0]
#     truth_z = df_event["truth_vtx_z"].iloc[0]
#     truth_vtx = (np.array([truth_x, truth_y, truth_z])
#                  if np.all(np.isfinite([truth_x, truth_y, truth_z])) else None)

#     results_vertices: List[Dict] = []
#     found_by_N = {n: False for n in comb_sizes}

#     for n in comb_sizes:
#         if len(tracks) < n:
#             continue
#         for combo in itertools.combinations(tracks, n):
#             combo_tracks = list(combo)
#             track_ids    = [t["track_id"] for t in combo_tracks]

#             # Use the average sigma across tracks in the combination.
#             # If all tracks share the same sigma_base this is exact;
#             # with MS corrections it is an approximation.
#             combo_sigma = np.mean([t["sigma"] for t in combo_tracks], axis=0)

#             # ---- Route-1 joint fit ----
#             vtx, fitted_dirs, chi2, cond, ndf = fit_vertex_route1(
#                 combo_tracks, combo_sigma, z_constraint=z_constraint
#             )

#             diag_by_N[n]["any_combo_considered"] = True
#             if not np.all(np.isfinite(vtx)):
#                 continue

#             r_vtx = math.hypot(vtx[0], vtx[1])

#             # Per-track residuals (for diagnostics and max_dca cut)
#             residuals, _, _ = compute_residuals_route1(
#                 vtx, combo_tracks, fitted_dirs, combo_sigma
#             )
#             add_truth_residuals(residuals, combo_tracks, fitted_dirs, truth_vtx)

#             if ndf <= 0 or not math.isfinite(chi2):
#                 continue

#             chi2ndf    = chi2 / ndf
#             max_dca_val = max(r["dca"] for r in residuals) if residuals else float("nan")
#             avg_dca_val = (sum(r["dca"] for r in residuals) / len(residuals)
#                            if residuals else float("nan"))

#             chi2_ok = (chi2ndf <= chi2ndf_max) if chi2ndf_max >= 0 else True
#             dca_ok  = (max_dca_val <= max_dca_cut) if max_dca_cut >= 0 else True
#             r_ok    = (r_vtx <= max_r_vertex)      if max_r_vertex >= 0 else True
#             accepted = chi2_ok and dca_ok and r_ok

#             if accepted:
#                 found_by_N[n] = True
#                 diag_by_N[n]["any_accepted"] = True
#             else:
#                 if not chi2_ok: diag_by_N[n]["any_chi2_fail"] = True
#                 if not dca_ok:  diag_by_N[n]["any_dca_fail"]  = True
#                 if not r_ok:    diag_by_N[n]["any_r_fail"]    = True

#             for tid in track_ids:
#                 m = track_meta.get(tid)
#                 if m is None: continue
#                 m["used_in_any_comb"] = True
#                 if accepted:
#                     m["used_in_any_accepted_raw"] = True
#                 else:
#                     if not chi2_ok: m["in_combo_fail_chi2"] = True
#                     if not dca_ok:  m["in_combo_fail_dca"]  = True
#                     if not r_ok:    m["in_combo_fail_r"]    = True

#             dv = vtx - truth_vtx if truth_vtx is not None else np.full(3, float("nan"))
#             results_vertices.append({
#                 "n_tracks":  n,
#                 "track_ids": track_ids,
#                 "vtx_x": float(vtx[0]), "vtx_y": float(vtx[1]), "vtx_z": float(vtx[2]),
#                 "vtx_r": float(r_vtx),
#                 "vtx_dx_true": float(dv[0]), "vtx_dy_true": float(dv[1]),
#                 "vtx_dz_true": float(dv[2]), "vtx_d3_true": float(np.linalg.norm(dv)),
#                 "chi2": float(chi2), "ndf": int(ndf), "chi2ndf": float(chi2ndf),
#                 "cond": float(cond),
#                 "max_dca": float(max_dca_val), "avg_dca": float(avg_dca_val),
#                 "accepted_raw": bool(accepted), "accepted": False,
#                 "chi2ndf_max_used":      float(chi2ndf_max),
#                 "max_dca_cut_used":      float(max_dca_cut),
#                 "max_r_vertex_cut_used": float(max_r_vertex),
#                 "fail_chi2": bool(not chi2_ok),
#                 "fail_dca":  bool(not dca_ok),
#                 "fail_r":    bool(not r_ok),
#                 "residuals": residuals,
#             })

#     # Overlap resolution: greedy, best chi2ndf first, no shared tracks
#     used_tracks: set = set()
#     for v in results_vertices:
#         v["accepted"] = False
#     for v in sorted([v for v in results_vertices if v["accepted_raw"]],
#                     key=lambda v: (-v["n_tracks"], v["chi2ndf"])):
#         if any(tid in used_tracks for tid in v["track_ids"]):
#             continue
#         v["accepted"] = True
#         used_tracks.update(v["track_ids"])

#     for m in track_meta.values():
#         m["used_in_any_accepted"] = False
#     for v in results_vertices:
#         if v["accepted"]:
#             for tid in v["track_ids"]:
#                 if tid in track_meta:
#                     track_meta[tid]["used_in_any_accepted"] = True

#     return {
#         "vertices":      results_vertices,
#         "found_by_N":    found_by_N,
#         "track_meta":    track_meta,
#         "n_valid_tracks": len(tracks),
#         "diag_by_N":     diag_by_N,
#     }


# # -----------------------------------------------------------------------
# # Control plots  (unchanged except chi2/ndf label clarification)
# # -----------------------------------------------------------------------

# THETA_RANGES_DEG = [
#     (40.0, 75.0), (75.0, 100.0), (100.0, 125.0), (125.0, 140.0), (87.0, 93.0),
# ]


# def _compute_theta_deg(ux, uy, uz) -> np.ndarray:
#     norm = np.sqrt(ux**2 + uy**2 + uz**2)
#     safe = norm > 0
#     cos_theta = np.where(safe, uz / np.where(safe, norm, 1.0), 0.0)
#     return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))


# def make_control_plots(output_root: str, output_pdf: str) -> None:
#     print(f"[control plots] Reading from {output_root} ...")

#     try:
#         with uproot.open(output_root) as f:
#             df_vtx = f["vertices"].arrays(library="pd")
#     except Exception as exc:
#         print(f"[control plots] Cannot read 'vertices' tree: {exc}"); return

#     try:
#         with uproot.open(output_root) as f:
#             chunks = [c for c in f["residuals"].iterate(library="pd", step_size=50_000)]
#             df_res = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
#     except Exception as exc:
#         print(f"[control plots] Cannot read 'residuals' tree: {exc}"); return

#     if df_vtx.empty:
#         print("[control plots] vertices tree is empty."); return

#     if "vertex_index" not in df_vtx.columns:
#         df_vtx["vertex_index"] = df_vtx.groupby("event_id").cumcount()
#     if not df_res.empty and "vertex_index" not in df_res.columns:
#         df_res["vertex_index"] = df_res.groupby(["event_id","n_tracks"]).cumcount()

#     with PdfPages(output_pdf) as pdf:

#         total_events = df_vtx["event_id"].nunique()
#         df_acc = df_vtx[df_vtx["accepted"] == 1]

#         labels = ["any"]
#         fracs  = [df_acc["event_id"].nunique() / total_events if total_events else 0.0]
#         for N in [2, 3, 4]:
#             ev_N = df_acc.loc[df_acc["n_tracks"] == N, "event_id"].nunique()
#             fracs.append(ev_N / total_events if total_events else 0.0)
#             labels.append(f"N={N}")
#         ev_gt4 = df_acc.loc[df_acc["n_tracks"] > 4, "event_id"].nunique()
#         fracs.append(ev_gt4 / total_events if total_events else 0.0)
#         labels.append("N>4")

#         fig, ax = plt.subplots(figsize=(6, 4))
#         ax.bar(range(len(labels)), fracs)
#         ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels)
#         ax.set_ylim(0, 1.05); ax.set_ylabel("Fraction of events")
#         ax.set_title("Fractions of events with accepted vertices")
#         for i, f in enumerate(fracs):
#             ax.text(i, f + 0.01, f"{f:.3f}", ha="center", fontsize=8)
#         ax.grid(True, axis="y", alpha=0.3); fig.tight_layout()
#         pdf.savefig(fig); plt.close(fig)

#         # chi2/ndf — now a proper statistical quantity (Route 1)
#         if "chi2ndf" in df_vtx.columns:
#             for sub_df, lbl in [(df_vtx, "all N"), (df_vtx[df_vtx["n_tracks"]>=3], "N≥3")]:
#                 chi_all = sub_df["chi2ndf"].replace([np.inf,-np.inf],np.nan).dropna()
#                 chi_acc = df_acc.loc[df_acc.index.isin(sub_df.index), "chi2ndf"]\
#                                 .replace([np.inf,-np.inf],np.nan).dropna()
#                 if len(chi_all) == 0: continue
#                 for xlim, sfx in [(None,"full range"),((0,4),"0–4")]:
#                     fig, ax = plt.subplots(figsize=(6,4))
#                     d_all = chi_all[(chi_all>=xlim[0])&(chi_all<=xlim[1])] if xlim else chi_all
#                     d_acc = chi_acc[(chi_acc>=xlim[0])&(chi_acc<=xlim[1])] if xlim else chi_acc
#                     ax.hist(d_all, bins=50, histtype="step", label="All",      density=True)
#                     if len(d_acc):
#                         ax.hist(d_acc, bins=50, histtype="step", ls="--",
#                                 label="Accepted", density=True)
#                     ax.axvline(1.0, color="red", ls=":", label="χ²/ndf = 1")
#                     if xlim: ax.set_xlim(*xlim)
#                     ax.set_xlabel("χ²/ndf (hit-space, Route 1)")
#                     ax.set_ylabel("Norm. counts")
#                     ax.set_title(f"Vertex χ²/ndf — {lbl} ({sfx})")
#                     ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
#                     pdf.savefig(fig); plt.close(fig)

#         if "n_tracks" in df_vtx.columns:
#             bins = range(int(df_vtx["n_tracks"].min()), int(df_vtx["n_tracks"].max())+2)
#             fig, ax = plt.subplots(figsize=(6,4))
#             ax.hist(df_vtx["n_tracks"], bins=bins, histtype="step", label="All")
#             if len(df_acc):
#                 ax.hist(df_acc["n_tracks"], bins=bins, histtype="step", ls="--", label="Accepted")
#             ax.set_xlabel("N tracks"); ax.set_ylabel("Counts")
#             ax.set_title("Vertex multiplicity")
#             ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
#             pdf.savefig(fig); plt.close(fig)

#         if "vtx_d3_true" in df_vtx.columns:
#             d3_all = df_vtx["vtx_d3_true"].replace([np.inf,-np.inf],np.nan).dropna()
#             d3_acc = df_acc["vtx_d3_true"].replace([np.inf,-np.inf],np.nan).dropna()
#             fig, ax = plt.subplots(figsize=(6,4))
#             if len(d3_all): ax.hist(d3_all, bins=50, histtype="step", label="All",      density=True)
#             if len(d3_acc): ax.hist(d3_acc, bins=50, histtype="step", ls="--", label="Accepted", density=True)
#             ax.set_xlabel("|v_reco - v_truth| [cm]"); ax.set_ylabel("Norm. counts")
#             ax.set_title("Vertex distance to truth")
#             ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
#             pdf.savefig(fig); plt.close(fig)

#         if df_res.empty:
#             print("[control plots] residuals tree empty — skipping track plots.")
#             return

#         if "is_signal" in df_res.columns:
#             df_res_merge = df_res.merge(
#                 df_vtx[["event_id","vertex_index","accepted"]],
#                 on=["event_id","vertex_index"], how="left"
#             )
#             df_acc_trk = df_res_merge[df_res_merge["accepted"]==1]
#             tot_sig = int((df_acc_trk["is_signal"]==1).sum())
#             tot_bg  = int((df_acc_trk["is_signal"]==0).sum())
#             tot = tot_sig + tot_bg
#             if tot > 0:
#                 fig, ax = plt.subplots(figsize=(5,4))
#                 ax.bar([0,1],[tot_sig/tot, tot_bg/tot])
#                 ax.set_xticks([0,1]); ax.set_xticklabels(["Signal","Background"])
#                 ax.set_ylim(0,1.05); ax.set_ylabel("Fraction of tracks in accepted vtx")
#                 ax.set_title("Signal vs background content")
#                 for x, y in zip([0,1],[tot_sig/tot, tot_bg/tot]):
#                     ax.text(x, y+0.02, f"{y:.2f}", ha="center", fontsize=8)
#                 ax.grid(True, axis="y", alpha=0.3); fig.tight_layout()
#                 pdf.savefig(fig); plt.close(fig)

#         needed = {"delta_x","delta_y","delta_z","dca","sigma_x","sigma_y","sigma_z",
#                   "u_x","u_y","u_z","dca_truth","delta_truth_x","delta_truth_y","delta_truth_z"}
#         if not needed.issubset(df_res.columns):
#             print("[control plots] Missing columns:", needed - set(df_res.columns))
#             return

#         df_res = df_res.copy()
#         df_res["theta_deg"] = _compute_theta_deg(
#             df_res["u_x"].to_numpy(), df_res["u_y"].to_numpy(), df_res["u_z"].to_numpy()
#         )
#         df_res["dca_residual"] = df_res["dca_truth"] - df_res["dca"]

#         fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
#         ax1.hist(df_res["dca_truth"].dropna(), bins=60, histtype="step", label="before (true)", density=True)
#         ax1.hist(df_res["dca"].dropna(),        bins=60, histtype="step", ls="--", label="after (fit)",  density=True)
#         ax1.set_xlabel("DCA [cm]"); ax1.set_ylabel("Norm. counts")
#         ax1.set_title("DCA before/after fit"); ax1.legend(); ax1.grid(True, alpha=0.3)
#         ax2.hist(df_res["dca_residual"].dropna(), bins=60, histtype="step", density=True)
#         ax2.set_xlabel("DCA_before - DCA_after [cm]"); ax2.set_ylabel("Norm. counts")
#         ax2.set_title("DCA residual"); ax2.grid(True, alpha=0.3)
#         fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

#         for lo, hi in THETA_RANGES_DEG:
#             sub = df_res[(df_res["theta_deg"]>=lo)&(df_res["theta_deg"]<hi)]
#             if sub.empty: continue
#             fig, axes = plt.subplots(1,3,figsize=(12,4))
#             for ax, col in zip(axes,["sigma_x","sigma_y","sigma_z"]):
#                 ax.hist(sub[col].dropna(), bins=50, histtype="step", density=True)
#                 ax.set_xlabel(f"{col} [cm]"); ax.set_ylabel("Norm. counts")
#                 ax.set_title(f"{col} theta [{lo:.0f}–{hi:.0f}]"); ax.grid(True,alpha=0.3)
#             fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

#         for lo, hi in THETA_RANGES_DEG:
#             sub = df_res[(df_res["theta_deg"]>=lo)&(df_res["theta_deg"]<hi)]
#             if sub.empty: continue
#             fig, ax = plt.subplots(figsize=(6,4))
#             ax.hist(sub["dca_residual"].dropna(), bins=60, histtype="step", density=True)
#             ax.set_xlabel("DELTA DCA [cm]"); ax.set_ylabel("Norm. counts")
#             ax.set_title(f"DELTA DCA theta [{lo:.0f}–{hi:.0f}]")
#             ax.grid(True,alpha=0.3); fig.tight_layout()
#             pdf.savefig(fig); plt.close(fig)

#         have_plane = "dca_plane_truth" in df_res.columns
#         theta_mids, rms_dca, rms_plane = [], [], []
#         for lo, hi in THETA_RANGES_DEG:
#             sub = df_res[(df_res["theta_deg"]>=lo)&(df_res["theta_deg"]<hi)]
#             if sub.empty: continue
#             vals       = sub["dca_truth"].dropna().to_numpy()
#             vals_plane = sub["dca_plane_truth"].dropna().to_numpy() if have_plane else np.array([])
#             if vals.size == 0 and vals_plane.size == 0: continue
#             fig, ax = plt.subplots(figsize=(6,4))
#             if vals.size:
#                 ax.hist(vals,       bins=60, histtype="step", density=True, label="minimal DCA_T(true)")
#             if vals_plane.size:
#                 ax.hist(vals_plane, bins=60, histtype="step", density=True, ls="--",
#                         label="plane-projected DCA_T(true)")
#             ax.set_xlabel("DCA_T to true vtx [cm]"); ax.set_ylabel("Norm. counts")
#             ax.set_title(f"DCA_T(true) theta [{lo:.0f}–{hi:.0f}]")
#             ax.legend(); ax.grid(True,alpha=0.3); fig.tight_layout()
#             pdf.savefig(fig); plt.close(fig)
#             theta_mids.append(0.5*(lo+hi))
#             rms_dca.append(float(np.sqrt(np.mean(vals**2)))        if vals.size       else float("nan"))
#             rms_plane.append(float(np.sqrt(np.mean(vals_plane**2))) if vals_plane.size else float("nan"))

#         if theta_mids:
#             fig, ax = plt.subplots(figsize=(6,4))
#             ax.plot(theta_mids, rms_dca,   marker="o", label="minimal DCA_T(true)")
#             if any(np.isfinite(rms_plane)):
#                 ax.plot(theta_mids, rms_plane, marker="s", label="plane-projected DCA_T(true)")
#             ax.set_xlabel("Theta [deg]"); ax.set_ylabel("RMS(DCA_T to true vtx) [cm]")
#             ax.set_title("RMS transverse DCA to true vtx vs theta")
#             ax.legend(); ax.grid(True,alpha=0.3); fig.tight_layout()
#             pdf.savefig(fig); plt.close(fig)

#         if {"vtx_x","vtx_y","vtx_z"}.issubset(df_vtx.columns):
#             df_acc_merged = df_res.merge(
#                 df_vtx.loc[df_vtx["accepted"]==1,
#                            ["event_id","vertex_index","vtx_x","vtx_y","vtx_z"]],
#                 on=["event_id","vertex_index"], how="inner"
#             )
#             for theta_spec in THETA_RANGES_DEG + [("all","all")]:
#                 if theta_spec == ("all","all"):
#                     sub = df_acc_merged; lbl = "all theta"
#                 else:
#                     lo, hi = theta_spec
#                     sub = df_acc_merged[(df_acc_merged["theta_deg"]>=lo)&(df_acc_merged["theta_deg"]<hi)]
#                     lbl = f"{lo:.0f}–{hi:.0f} deg"
#                 if sub.empty: continue
#                 fig, axes = plt.subplots(1,3,figsize=(12,4))
#                 for ax, comp in zip(axes,["x","y","z"]):
#                     ax.hist(sub[f"delta_{comp}"].dropna(),       bins=60, histtype="step",
#                             label="track - fitted vtx",  density=True)
#                     ax.hist(sub[f"delta_truth_{comp}"].dropna(), bins=60, histtype="step",
#                             ls="--", label="track - true vtx", density=True)
#                     ax.set_xlabel(f"Delta {comp} [cm]"); ax.set_ylabel("Norm. counts")
#                     ax.set_title(f"{comp}-residuals, {lbl}"); ax.legend(); ax.grid(True,alpha=0.3)
#                 fig.suptitle(f"Track residuals ({lbl})", y=1.03)
#                 fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

#     print(f"[control plots] Wrote {output_pdf}")


# # -----------------------------------------------------------------------
# # Main
# # -----------------------------------------------------------------------

# def main():
#     parser = argparse.ArgumentParser(
#         description="HIBEAM Route-1 vertex finder: hit-space chi2 with proper ndf."
#     )
#     parser.add_argument("input_root")
#     parser.add_argument("--tree",         dest="tree_name",    default="digitizedHits")
#     parser.add_argument("--output-root",  dest="output_root",  default="vertex_products.root")
#     parser.add_argument("--output-pdf",   dest="output_pdf",   default="vertex_diagnostics.pdf")

#     parser.add_argument("--tpc-sigma-x",  dest="tpc_sigma_x",  type=float, default=0.05)
#     parser.add_argument("--tpc-sigma-y",  dest="tpc_sigma_y",  type=float, default=0.05)
#     parser.add_argument("--tpc-sigma-z",  dest="tpc_sigma_z",  type=float, default=0.10)

#     parser.add_argument("--chi2ndf-max",  dest="chi2ndf_max",  type=float, default=3.0)
#     parser.add_argument("--max-dca",      dest="max_dca",      type=float, default=3.0)
#     parser.add_argument("--max-r-vertex", dest="max_r_vertex", type=float, default=20.0)

#     parser.add_argument("--min-comb",     dest="min_comb",     type=int,   default=3)
#     parser.add_argument("--max-comb",     dest="max_comb",     type=int,   default=5)
#     parser.add_argument("--z-constraint", dest="z_constraint", type=float, default=0.0)

#     parser.add_argument("--use-tpc-smearing", dest="use_tpc_smearing",
#                         action="store_true", default=False)
#     parser.add_argument("--use-ms-sigma",     dest="use_ms_sigma",
#                         action="store_true", default=True)

#     parser.add_argument("--ms-density",      dest="ms_density",      type=float, default=2.7)
#     parser.add_argument("--ms-radlen-gcm2",  dest="ms_radlen_gcm2",  type=float, default=24.0)
#     parser.add_argument("--ms-thickness-cm", dest="ms_thickness_cm", type=float, default=0.5)
#     parser.add_argument("--ms-lever-arm-cm", dest="ms_lever_arm_cm", type=float, default=20.0)
#     parser.add_argument("--ms-p-MeV",        dest="ms_p_MeV",        type=float, default=200.0)
#     parser.add_argument("--ms-min-cos",       dest="ms_min_cos",      type=float, default=0.05)

#     parser.add_argument("--seed",       type=int, default=None)
#     parser.add_argument("--card",       type=str, default=None)
#     parser.add_argument("--max-events", dest="max_events", type=int, default=0)

#     args = parser.parse_args()

#     if args.card:
#         apply_card_overrides(args, load_steering_card(args.card))
#     if args.seed is not None:
#         np.random.seed(args.seed)

#     comb_sizes = list(range(args.min_comb, args.max_comb + 1))
#     print("[config] Route-1 hit-space chi2 enabled")
#     print("[config] comb_sizes =", comb_sizes)
#     print(f"[config] chi2ndf_max={args.chi2ndf_max}  max_dca={args.max_dca}"
#           f"  max_r_vertex={args.max_r_vertex}")
#     print(f"[config] tpc_sigma = ({args.tpc_sigma_x}, {args.tpc_sigma_y},"
#           f" {args.tpc_sigma_z}) cm")
#     print(f"[config] ndf formula: 2N - 2  (6N obs - (2 + 4N) params)")

#     total_events = 0
#     found_counts = {n: 0 for n in comb_sizes}
#     total_tracks = idx_p = idx_q = idx_r = idx_s = idx_t = 0
#     event_diag = {
#         n: dict(too_few_tracks=0, geN_no_accepted=0,
#                 fail_chi2=0, fail_dca=0, fail_r=0, other_fail=0)
#         for n in comb_sizes
#     }

#     writer = RootWriter(args.output_root)

#     for event_id, df_event in iter_events_from_root(args.input_root, args.tree_name):
#         total_events += 1
#         if args.max_events > 0 and total_events > args.max_events:
#             break

#         ev_result = process_event(
#             df_event=df_event,
#             comb_sizes=comb_sizes,
#             tpc_sigma_x=args.tpc_sigma_x, tpc_sigma_y=args.tpc_sigma_y,
#             tpc_sigma_z=args.tpc_sigma_z,
#             chi2ndf_max=args.chi2ndf_max,  max_dca_cut=args.max_dca,
#             max_r_vertex=args.max_r_vertex, z_constraint=args.z_constraint,
#             use_tpc_smearing=args.use_tpc_smearing, use_ms_sigma=args.use_ms_sigma,
#             ms_density=args.ms_density, ms_radlen_gcm2=args.ms_radlen_gcm2,
#             ms_thickness_cm=args.ms_thickness_cm, ms_lever_arm_cm=args.ms_lever_arm_cm,
#             ms_p_MeV=args.ms_p_MeV, ms_min_cos=args.ms_min_cos,
#         )

#         vtx_rows: List[Dict] = []
#         res_rows: List[Dict] = []

#         for vtx_index, v in enumerate(ev_result["vertices"]):
#             vtx_rows.append({
#                 "event_id": event_id, "vertex_index": vtx_index,
#                 "n_tracks": v["n_tracks"],
#                 "vtx_x": v["vtx_x"], "vtx_y": v["vtx_y"], "vtx_z": v["vtx_z"],
#                 "vtx_r": v["vtx_r"],
#                 "vtx_dx_true": v["vtx_dx_true"], "vtx_dy_true": v["vtx_dy_true"],
#                 "vtx_dz_true": v["vtx_dz_true"], "vtx_d3_true": v["vtx_d3_true"],
#                 "chi2": v["chi2"], "ndf": v["ndf"], "chi2ndf": v["chi2ndf"],
#                 "cond": v["cond"],
#                 "max_dca": v["max_dca"], "avg_dca": v["avg_dca"],
#                 "accepted_raw": int(v["accepted_raw"]), "accepted": int(v["accepted"]),
#                 "chi2ndf_max_used":      v["chi2ndf_max_used"],
#                 "max_dca_cut_used":      v["max_dca_cut_used"],
#                 "max_r_vertex_cut_used": v["max_r_vertex_cut_used"],
#                 "fail_chi2": int(v["fail_chi2"]),
#                 "fail_dca":  int(v["fail_dca"]),
#                 "fail_r":    int(v["fail_r"]),
#             })
#             for res in v["residuals"]:
#                 res_rows.append({
#                     "event_id": event_id, "vertex_index": vtx_index,
#                     "n_tracks": v["n_tracks"],
#                     "track_id":  res.get("track_id", -1),
#                     "is_signal": res.get("is_signal", 1),
#                     "delta_x": res["delta_x"], "delta_y": res["delta_y"],
#                     "delta_z": res["delta_z"],
#                     "dca":             res["dca"],
#                     "dca_truth":       res.get("dca_truth",       float("nan")),
#                     "dca_truth_3d":    res.get("dca_truth_3d",    float("nan")),
#                     "dca_plane_truth": res.get("dca_plane_truth", float("nan")),
#                     "delta_truth_x":   res.get("delta_truth_x",   float("nan")),
#                     "delta_truth_y":   res.get("delta_truth_y",   float("nan")),
#                     "delta_truth_z":   res.get("delta_truth_z",   float("nan")),
#                     "sigma_x": res.get("sigma_x", float("nan")),
#                     "sigma_y": res.get("sigma_y", float("nan")),
#                     "sigma_z": res.get("sigma_z", float("nan")),
#                     "p_x": res.get("p_x", float("nan")),
#                     "p_y": res.get("p_y", float("nan")),
#                     "p_z": res.get("p_z", float("nan")),
#                     "u_x": res.get("u_x", float("nan")),
#                     "u_y": res.get("u_y", float("nan")),
#                     "u_z": res.get("u_z", float("nan")),
#                 })

#         writer.flush(vtx_rows, res_rows)

#         n_valid    = ev_result["n_valid_tracks"]
#         diag_by_N  = ev_result["diag_by_N"]
#         for n in comb_sizes:
#             if n_valid < n:
#                 event_diag[n]["too_few_tracks"] += 1
#             else:
#                 if not ev_result["found_by_N"].get(n, False):
#                     event_diag[n]["geN_no_accepted"] += 1
#                     dN = diag_by_N.get(n, {})
#                     if dN.get("any_chi2_fail"): event_diag[n]["fail_chi2"] += 1
#                     if dN.get("any_dca_fail"):  event_diag[n]["fail_dca"]  += 1
#                     if dN.get("any_r_fail"):    event_diag[n]["fail_r"]    += 1
#                     if dN.get("any_combo_considered") and not (
#                         dN.get("any_chi2_fail") or dN.get("any_dca_fail") or dN.get("any_r_fail")
#                     ):
#                         event_diag[n]["other_fail"] += 1
#             if ev_result["found_by_N"].get(n, False):
#                 found_counts[n] += 1

#         for _, meta in ev_result["track_meta"].items():
#             total_tracks += 1
#             if meta["invalid_exit"] and not meta["used_in_any_accepted"]:     idx_p += 1
#             if (not meta["invalid_exit"]) and (not meta["used_in_any_comb"]) \
#                     and (not meta["used_in_any_accepted"]):                    idx_q += 1
#             if meta["used_in_any_comb"] and not meta["used_in_any_accepted"]:
#                 if meta["in_combo_fail_chi2"]: idx_r += 1
#                 if meta["in_combo_fail_dca"]:  idx_s += 1
#                 if meta["in_combo_fail_r"]:    idx_t += 1

#         if total_events % 500 == 0:
#             print(f"  ... processed {total_events} events")

#     writer.close()

#     print(f"\n=== Summary: {total_events} events processed ===")
#     for n in comb_sizes:
#         frac = found_counts[n] / total_events if total_events else 0.0
#         print(f"  Accepted {n}-track vertices: {found_counts[n]}/{total_events} ({frac:.3f})")
#     print(f"\nTrack indices: total={total_tracks}, p={idx_p}, q={idx_q}, "
#           f"r={idx_r}, s={idx_s}, t={idx_t}")
#     for n in comb_sizes:
#         d = event_diag[n]
#         print(f"\nN={n}: too_few={d['too_few_tracks']}, no_accepted={d['geN_no_accepted']} "
#               f"[chi2={d['fail_chi2']}, dca={d['fail_dca']}, r={d['fail_r']}, "
#               f"other={d['other_fail']}]")

#     try:
#         make_control_plots(args.output_root, args.output_pdf)
#     except Exception as exc:
#         print(f"[control plots] Warning: {exc}")
#         import traceback; traceback.print_exc()


# if __name__ == "__main__":
#     main()