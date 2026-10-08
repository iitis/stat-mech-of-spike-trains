#!/usr/bin/env python3
"""
Analyze phase diagram results saved by run_phase_diagram.m

Now supports merging multiple .mat files with identical format, split across beta ranges.

Usage:
    python plotPhaseMulti.py Gphase30x10_beta0to2_h12_phase.mat Gphase30x11_beta2to4_h12_phase.mat --outdir results_h12
"""

from __future__ import annotations
import argparse
import os
from dataclasses import dataclass
from typing import Optional, Dict, Tuple, List
from matplotlib.colors import ListedColormap, BoundaryNorm

import numpy as np
import matplotlib.pyplot as plt

def fmt_sci_math(x: float) -> str:
    """
    Format number as mathtext: 2×10^{5}, 10^{5}, or plain for small exponents.
    """
    if x == 0 or not np.isfinite(x):
        return "0"
    ax = abs(x)
    exp = int(np.floor(np.log10(ax)))
    mant = x / (10 ** exp)

    # dla "normalnych" wartości nie róbmy cyrku z wykładnikiem
    if -2 <= exp <= 3:
        return f"{x:g}"

    # mantysa blisko 1 -> samo 10^{exp}
    if np.isclose(abs(mant), 1.0, rtol=1e-6, atol=1e-12):
        sign = "-" if x < 0 else ""
        return rf"${sign}10^{{{exp}}}$"

    # standard: a×10^{b}
    return rf"${mant:.2g}\times 10^{{{exp}}}$"



def load_mat_any(path: str) -> dict:
    """
    Loads MATLAB .mat that may be v7.3 (HDF5) or older.
    Returns a dict with numpy arrays.
    """
    # Try scipy.io.loadmat first (works for v7.2 and earlier)
    try:
        import scipy.io as sio
        d = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
        d = {k: v for k, v in d.items() if not k.startswith("__")}
        return d
    except Exception:
        pass

    # Fall back to h5py for v7.3
    try:
        import h5py
        out = {}
        with h5py.File(path, "r") as f:
            for k in f.keys():
                obj = f[k]
                if isinstance(obj, h5py.Dataset):
                    out[k] = obj[()]

        # Heuristic transpose fixes
        if "betas" in out:
            out["betas"] = np.array(out["betas"]).reshape(-1)
        if "sdels" in out:
            out["sdels"] = np.array(out["sdels"]).reshape(-1)

        if "x_grid" in out:
            xg = np.array(out["x_grid"])
            nb = int(np.asarray(out["betas"]).reshape(-1).shape[0])
            ns = int(np.asarray(out["sdels"]).reshape(-1).shape[0])

            T_file = None
            if "T" in out:
                try:
                    T_file = int(np.asarray(out["T"]).squeeze())
                except Exception:
                    T_file = None

            if xg.ndim != 3:
                raise ValueError(f"x_grid ma złe ndim={xg.ndim}, shape={xg.shape}")

            if T_file is not None and xg.shape[0] == T_file and xg.shape[1] == ns and xg.shape[2] == nb:
                out["x_grid"] = np.transpose(xg, (2, 1, 0))
            elif T_file is not None and xg.shape[0] == T_file and xg.shape[1] == nb and xg.shape[2] == ns:
                out["x_grid"] = np.transpose(xg, (1, 2, 0))
            else:
                candidates = [
                    xg,
                    np.transpose(xg, (2, 1, 0)),
                    np.transpose(xg, (1, 2, 0)),
                    np.transpose(xg, (2, 0, 1)),
                    np.transpose(xg, (0, 2, 1)),
                    np.transpose(xg, (1, 0, 2)),
                ]
                fixed = None
                for c in candidates:
                    if c.ndim == 3 and c.shape[0] == nb and c.shape[1] == ns:
                        fixed = c
                        break
                if fixed is None:
                    raise ValueError(
                        f"Nie mogę dopasować wymiarów x_grid. Odczytano shape={xg.shape}, "
                        f"spodziewam się (nb={nb}, ns={ns}, T)."
                    )
                out["x_grid"] = fixed

        if "E_final" in out:
            Ef = np.array(out["E_final"])
            nb = out["betas"].shape[0]
            ns = out["sdels"].shape[0]
            if Ef.ndim == 2 and Ef.shape != (nb, ns):
                if Ef.T.shape == (nb, ns):
                    Ef = Ef.T
            out["E_final"] = Ef

        return out
    except Exception as e:
        raise RuntimeError(
            f"Nie udało się wczytać pliku .mat ani przez scipy, ani przez h5py.\n"
            f"Plik: {path}\nBłąd: {e}"
        )


def normalized_autocorr_peak(x: np.ndarray, lag_max: Optional[int] = None) -> float:
    x = np.asarray(x, dtype=np.float64)
    T = x.size
    if T < 3:
        return 0.0

    if lag_max is None:
        lag_max = T // 2
    lag_max = int(max(1, min(lag_max, T - 1)))

    x = x - x.mean()
    var = np.dot(x, x)
    if var <= 1e-12:
        return 0.0

    n = 1 << (2 * T - 1).bit_length()
    fx = np.fft.rfft(x, n=n)
    ac = np.fft.irfft(fx * np.conjugate(fx), n=n)[:T]
    ac /= ac[0] if ac[0] != 0 else 1.0

    peak = float(np.max(ac[1:lag_max + 1]))
    if not np.isfinite(peak):
        peak = 0.0
    return max(0.0, min(1.0, peak))


def ipr_concentration(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    s = x.sum()
    if s <= 0:
        return 0.0
    p = x / s
    return float(np.sum(p * p))

def temporal_entropy(x: np.ndarray, normalize: bool = True) -> float:
    """
    Entropy of normalized population activity x(t).

    H = -sum_t p_t log(p_t), where p_t = x_t / sum_t x_t.
    If normalize=True, returns H/log(T), so the range is approximately [0,1].
    """
    x = np.asarray(x, dtype=np.float64)
    T = x.size
    s = x.sum()

    if s <= 0 or T <= 1:
        return np.nan

    p = x / s
    p = p[p > 0]

    H = -float(np.sum(p * np.log(p)))

    if normalize:
        H = H / np.log(T)

    return H

@dataclass
class PhaseParams:
    osc_thresh: float = 0.25
    ipr_thresh: float = 0.01
    lag_max_frac: float = 0.5


def classify_phase(op_osc: float, op_ipr: float, params: PhaseParams) -> int:
    if op_ipr > params.ipr_thresh:
        return 2
    if op_osc > params.osc_thresh:
        return 1
    return 0


def _choose_beta_ticks(betas: np.ndarray, max_ticks: int = 12) -> Tuple[np.ndarray, List[str]]:
    """
    Returns tick positions (indices) and labels for beta axis.
    - fewer ticks (<= max_ticks)
    - strong digit truncation via format
    """
    betas = np.asarray(betas, dtype=float).reshape(-1)
    nb = betas.size
    if nb == 0:
        return np.array([], dtype=int), []

    if nb <= max_ticks:
        idx = np.arange(nb, dtype=int)
    else:
        idx = np.unique(np.round(np.linspace(0, nb - 1, max_ticks)).astype(int))

    # "mocne obcięcie": 3 significant digits usually enough
    labels = [fmt_sci_math(betas[i]) for i in idx]
    return idx, labels


def heatmap(betas, sdels, Z, title, outpath, is_discrete=False, max_beta_ticks: int = 6):
    betas = np.asarray(betas, dtype=float).reshape(-1)
    sdels = np.asarray(sdels, dtype=float).reshape(-1)
    Z = np.asarray(Z)

    # Z is (nb, ns) -> plot as (ns, nb)
    C = Z.T

    plt.figure()
    if is_discrete:
    # 3 stałe kolory dla 0/1/2
        cmap = ListedColormap(["#1f77b4", "#8c564b", "#17becf"])  # możesz zmienić
        norm = BoundaryNorm(boundaries=[-0.5, 0.5, 1.5, 2.5], ncolors=cmap.N)

        im = plt.imshow(C, origin="lower", aspect="auto", cmap=cmap, norm=norm)
        cb = plt.colorbar(im, ticks=[0, 1, 2])
        cb.ax.set_yticklabels(["0", "1", "2"])
        cb.set_label("phase (0=uniform, 1=osc, 2=collapsed)")
    else:
        im = plt.imshow(C, origin="lower", aspect="auto")
        cb = plt.colorbar(im)
        cb.set_label("value")

    # beta ticks: fewer + truncated labels
    xt_idx, xt_lbl = _choose_beta_ticks(betas, max_ticks=max_beta_ticks)
    plt.xticks(xt_idx, xt_lbl, rotation=0)
    ax = plt.gca()
    ax.xaxis.get_offset_text().set_visible(False)
    plt.xticks(rotation=30, ha="right")

    plt.yticks(np.arange(sdels.size), [f"{s:g}" for s in sdels])

    plt.xlabel("beta")
    plt.ylabel("delay std (sdel)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def merge_phase_mats(paths: List[str]) -> Dict[str, np.ndarray]:
    """
    Merge multiple phase .mat files that share identical sdels and T,
    but have different (possibly disjoint) beta ranges.

    Returns dict with: betas_all, sdels, x_grid_all, (optional) E_final_all
    """
    if len(paths) == 0:
        raise ValueError("Podaj przynajmniej jeden plik .mat.")

    parts = []
    sdels_ref = None
    T_ref = None
    ns_ref = None

    for p in paths:
        d = load_mat_any(p)
        if "betas" not in d or "sdels" not in d or "x_grid" not in d:
            raise KeyError(f"Plik {p} musi zawierać betas, sdels, x_grid.")

        betas = np.asarray(d["betas"], dtype=np.float64).reshape(-1)
        sdels = np.asarray(d["sdels"], dtype=np.float64).reshape(-1)
        x_grid = np.asarray(d["x_grid"])
        if x_grid.ndim != 3:
            raise ValueError(f"{p}: x_grid powinno mieć (nb,ns,T), ma {x_grid.shape}")

        nb, ns, T = x_grid.shape
        if nb != betas.size or ns != sdels.size:
            raise ValueError(f"{p}: Niezgodne wymiary: x_grid={x_grid.shape}, betas={betas.size}, sdels={sdels.size}")

        if sdels_ref is None:
            sdels_ref = sdels
            ns_ref = ns
            T_ref = T
        else:
            if ns != ns_ref or T != T_ref:
                raise ValueError(f"{p}: ns/T niezgodne z referencją: (ns,T)=({ns},{T}) vs ({ns_ref},{T_ref})")
            if not np.allclose(sdels, sdels_ref, rtol=0, atol=0):
                raise ValueError(f"{p}: sdels różne od referencji. Nie scalę tego w sensowny sposób.")

        Ef = None
        if "E_final" in d:
            Ef = np.asarray(d["E_final"])
            if Ef.ndim != 2 or Ef.shape != (nb, ns):
                raise ValueError(f"{p}: E_final ma złe wymiary: {Ef.shape}, oczekuję ({nb},{ns})")

        parts.append({"path": p, "betas": betas, "x_grid": x_grid, "E_final": Ef})

    # Build global beta axis
    betas_all = np.unique(np.concatenate([pt["betas"] for pt in parts]))
    betas_all.sort()
    nb_all = betas_all.size
    ns = int(ns_ref)
    T = int(T_ref)

    # Allocate merged arrays
    x_grid_all = np.full((nb_all, ns, T), np.nan, dtype=np.float64)
    have_E = any(pt["E_final"] is not None for pt in parts)
    E_final_all = np.full((nb_all, ns), np.nan, dtype=np.float64) if have_E else None

    filled = np.zeros((nb_all,), dtype=bool)

    # Fill from parts
    for pt in parts:
        b = pt["betas"]
        x = pt["x_grid"]

        # map each beta in this file to global index
        # exact match assumed because betas are generated deterministically
        idx = np.searchsorted(betas_all, b)
        if not np.all(betas_all[idx] == b):
            raise ValueError(f"{pt['path']}: beta values do not match global grid exactly (floating mismatch).")

        # handle overlaps: keep first, warn if conflict
        for k, gi in enumerate(idx):
            if filled[gi]:
                # overlap: check if identical-ish
                prev = x_grid_all[gi, :, :]
                cur = x[k, :, :].astype(np.float64)
                if not np.allclose(prev, cur, equal_nan=True):
                    print(f"[WARN] Overlap beta={betas_all[gi]:.6g}: keeping first, {pt['path']} differs.")
                continue
            x_grid_all[gi, :, :] = x[k, :, :].astype(np.float64)
            filled[gi] = True

            if have_E and pt["E_final"] is not None:
                E_final_all[gi, :] = pt["E_final"][k, :].astype(np.float64)

    if not np.all(filled):
        missing = np.where(~filled)[0]
        miss_vals = betas_all[missing]
        raise RuntimeError(f"Po scaleniu brakuje danych dla {missing.size} wartości beta, np. {miss_vals[:10]} ...")

    # Optionally cast back to original-ish (not necessary for analysis)
    return {
        "betas": betas_all,
        "sdels": sdels_ref,
        "x_grid": x_grid_all,
        "E_final": E_final_all,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matfiles", nargs="+", help="One or more *_phase.mat produced by run_phase_diagram (split by beta ranges)")
    ap.add_argument("--outdir", default="phase_results", help="Output directory for plots and npz")
    ap.add_argument("--osc-thresh", type=float, default=0.25, help="Threshold for oscillation order parameter")
    ap.add_argument("--ipr-thresh", type=float, default=0.01, help="Threshold for IPR (collapse) order parameter")
    ap.add_argument("--lag-max-frac", type=float, default=0.5, help="Max lag fraction of T for autocorr peak (0..1]")
    ap.add_argument("--auto-thresholds", action="store_true",
                    help="Set thresholds from data quantiles (useful for quick first pass).")
    ap.add_argument("--max-beta-ticks", type=int, default=12,
                    help="Max number of tick labels on beta axis (avoid overlapping labels).")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    merged = merge_phase_mats(args.matfiles)
    betas = np.asarray(merged["betas"], dtype=np.float64).reshape(-1)
    sdels = np.asarray(merged["sdels"], dtype=np.float64).reshape(-1)
    x_grid = np.asarray(merged["x_grid"], dtype=np.float64)

    nb, ns, T = x_grid.shape

    params = PhaseParams(
        osc_thresh=float(args.osc_thresh),
        ipr_thresh=float(args.ipr_thresh),
        lag_max_frac=float(args.lag_max_frac),
    )
    lag_max = int(max(1, min(T - 1, np.floor(params.lag_max_frac * T))))

    op1_osc = np.zeros((nb, ns), dtype=np.float64)
    op2_ipr = np.zeros((nb, ns), dtype=np.float64)
    op3_entropy = np.zeros((nb, ns), dtype=np.float64)

    for ib in range(nb):
        for is_ in range(ns):
            x = x_grid[ib, is_].astype(np.float64)
            op1_osc[ib, is_] = normalized_autocorr_peak(x, lag_max=lag_max)
            op2_ipr[ib, is_] = ipr_concentration(x)
            op3_entropy[ib, is_] = temporal_entropy(x, normalize=True)

    if args.auto_thresholds:
        params.osc_thresh = float(np.quantile(op1_osc, 0.75))
        params.ipr_thresh = float(np.quantile(op2_ipr, 0.75))
        print(f"[auto] osc_thresh={params.osc_thresh:.4g}, ipr_thresh={params.ipr_thresh:.4g}")

    phase = np.zeros((nb, ns), dtype=np.int8)
    for ib in range(nb):
        for is_ in range(ns):
            phase[ib, is_] = classify_phase(op1_osc[ib, is_], op2_ipr[ib, is_], params)

    np.savez_compressed(
        os.path.join(args.outdir, "phase_analysis.npz"),
        betas=betas,
        sdels=sdels,
        op1_osc=op1_osc,
        op2_ipr=op2_ipr,
        op3_entropy=op3_entropy,
        phase=phase,
        osc_thresh=params.osc_thresh,
        ipr_thresh=params.ipr_thresh,
        lag_max=lag_max,
        source_mats=[os.path.abspath(p) for p in args.matfiles],
    )

    heatmap(
        betas, sdels, op1_osc,
        title="OP1: oscillation strength (max normalized autocorr peak)",
        outpath=os.path.join(args.outdir, "heatmap_op1_osc.png"),
        is_discrete=False,
        max_beta_ticks=args.max_beta_ticks,
    )
    heatmap(
        betas, sdels, op2_ipr,
        title="OP2: concentration (IPR of normalized x(t))",
        outpath=os.path.join(args.outdir, "heatmap_op2_ipr.png"),
        is_discrete=False,
        max_beta_ticks=args.max_beta_ticks,
    )

    heatmap(
        betas, sdels, op3_entropy,
        title="OP3: normalized temporal entropy H/log(T)",
        outpath=os.path.join(args.outdir, "heatmap_op3_entropy.png"),
        is_discrete=False,
        max_beta_ticks=args.max_beta_ticks,
    )

    heatmap(
        betas, sdels, phase,
        title="Phase map (0=uniform, 1=oscillatory, 2=collapsed)",
        outpath=os.path.join(args.outdir, "heatmap_phase.png"),
        is_discrete=True,
        max_beta_ticks=args.max_beta_ticks,
    )

    print(f"Zapisano wyniki do: {args.outdir}")
    print("Pliki:")
    print("- heatmap_op1_osc.png")
    print("- heatmap_op2_ipr.png")
    print("- heatmap_op3_entropy.png")
    print("- heatmap_phase.png")
    print("- phase_analysis.npz")
    

if __name__ == "__main__":
    main()

