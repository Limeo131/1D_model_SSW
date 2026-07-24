#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gray-box calibration: optimize (beta, s_tau) to match ERA5 10-hPa seasonal cycle.
- Calls Fortran model via run.sh (foreground, auto-compiles)
- Reads uout_real_urad_<tag>.nc output (time, height)
- Computes RMSE of 10-hPa multi-year seasonal cycle vs ERA5
- Uses coarse grid then Powell optimization













  python optimize_beta_stau.py --exe ./d.out --runner direct
"""

import os, math, json, argparse, subprocess
import numpy as np
from pathlib import Path
from netCDF4 import Dataset
from scipy.optimize import minimize

# -------------------- Default parameters --------------------
DEFAULT_EXE = "./run.sh"       # can also use ./d.out directly
DEFAULT_RUNNER = "bash"        # "bash": call run.sh via bash -lc; "direct": call executable directly
WORK_DIR  = "."
YEARS     = 42
VAR_NAME  = "ubar"
OUT_BASENAME = "uout_real_urad"
ERA_CLIM_PATH = "era_10hpa_clim.npy"
TIMEOUT_SEC = 43200            # 12 hours

# Vertical grid parameters (consistent with Fortran)
DZ  = 100.0
H   = 7000.0
Z0  = 10000.0
KMAX = 1001

# --------------------  --------------------
def k_index_for_10hpa(dz=DZ, h=H, z0=Z0, kmax=KMAX):
    z10 = h * math.log(1000.0/10.0)
    k10 = int(round((z10 - z0)/dz)) + 1
    return max(2, min(kmax-1, k10))  # 1-based

def tag_from_params(beta, s_tau):
    b = f"{beta:.3f}".replace(".", "p")
    s = f"{s_tau:.3f}".replace(".", "p")
    return f"b{b}_s{s}"

def run_once(exe, runner, beta, s_tau, outtag=None, workdir=WORK_DIR, timeout_sec=TIMEOUT_SEC):
    """Call Fortran: supports two runner modes
       - runner='bash':  bash -lc './run.sh beta s_tau tag'
       - runner='direct':  './d.out beta s_tau tag'
    """
    if outtag is None:
        outtag = tag_from_params(beta, s_tau)

    if runner == "bash":
        cmdline = f'{exe} {beta} {s_tau} {outtag}'
        res = subprocess.run(
            ["bash", "-lc", cmdline],
            cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=timeout_sec
        )
    elif runner == "direct":
        res = subprocess.run(
            [exe, str(beta), str(s_tau), outtag],
            cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=timeout_sec
        )
    else:
        raise ValueError(f"Unknown runner: {runner}")

    if res.returncode != 0:
        print("=== STDOUT ===\n", res.stdout)
        print("=== STDERR ===\n", res.stderr)
        raise RuntimeError(f"Model run failed (beta={beta}, s_tau={s_tau}).")

    # “ tag”“ tag” Fortran
    nc_tag = Path(workdir) / f"{OUT_BASENAME}_{outtag}.nc"
    if nc_tag.exists():
        return nc_tag
    nc_plain = Path(workdir) / f"{OUT_BASENAME}.nc"
    if nc_plain.exists():
        return nc_plain

    print("=== STDOUT ===\n", res.stdout)
    print("=== STDERR ===\n", res.stderr)
    raise FileNotFoundError(f"Expected output not found: {nc_tag} or {nc_plain}")

def load_model_clim_10hpa(nc_path, years=YEARS, var=VAR_NAME, k10=None):
    if k10 is None:
        k10 = k_index_for_10hpa() - 1   # 0-based
    with Dataset(nc_path, "r") as ds:
        arr = ds.variables[var][:]
    if arr.ndim != 2:
        raise ValueError(f"{nc_path} variable {var} has ndim={arr.ndim}, expected 2")

    #  (height,time) 
    if arr.shape[0] != years and arr.shape[1] == years:
        arr = arr.T

    time_len, kmax_found = arr.shape
    if kmax_found != KMAX:
        print(f"[WARN] KMAX mismatch in file: {kmax_found} vs expected {KMAX}")
    assert time_len % years == 0, f"time_len={time_len} not divisible by years={years}"
    n_per_year = time_len // years
    arr3 = arr.reshape(years, n_per_year, kmax_found)
    clim = arr3.mean(axis=0)            # (n_per_year, k)
    return clim[:, k10], n_per_year

def load_era_target(n_per_year, path=ERA_CLIM_PATH):
    y = np.load(path)
    if len(y) == n_per_year:
        return y
    x0 = np.linspace(0.0, 1.0, len(y),     endpoint=False)
    x1 = np.linspace(0.0, 1.0, n_per_year, endpoint=False)
    return np.interp(x1, x0, y)

# --------------------  --------------------
class GreyBoxBetaSTau:
    def __init__(self, exe=DEFAULT_EXE, runner=DEFAULT_RUNNER, years=YEARS, var=VAR_NAME, cache=True):
        self.exe = exe
        self.runner = runner
        self.years = years
        self.var   = var
        self.use_cache = cache
        self._cache = {}

    def simulate(self, beta, s_tau):
        key = (round(beta,6), round(s_tau,6))
        if self.use_cache and key in self._cache:
            return self._cache[key]
        outtag = tag_from_params(beta, s_tau)
        u_nc   = run_once(self.exe, self.runner, beta, s_tau, outtag=outtag)
        clim10, n_per_year = load_model_clim_10hpa(u_nc, years=self.years, var=self.var)
        if self.use_cache:
            self._cache[key] = (clim10, n_per_year, outtag, str(u_nc))
        return clim10, n_per_year, outtag, str(u_nc)

    def residuals(self, x, era_target_override=None):
        beta, s_tau = float(x[0]), float(x[1])
        beta  = min(1.5, max(0.2, beta))
        s_tau = min(3.0, max(0.5, s_tau))
        clim10, n_per_year, _, _ = self.simulate(beta, s_tau)
        tgt = era_target_override if era_target_override is not None else load_era_target(n_per_year)
        return float(np.sqrt(np.mean((clim10 - tgt)**2)))

    def estimate_parameters(self, x0, bounds=((0.2, 1.5), (0.5, 3.0)), method="Powell", era_target_override=None, maxiter=60):
        res = minimize(
            fun=lambda x: self.residuals(x, era_target_override=era_target_override),
            x0=np.array(x0, dtype=float),
            method=method,
            bounds=bounds,
            options={"maxiter": maxiter, "disp": True}
        )
        self.params_ = res.x
        self.fun_    = res.fun
        self.result_ = res
        return res

# --------------------  --------------------
def main():

    global ERA_CLIM_PATH
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=DEFAULT_EXE, help="run.sh  d.out")
    ap.add_argument("--runner", default=DEFAULT_RUNNER, choices=["bash","direct"], help="bash: bash -lc  run.shdirect:  exe")
    ap.add_argument("--years", type=int, default=YEARS)
    ap.add_argument("--var", default=VAR_NAME)
    ap.add_argument("--era", default=ERA_CLIM_PATH)
    ap.add_argument("--bounds", default=None, help=" '(0.3,1.4),(0.8,2.8)'")
    ap.add_argument("--method", default="Powell")
    ap.add_argument("--maxiter", type=int, default=50)

    ap.add_argument("--grid-beta", default="0.6,1.1,6", help="beta : start,end,count")
    ap.add_argument("--grid-stau", default="1.0,2.2,7", help="s_tau : start,end,count")
    ap.add_argument("--smoke", action="store_true", help="")
    ap.add_argument("--beta", type=float, default=0.8)
    ap.add_argument("--stau", type=float, default=1.5)
    ap.add_argument("--out", default="fit_result.json")
    args = ap.parse_args()

    ERA_CLIM_PATH = args.era

    print("[INFO] CWD:", os.getcwd())
    print("[INFO] exe:", args.exe, "| runner:", args.runner)

    gb = GreyBoxBetaSTau(exe=args.exe, runner=args.runner, years=args.years, var=args.var, cache=True)

    #  & 
    if args.smoke:
        try:
            clim10, n_per_year, tag, nc_path = gb.simulate(args.beta, args.stau)
            print(f"[SMOKE] OK. tag={tag}, n_per_year={n_per_year}, output={nc_path}")
            return 0
        except Exception as e:
            print("[SMOKE] FAILED:", repr(e))
            return 2

    # 
    b0,b1,bn = map(float, args.grid_beta.split(","))
    s0,s1,sn = map(float, args.grid_stau.split(","))
    bn = int(bn); sn = int(sn)

    beta_grid  = np.linspace(b0, b1, bn)
    stau_grid  = np.linspace(s0, s1, sn)

    best = (1e9, None, None)
    for b in beta_grid:
        for s in stau_grid:
            try:
                val = gb.residuals([b, s])
                if val < best[0]:
                    best = (val, b, s)
                print(f"[grid] RMSE={val:.5f} at beta={b:.3f}, s_tau={s:.3f}")
            except Exception as e:
                print(f"[grid] FAILED at beta={b}, s_tau={s} -> {repr(e)}")
    x0 = [best[1], best[2]]
    print(f"\n[grid best] x0 = beta={x0[0]:.3f}, s_tau={x0[1]:.3f}, RMSE={best[0]:.5f}\n")

    # 
    bounds = ((0.3, 1.4), (0.8, 2.8)) if args.bounds is None else eval(args.bounds)
    res = gb.estimate_parameters(x0=x0, bounds=bounds, method=args.method, maxiter=args.maxiter)
    b_opt, s_opt = map(float, res.x)

    # 
    _, n_per_year, tag, nc_path = gb.simulate(b_opt, s_opt)

    result = {
        "beta": b_opt,
        "s_tau": s_opt,
        "RMSE": float(res.fun),
        "status": str(res.message),
        "n_per_year": int(n_per_year),
        "tag": tag,
        "uout_nc": nc_path,
        "bounds": bounds,
        "method": args.method,
        "years": args.years
    }
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)
    print("\n===== BEST =====")
    print(json.dumps(result, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
