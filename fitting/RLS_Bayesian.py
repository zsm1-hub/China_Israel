#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# ============================================================
# Fk_fitting_SF3_Bayesian_RLS_Lcurve_uncertainty4
# ============================================================
import numpy as np
import xarray as xr
from scipy.special import jv
from scipy.linalg import inv
import os
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.linalg import solve

def get_optimal_prior_from_Lcurve(
    Y, W, A, kf_final,
    fac0_range=np.logspace(-20, -5, 80),
    po_range=np.logspace(-20, -5, 80)
):
    optimal_fac, optimal_po, Lcurve = L_curve_analysis(
        Y, W, A, kf_final, fac0_range, po_range
    )

    print("L-curve prior：")
    print(f"  fac0 = {optimal_fac}")
    print(f"  po   = {optimal_po}")

    return optimal_fac, optimal_po, Lcurve   


def Fk_fitting_SF3_Bayesian_RLS_Lcurve_uncertainty4(
    fname, mindist, maxdist, dot,
    n_bootstrap=1, bootstrap_size=900
):
    data = sio.loadmat(fname)

    SF3 = np.asarray(data["SF3"])
    dist_axis = np.asarray(data["dist_axis"]).ravel()

    n_total_samples = SF3.shape[1]

    # -------- 波数 --------
    kf_temp = np.logspace(np.log10(1/maxdist), np.log10(1/mindist), dot)*2*np.pi 
    kf_temp = kf_temp.ravel()
    dkf_temp = np.diff(kf_temp)
    kf_final = 0.5 * (kf_temp[:-1] + kf_temp[1:])

    n_kf = len(kf_final)
    n_R = len(dist_axis)

    # -------- Bootstrap 容器 --------
    class BootstrapResults:
        pass

    res = BootstrapResults()
    res.eps_all = np.zeros((n_kf + 1, n_bootstrap))
    res.Fk_all = np.zeros((n_kf, n_bootstrap))
    res.Vt_all = np.zeros((n_R, n_bootstrap))
    res.residual = np.zeros((n_R, n_bootstrap))
    res.eps_error = np.zeros((n_kf + 1, n_bootstrap))
    res.Fk_error = np.zeros((n_kf, n_bootstrap))
    res.optimal_fac_all = np.zeros(n_bootstrap)
    res.optimal_po_all = np.zeros(n_bootstrap)
    res.kf_final = kf_final.copy()

    # =========================================================
    # ✅ Step 1：只用原始数据做一次 L-curve（不 bootstrap）
    # =========================================================
    Y_full = np.nanmean(SF3, axis=1)
    Y_std_full = np.nanstd(SF3, axis=1)
    Y_std_full = np.clip(Y_std_full, 1e-12, None)

    W_full = np.diag(Y_std_full**2)
    A = defA(dist_axis, kf_final, dkf_temp)

    # optimal_fac, optimal_po = get_optimal_prior_from_Lcurve(
    #     Y_full, W_full, A, kf_final
    # )
    optimal_fac, optimal_po, Lcurve = get_optimal_prior_from_Lcurve(
        Y_full, W_full, A, kf_final
    )

    # =========================================================
    # ✅ Step 2：Bootstrap 中直接用这两个常数
    # =========================================================
    for kk in range(n_bootstrap):
        print(f"Bootstrap {kk+1}/{n_bootstrap}")

        idx = np.random.choice(n_total_samples, bootstrap_size, replace=True)
        SF3_boot = SF3[:, idx]

        Y = np.nanmean(SF3_boot, axis=1)
        Y_std = np.nanstd(SF3_boot, axis=1)
        Y_std = np.clip(Y_std, 1e-12, None)

        W = np.diag(Y_std**2)

        # ✅ 不再扫描
        P_diag = np.hstack([[optimal_fac], optimal_po * np.ones(n_kf)])
        P = np.diag(P_diag)

        eps, Vt, residual, Cxx = RLS(Y, W, P, A)
        Fk = calcFk(eps, kf_final, dkf_temp)
        eps_error = np.sqrt(np.diag(Cxx))

        H = defH(kf_final, dkf_temp);
        Fxx = errorsFlux(Cxx, H);
        Fk_error=np.sqrt(np.diag(Fxx));

        res.Vt_all[:, kk] = Vt
        res.eps_all[:, kk] = eps
        res.Fk_all[:, kk] = Fk
        res.residual[:,kk] = residual
        res.eps_error[:,kk] = eps_error
        res.Fk_error[:,kk] = Fk_error
        
        res.optimal_fac_all[kk] = optimal_fac
        res.optimal_po_all[kk] = optimal_po
        res.Lcurve = Lcurve
        

    return dist_axis, Y_full, Y_std_full, optimal_fac, optimal_po, res

def L_curve_analysis(Y, W, A, kbins, fac0_range, po_range):
    nk = len(kbins)
    n_po = len(po_range)
    n_fac0 = len(fac0_range)

    norm2_eps = np.zeros((n_po, n_fac0))
    norm2_res = np.zeros((n_po, n_fac0))

    wgts = np.sqrt(np.diag(W))

    for n in range(n_po):
        for m in range(n_fac0):
            P_diag = np.hstack([[fac0_range[m]], po_range[n] * np.ones(nk)])
            P_test = np.diag(P_diag)

            try:
                eps_test, Vt_test, res_test, _ = RLS(Y, W, P_test, A)

                norm2_eps[n, m] = np.linalg.norm(eps_test[1:], 2)
                norm2_res[n, m] = np.linalg.norm(res_test / wgts, 2)

            except Exception:
                norm2_eps[n, m] = np.inf
                norm2_res[n, m] = np.inf

    po_idx, fac0_idx = find_L_curve_corner(
        norm2_res, norm2_eps, n_po, n_fac0
    )

    # ✅ 打包 L-curve 数据（与 return 同级）
    Lcurve = {
        "residual_norms": norm2_res,
        "solution_norms": norm2_eps,
        "po_range": po_range,
        "fac0_range": fac0_range,
        "po_idx": po_idx,
        "fac0_idx": fac0_idx
    }

    return fac0_range[fac0_idx], po_range[po_idx], Lcurve


def find_L_curve_corner(residual_norms, solution_norms, n_po, n_fac0):
    log_res = np.log10(residual_norms.flatten())
    log_sol = np.log10(solution_norms.flatten())

    valid = np.isfinite(log_res) & np.isfinite(log_sol)
    log_res = log_res[valid]
    log_sol = log_sol[valid]

    if len(log_res) < 3:
        return n_po // 2, n_fac0 // 2

    curvature = np.zeros_like(log_res)

    for i in range(1, len(log_res) - 1):
        dx1 = log_res[i] - log_res[i - 1]
        dy1 = log_sol[i] - log_sol[i - 1]
        dx2 = log_res[i + 1] - log_res[i]
        dy2 = log_sol[i + 1] - log_sol[i]

        denom = (dx1**2 + dy1**2) ** 1.5
        if denom > 0:
            curvature[i] = abs(dx1 * dy2 - dx2 * dy1) / denom

    max_idx = np.argmax(curvature)

    po_idx = max_idx % n_po
    fac0_idx = max_idx // n_po

    return po_idx, fac0_idx

def calcFk(eps, k, dk):
    eps = np.asarray(eps).ravel()
    k = np.asarray(k).ravel()
    dk = np.asarray(dk).ravel()

    nk = len(k)
    Fk = np.zeros(nk)

    Fk[0] = -eps[0]
    for jj in range(nk - 1):
        Fk[jj + 1] = Fk[jj] + eps[jj + 1] * dk[jj]

    return Fk

from scipy.linalg import inv

def RLS(Y, W, P, X):
    Y = np.asarray(Y).ravel()
    X = np.asarray(X)

    # -------- W --------
    if np.ndim(W) == 1:
        W = np.diag(W)
    W = np.asarray(W)

    if W.shape[0] != W.shape[1]:
        raise ValueError("W 不是方阵")
    if W.shape[0] != X.shape[0]:
        raise ValueError("diag(W) 长度与 X 行数不匹配")

    # -------- P --------
    if np.ndim(P) == 0:
        P = P * np.eye(X.shape[1])
    elif np.ndim(P) == 1:
        if len(P) == 1:
            P = P * np.eye(X.shape[1])
        elif len(P) == X.shape[1]:
            P = np.diag(P)
        else:
            raise ValueError("diag(P) 长度与 X 列数不匹配")
    elif np.ndim(P) == 2:
        if P.shape[0] != P.shape[1]:
            raise ValueError("P 不是方阵")
        if P.shape[0] != X.shape[1]:
            raise ValueError("diag(P) 长度与 X 列数不匹配")
    else:
        raise ValueError("P 维度非法")

    # -------- Y --------
    if len(Y) != X.shape[0]:
        raise ValueError("Y 长度与 X 行数不匹配")

    # -------- 核心计算（完全 MATLAB 行为）--------
    if np.all(P == 0):
        Pinv = P
    else:
        Pinv = inv(P)

    Winv = inv(W)

    M = X.T @ Winv @ X
    Minv = inv(M + Pinv)
    N = Minv @ X.T @ Winv

    x0 = N @ Y
    y0 = X @ x0
    n0 = Y - y0
    Cxx = Minv

    return x0, y0, n0, Cxx

def defA(r, k, dk):
    r = np.asarray(r).ravel()
    k = np.asarray(k).ravel()
    dk = np.asarray(dk).ravel()

    nr = len(r)
    nk = len(k)

    A = np.zeros((nr, nk + 1))

    for jj in range(nk):
        for ii in range(nr):
            A[ii, jj + 1] = -4.0 * jv(1, k[jj] * r[ii]) / k[jj] * dk[jj]

    A[:, 0] = 2.0 * r
    return A

def calc_radial(S3L1, S3T1, N, xscale):
    """
    MATLAB calc_radial 的 Python/Jupyter 等价实现
    """
    x = np.arange(-N//2, N//2)
    X, Y = np.meshgrid(x, x)

    R = np.sqrt(X**2 + Y**2)

    r_max = R.max()
    dr = 1.0
    r_edges = np.arange(0, r_max + dr, dr)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])

    S3L_radial = np.zeros_like(r_centers)
    S3T_radial = np.zeros_like(r_centers)
    counts = np.zeros_like(r_centers)

    for i in range(len(r_centers)):
        mask = (R >= r_edges[i]) & (R < r_edges[i + 1])
        if np.any(mask):
            S3L_radial[i] = np.mean(S3L1[mask])
            S3T_radial[i] = np.mean(S3T1[mask])
            counts[i] = np.sum(mask)

    valid_mask = counts > 0
    r = r_centers[valid_mask]
    S3L1 = S3L_radial[valid_mask]
    S3T1 = S3T_radial[valid_mask]
    SF3 = S3L1 + S3T1
    r = r * xscale

    return r, SF3, S3L1, S3T1

def defH(k, dk):
    """
    从 KE 注入率误差 → KE 通量误差 的变换矩阵
    """
    k = np.asarray(k).ravel()
    dk = np.asarray(dk).ravel()

    nk = len(k)
    H = np.zeros((nk, nk + 1))
    Hlog = np.zeros((nk, nk))

    # Heaviside step: k >= k(i)
    for i in range(nk):
        Hlog[:, i] = (k >= k[i]).astype(float)

    # 构造 H
    for jj in range(nk):
        for ii in range(nk):
            H[jj, ii + 1] = Hlog[jj, ii] * dk[jj]

    H[:, 0] = -1.0
    return H

def errorsFlux(Cxx, H):
    """
    从 KE 注入率误差协方差 → KE 通量误差协方差
    """
    Cxx = np.asarray(Cxx)
    H = np.asarray(H)

    # 通量协方差矩阵
    Fxx = H @ Cxx @ H.T
    return Fxx

def bin_SF3_by_dist_axis(dist_axis, r, SF3):
    """
    用 dist_axis 作为分箱中心，对 SF3 做径向平均
    """
    dist_axis = np.asarray(dist_axis)
    r = np.asarray(r)
    SF3 = np.asarray(SF3)

    # 构造 bin edges（假设 dist_axis 均匀）
    dr = np.mean(np.diff(dist_axis))
    edges = np.hstack([
        dist_axis[0] - dr/2,
        dist_axis[:-1] + np.diff(dist_axis)/2,
        dist_axis[-1] + dr/2
    ])

    # 每个 r 属于哪个 bin
    bin_idx = np.digitize(r, edges) - 1

    # 初始化
    SF3_binned = np.full(len(dist_axis), np.nan)
    counts = np.zeros(len(dist_axis), dtype=int)

    # 分箱平均
    for i in range(len(dist_axis)):
        mask = bin_idx == i
        if np.any(mask):
            SF3_binned[i] = np.mean(SF3[mask])
            counts[i] = np.sum(mask)

    return dist_axis, SF3_binned, counts


import numpy as np
from scipy.special import jv
from scipy.optimize import lsq_linear
from scipy.linalg import lstsq


import numpy as np
from scipy.special import jv
from scipy.linalg import lstsq

def RLS_Tikhonov(
    SF3,                  # (Nr, Nt)
    dist_axis,             # (Nr,)
    mindist,
    maxdist,
    inv_style,
    lam,
    npoints
):
    """
    Python version of Fk_fitting_SF3 (MATLAB)
    Returns:
        SpecFlux : (Nk, Nt)
        Vt       : (NR, Nt)
        ebs      : (Nk+1, Nt)
        kf       : (Nk,)
        lf       : R
    """

    SF3 = np.asarray(SF3)
    dist_axis = np.asarray(dist_axis).ravel()

    nsamps = SF3.shape[1]
    r = dist_axis
    Nr = len(r)

    # ---- select radial range ----
    ns = np.where(r >= mindist)[0]
    ne = np.where(r <= maxdist)[0]

    if ns.size == 0 or ne.size == 0:
        raise ValueError("mindist / maxdist 超出 dist_axis 范围")

    ns = ns[0]
    ne = ne[-1]

    R = r[ns:ne]
    NR = len(R)
    lf = R.copy()

    # ---- wavenumber grid ----
    # if kftype == 'log':
    #     kf = np.logspace(
    #         np.log10(1 / R.max()),
    #         np.log10(1 / R.min()),
    #         NR - 1
    #     )
    # else:
    #     kf = kf1.copy()
    kf=np.logspace(np.log10(1/maxdist),np.log10(1/mindist),npoints)*2*np.pi
    dk = np.diff(kf)
    kf = 0.5 * (kf[:-1] + kf[1:])

    # MATLAB code assumes kf decreasing
    kf = np.flip(kf)
    dk = np.flip(dk)

    Nk = len(kf)

    # ---- allocate ----
    ebs = np.zeros((Nk + 1, nsamps))
    Vt = np.zeros((NR, nsamps))
    SpecFlux = np.zeros((Nk, nsamps))

    norm_flag = 1

    # ---- loop over time / samples ----
    for n in range(nsamps):

        S = SF3[:, n]

        # select radial window
        V = S[ns:ne]

        # design matrix A (NR × (Nk+1))
        A = np.zeros((NR, Nk + 1))
        for j in range(Nk):
            A[:, j] = (-4 / kf[j]) * jv(1, kf[j] * R) * dk[j]

        A[:, -1] = 2 * R

        # ---- normalize by r ----
        if norm_flag == 1:
            A = A / R[:, None]
            V = V / R

        # ---- inversion ----
        if inv_style == 'NNLS':
            # 如需 NNLS 请安装: pip install scipy
            from scipy.optimize import lsq_linear
            res = lsq_linear(A, V, bounds=(0, np.inf))
            ebs[:, n] = res.x

        elif inv_style == 'LS':
            sol, _, _, _ = lstsq(A, V)
            ebs[:, n] = sol

        elif inv_style == 'RLS':
            V = V.reshape(-1)
            n_cols = A.shape[1]

            A_aug = np.vstack([
                A,
                np.sqrt(lam) * np.eye(n_cols)
            ])
            V_aug = np.concatenate([V, np.zeros(n_cols)])

            sol, _, _, _ = lstsq(A_aug, V_aug)
            ebs[:, n] = sol

        # ---- reconstruct Vt ----
        Vt[:, n] = 2 * ebs[-1, n] * R
        for j in range(Nk):
            Vt[:, n] -= (
                4 * ebs[j, n] / kf[j]
                * jv(1, kf[j] * R) * dk[j]
            )

        # ---- spectral flux (✅ 修复了 IndexError 的核心部分) ----
        SpecFlux[0, n] = -ebs[-1, n]

        for j in range(1, Nk):
            SpecFlux[j, n] = (
                SpecFlux[j - 1, n]
                + ebs[-(j + 1), n] * dk[-j]
            )

    SpecFlux = np.flipud(SpecFlux)

    return SpecFlux, Vt, ebs, kf, lf

def plot_L_curve(Lcurve, savepath=None):
    """
    Plot L-curve from precomputed Lcurve dictionary.
    Does NOT perform any computation.

    Parameters
    ----------
    Lcurve : dict
        Must contain:
        - residual_norms
        - solution_norms
        - po_range
        - fac0_range
        - po_idx
        - fac0_idx
    savepath : str, optional
        If given, save figure to path
    """
    residual_norms = Lcurve["residual_norms"]
    solution_norms = Lcurve["solution_norms"]
    po_range = Lcurve["po_range"]
    fac0_range = Lcurve["fac0_range"]
    po_idx = Lcurve["po_idx"]
    fac0_idx = Lcurve["fac0_idx"]

    fig, ax = plt.subplots(figsize=(6, 5))

    # ---- Scatter (L-curve) ----
    X = residual_norms.flatten()
    Y = solution_norms.flatten()

    sc = ax.scatter(
        X, Y,
        c=np.log10(Y),
        cmap="viridis",
        s=20,
        alpha=0.8
    )

    # ---- Corner point ----
    opt_res = residual_norms[po_idx, fac0_idx]
    opt_sol = solution_norms[po_idx, fac0_idx]

    ax.scatter(
        opt_res, opt_sol,
        c="red", s=150, marker="*",
        label="Selected corner"
    )

    # ---- Labels ----
    ax.set_xlabel(r"$||W^{-1/2}(Y - A\epsilon)||_2$")
    ax.set_ylabel(r"$||\epsilon||_2$")
    ax.set_title("L-curve for Bayesian RLS prior selection")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle="--", alpha=0.5)

    # ========== 关键修正：同时设置 X 和 Y 轴 ==========
    # 1. 裁剪 X 轴（最重要！）
    # 你的红点在 X=1~2 之间，把右边的大值切掉
    # ax.set_xlim(left=0.8, right=5)  # <-- 重点在这里

    # 2. 修正 Y 轴（合并为一行）
    # 根据你的数值 1e-8 和 1e-5，这个范围最合适
    # ax.set_ylim(bottom=1e-5, top=1e-3)  # <-- 去掉了重复的设置
    # ==================================================

    ax.legend(fontsize=9)
    plt.tight_layout()

    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
    plt.show()

