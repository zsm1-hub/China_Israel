#!/usr/bin/env python
# coding: utf-8

# # `homo_gridblock_hit.m` — Python/Jupyter port
# Numerically matched port using `gridblock_jupyter.py`. It retains raw Euclidean HIT pair distances (no periodic minimum-image correction, matching the MATLAB helper), fixed 4×4 blocks, per-block particle limiting, streaming accumulators, and MAT field names.

# In[8]:


import time
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.io import savemat

from tqdm import tqdm

import gridblock_jupyter as gb


def run_hit_eulerian_scales(cfg):
    """
    HIT Lagrangian statistics evaluated at Eulerian reference scales.

    Required cfg keys:
        input_dir
        nparticles_file
        num_to_select
        seconds
        dt
        timerange_matlab
        nblock_x
        nblock_y
        min_pairs
        min_valid_blocks
        r_requested
    """

    t0 = time.perf_counter()

    rng = np.random.default_rng(
        cfg.get("random_seed", None)
    )

    tr = (
        np.asarray(
            cfg["timerange_matlab"],
            dtype=int
        )
        - 1
    )

    nparticles_file = int(
        cfg["nparticles_file"]
    )

    num_to_select = int(
        cfg["num_to_select"]
    )

    r_requested = np.asarray(
        cfg["r_requested"],
        dtype=float
    ).ravel()

    r_requested = np.sort(
        np.unique(
            r_requested[
                np.isfinite(r_requested)
                & (r_requested > 0)
            ]
        )
    )

    if r_requested.size < 2:
        raise ValueError(
            "r_requested must contain at least two positive scales."
        )

    # --------------------------------------------------------
    # Eulerian scales define logarithmic pair-distance bins
    # --------------------------------------------------------

    log_r = np.log(r_requested)

    log_edges = np.empty(
        r_requested.size + 1
    )

    log_edges[1:-1] = (
        0.5
        * (log_r[:-1] + log_r[1:])
    )

    log_edges[0] = (
        log_r[0]
        - 0.5
        * (log_r[1] - log_r[0])
    )

    log_edges[-1] = (
        log_r[-1]
        + 0.5
        * (log_r[-1] - log_r[-2])
    )

    dist_bin = np.exp(log_edges)
    dist_axis = r_requested
    nscale = len(dist_axis)

    # --------------------------------------------------------
    # Read HIT trajectory data
    # --------------------------------------------------------

    input_dir = Path(
        cfg["input_dir"]
    )

    nc_file = (
        input_dir
        / f'HIT2d_pars_P{nparticles_file}'
          f'T{cfg["seconds"]:.1f}seconds.nc'
    )

    print("Opening:")
    print(nc_file)

    with xr.open_dataset(
        nc_file,
        decode_times=False
    ) as ds:

        lon_all = gb.matlab_layout(
            ds.lon.values,
            nparticles_file,
            "lon"
        )

        lat_all = gb.matlab_layout(
            ds.lat.values,
            nparticles_file,
            "lat"
        )

        u_all = gb.matlab_layout(
            ds.ue.values,
            nparticles_file,
            "ue"
        )

        v_all = gb.matlab_layout(
            ds.ve.values,
            nparticles_file,
            "ve"
        )

        available = lon_all.shape[1]

        if num_to_select > available:
            raise ValueError(
                f"num_to_select={num_to_select} "
                f"but only {available} particles are available."
            )

        selected = rng.choice(
            available,
            size=num_to_select,
            replace=False
        )

        ix = np.ix_(
            tr,
            selected
        )

        lon = lon_all[ix]
        lat = lat_all[ix]
        u = u_all[ix]
        v = v_all[ix]

    ntime = len(tr)

    # --------------------------------------------------------
    # Spatial blocks
    # --------------------------------------------------------

    x_min = np.nanmin(lon)
    x_max = np.nanmax(lon)
    y_min = np.nanmin(lat)
    y_max = np.nanmax(lat)

    x_edges = np.linspace(
        x_min,
        x_max,
        cfg["nblock_x"] + 1
    )

    y_edges = np.linspace(
        y_min,
        y_max,
        cfg["nblock_y"] + 1
    )

    nblock = (
        cfg["nblock_x"]
        * cfg["nblock_y"]
    )

    # --------------------------------------------------------
    # Accumulators
    # --------------------------------------------------------

    sums = [
        np.zeros(nscale),  # Dl
        np.zeros(nscale),  # Dt
        np.zeros(nscale),  # Dll
        np.zeros(nscale),  # Dtt
        np.zeros(nscale),  # Dlll
        np.zeros(nscale),  # Dltt
    ]

    local_sums = [
        np.zeros((nscale, nblock)),
        np.zeros((nscale, nblock)),
        np.zeros((nscale, nblock)),
        np.zeros((nscale, nblock)),
        np.zeros((nscale, nblock)),
        np.zeros((nscale, nblock)),
    ]

    total_count = np.zeros(
        nscale,
        dtype=np.int64
    )

    block_count = np.zeros(
        (nscale, nblock),
        dtype=np.int64
    )

    # --------------------------------------------------------
    # Progress iterator
    # --------------------------------------------------------

    if tqdm is not None:
        iterator = tqdm(
            range(ntime),
            total=ntime,
            desc="HIT Lagrangian pair statistics",
            unit="time"
        )
    else:
        iterator = range(ntime)

    # --------------------------------------------------------
    # Main calculation
    # --------------------------------------------------------

    for it in iterator:

        valid = (
            np.isfinite(lon[it])
            & np.isfinite(lat[it])
            & np.isfinite(u[it])
            & np.isfinite(v[it])
        )

        ids = np.flatnonzero(valid)

        if ids.size < 2:
            continue

        x = lon[it, ids]
        y = lat[it, ids]
        uu = u[it, ids]
        vv = v[it, ids]

        # HIT coordinates are nondimensional
        ii, jj, distance, dl, dt = gb.condensed_pairs(
            x,
            y,
            uu,
            vv,
            geographic=False
        )

        # Pair midpoint determines spatial block
        mid_x = 0.5 * (
            x[ii] + x[jj]
        )

        mid_y = 0.5 * (
            y[ii] + y[jj]
        )

        bx = gb.discretize(
            mid_x,
            x_edges
        )

        by = gb.discretize(
            mid_y,
            y_edges
        )

        block = (
            bx * cfg["nblock_y"]
            + by
        )

        scale_bin = gb.discretize(
            distance,
            dist_bin
        )

        good = (
            (scale_bin >= 0)
            & (block >= 0)
            & (block < nblock)
            & np.isfinite(dl)
            & np.isfinite(dt)
        )

        scale_bin = scale_bin[good]
        block = block[good]
        dl = dl[good]
        dt = dt[good]

        if len(dl) == 0:
            continue

        values = [
            dl,
            dt,
            dl ** 2,
            dt ** 2,
            dl ** 3,
            dl * dt ** 2,
        ]

        # Global statistics
        for target, bins, values_i in zip(
            sums,
            [scale_bin] * len(values),
            values
        ):
            gb.add_by_bin(
                target,
                bins,
                values_i
            )

        gb.add_by_bin(
            total_count,
            scale_bin,
            np.ones(len(scale_bin))
        )

        # Block statistics
        flat_index = (
            scale_bin * nblock
            + block
        )

        for target, values_i in zip(
            local_sums,
            values
        ):
            np.add.at(
                target.ravel(),
                flat_index,
                values_i
            )

        np.add.at(
            block_count.ravel(),
            flat_index,
            1
        )

        if tqdm is None and (
            (it + 1) % 10 == 0
            or it == ntime - 1
        ):
            print(
                f"time {it + 1}/{ntime}"
            )

    # --------------------------------------------------------
    # Global structure functions
    # --------------------------------------------------------

    global_sf = []

    for total in sums:

        value = np.full(
            nscale,
            np.nan
        )

        valid = total_count > 0

        value[valid] = (
            total[valid]
            / total_count[valid]
        )

        global_sf.append(value)

    (
        Dl,
        Dt,
        Dll,
        Dtt,
        Dlll,
        Dltt
    ) = global_sf

    D1 = Dl + Dt
    D2 = Dll + Dtt
    D3 = Dlll + Dltt

    # --------------------------------------------------------
    # Local block structure functions
    # --------------------------------------------------------

    local_sf = []

    valid_block = (
        block_count
        >= cfg["min_pairs"]
    )

    for total in local_sums:

        value = np.full(
            (nscale, nblock),
            np.nan
        )

        value[valid_block] = (
            total[valid_block]
            / block_count[valid_block]
        )

        local_sf.append(value)

    (
        local_Dl,
        local_Dt,
        local_Dll,
        local_Dtt,
        local_Dlll,
        local_Dltt
    ) = local_sf

    # --------------------------------------------------------
    # Homogeneity metrics
    # --------------------------------------------------------

    def calculate_H(local_values):

        H = np.full(
            nscale,
            np.nan
        )

        mean = np.full(
            nscale,
            np.nan
        )

        std = np.full(
            nscale,
            np.nan
        )

        rms = np.full(
            nscale,
            np.nan
        )

        nvalid = np.zeros(
            nscale,
            dtype=int
        )

        for ir in range(nscale):

            valid = (
                block_count[ir]
                >= cfg["min_pairs"]
            ) & np.isfinite(
                local_values[ir]
            )

            values = local_values[
                ir,
                valid
            ]

            nvalid[ir] = len(values)

            if len(values) < cfg[
                "min_valid_blocks"
            ]:
                continue

            mean[ir] = np.mean(
                values
            )

            std[ir] = np.std(
                values,
                ddof=0
            )

            rms[ir] = np.sqrt(
                np.mean(values ** 2)
            )

            denominator = (
                abs(mean[ir])
                + rms[ir]
            )

            if denominator > 0:
                H[ir] = (
                    std[ir]
                    / denominator
                )

        return H, mean, std, rms, nvalid

    H1L, mean_Dl, std_Dl, rms_Dl, nvalid_Dl = (
        calculate_H(local_Dl)
    )

    H1T, mean_Dt, std_Dt, rms_Dt, nvalid_Dt = (
        calculate_H(local_Dt)
    )

    H2L, mean_Dll, std_Dll, rms_Dll, nvalid_Dll = (
        calculate_H(local_Dll)
    )

    H2T, mean_Dtt, std_Dtt, rms_Dtt, nvalid_Dtt = (
        calculate_H(local_Dtt)
    )

    H3L, mean_Dlll, std_Dlll, rms_Dlll, nvalid_Dlll = (
        calculate_H(local_Dlll)
    )

    H3LTT, mean_Dltt, std_Dltt, rms_Dltt, nvalid_Dltt = (
        calculate_H(local_Dltt)
    )

    H2, _, _, _, _ = calculate_H(
        local_Dll + local_Dtt
    )

    H3, _, _, _, _ = calculate_H(
        local_Dlll + local_Dltt
    )

    runtime = time.perf_counter() - t0

    # --------------------------------------------------------
    # Save output
    # --------------------------------------------------------

    output = {
        "Case": "HIT2d",
        "nparticles_file": nparticles_file,
        "num_to_select": num_to_select,
        "selected_indices": selected + 1,
        "seconds": cfg["seconds"],
        "dt": cfg["dt"],
        "timerange": tr + 1,

        # Eulerian reference scales
        "r_requested": r_requested,
        "dist_axis": dist_axis,
        "dist_bin": dist_bin,

        "x_domain": np.array([x_min, x_max]),
        "y_domain": np.array([y_min, y_max]),
        "x_edges": x_edges,
        "y_edges": y_edges,

        "nblock_x": cfg["nblock_x"],
        "nblock_y": cfg["nblock_y"],
        "nBlock": nblock,

        "min_pairs": cfg["min_pairs"],
        "min_valid_blocks": cfg["min_valid_blocks"],

        "Dl": Dl,
        "Dt": Dt,
        "D1": D1,
        "Dll": Dll,
        "Dtt": Dtt,
        "D2": D2,
        "Dlll": Dlll,
        "Dltt": Dltt,
        "D3": D3,

        "local_Dl": local_Dl,
        "local_Dt": local_Dt,
        "local_Dll": local_Dll,
        "local_Dtt": local_Dtt,
        "local_Dlll": local_Dlll,
        "local_Dltt": local_Dltt,

        "pair_count": block_count,
        "count_total": total_count,

        "H1L": H1L,
        "H1T": H1T,
        "H2L": H2L,
        "H2T": H2T,
        "H2": H2,
        "H3L": H3L,
        "H3LTT": H3LTT,
        "H3": H3,

        "mean_Dl": mean_Dl,
        "std_Dl": std_Dl,
        "rms_Dl": rms_Dl,
        "nvalid_Dl": nvalid_Dl,

        "mean_Dt": mean_Dt,
        "std_Dt": std_Dt,
        "rms_Dt": rms_Dt,
        "nvalid_Dt": nvalid_Dt,

        "mean_Dll": mean_Dll,
        "std_Dll": std_Dll,
        "rms_Dll": rms_Dll,
        "nvalid_Dll": nvalid_Dll,

        "mean_Dtt": mean_Dtt,
        "std_Dtt": std_Dtt,
        "rms_Dtt": rms_Dtt,
        "nvalid_Dtt": nvalid_Dtt,

        "mean_Dlll": mean_Dlll,
        "std_Dlll": std_Dlll,
        "rms_Dlll": rms_Dlll,
        "nvalid_Dlll": nvalid_Dlll,

        "mean_Dltt": mean_Dltt,
        "std_Dltt": std_Dltt,
        "rms_Dltt": rms_Dltt,
        "nvalid_Dltt": nvalid_Dltt,

        "runtime_seconds": runtime,
    }

    output_path = (
        input_dir
        / f"HIT2d_pars_P{num_to_select}"
          f"T{tr[-1] + 1}_EulerianScalesHL.mat"
    )

    savemat(
        output_path,
        output,
        do_compression=True,
        oned_as="row"
    )

    print("\nFinished.")
    print(f"Runtime: {runtime / 60:.2f} min")
    print("Saved:")
    print(output_path)

    return output, output_path


# In[9]:


from pathlib import Path
import numpy as np

# 使用 notebook 中重新定义的函数
# run_hit_eulerian_scales

dx_hit = 2.0 * np.pi / 512.0

# Eulerian reference scales
r_requested_hit = np.geomspace(
    2.0 * dx_hit,
    np.pi,
    25
)



cfg_hit = dict(
    input_dir=Path(
        "/meddy/simingzhang/Data/Parcels_data/"
        "HIT2d_rough"
    ),

    nparticles_file=62500,
    num_to_select=30000,

    seconds=20.0,
    dt=0.1,

    # 第 50–150 个时间戳，包含两端
    timerange_matlab=np.arange(
        50,
        151
    ),

    nblock_x=4,
    nblock_y=4,

    min_pairs=500,
    min_valid_blocks=8,

    random_seed=1,

    r_requested=r_requested_hit
)

result_hit, output_hit = run_hit_eulerian_scales(
    cfg_hit
)

result_hit, output_hit = run_hit_eulerian_scales(
    cfg_hit
)

print("HIT output:")
print(output_hit)


# ## Random-selection note
# The seed is retained, but NumPy and MATLAB use different random-number/permutation implementations, so `selected_indices` will not be identical solely from the same seed. To compare every output value exactly, load MATLAB's saved `selected_indices` in Python; the deterministic calculations after selection use the same formulas and bin rules.
