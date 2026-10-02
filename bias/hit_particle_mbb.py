import time
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.io import savemat
from tqdm import tqdm

import gridblock_jupyter as gb

from third_order_diagnostics import (
    init_third_order,
    update_third_order,
    finish_third_order,
)

def moving_time_block_bootstrap(
    time_sum,
    time_count,
    block_length,
    num_boot=1000,
    confidence=0.95,
    random_seed=12345,
):
    """Circular moving-block bootstrap of a count-weighted time mean.

    time_sum and time_count have shape (ntime, nscale). Each time record may
    contain a different number of particle pairs. Continuous time records are
    resampled in circular blocks, and each replicate is sum/sum(count), not an
    unweighted mean of per-time means.
    """
    time_sum=np.asarray(time_sum,float); time_count=np.asarray(time_count,float)
    if time_sum.shape != time_count.shape or time_sum.ndim != 2:
        raise ValueError("time_sum and time_count must have identical (ntime,nscale) shape")
    ntime,nscale=time_sum.shape; lengths=np.asarray(block_length,int).ravel()
    if lengths.size==1: lengths=np.full(nscale,int(lengths[0]))
    if lengths.size != nscale: raise ValueError("block_length must be scalar or one value per scale")
    rng=np.random.default_rng(random_seed)
    estimates=np.full((nscale,num_boot),np.nan)
    effective_blocks=np.zeros(nscale,int)
    for ir in range(nscale):
        length=int(np.clip(lengths[ir],1,ntime))
        effective_blocks[ir]=ntime//length
        # Fewer than two decorrelation blocks cannot provide an uncertainty estimate.
        if effective_blocks[ir] < 2 or np.sum(time_count[:,ir]) <= 0:
            continue
        n_draw=int(np.ceil(ntime/length))
        offsets=np.arange(length)
        for ib in range(num_boot):
            starts=rng.integers(0,ntime,size=n_draw)
            indices=((starts[:,None]+offsets[None,:])%ntime).ravel()[:ntime]
            denominator=np.sum(time_count[indices,ir])
            if denominator>0:
                estimates[ir,ib]=np.sum(time_sum[indices,ir])/denominator
    alpha=0.5*(1.0-confidence)
    return {
        "D3_mbb": estimates,
        "D3_mbb_mean": np.nanmean(estimates,axis=1),
        "D3_mbb_stderr": np.nanstd(estimates,axis=1,ddof=1),
        "D3_mbb_ci_low": np.nanpercentile(estimates,100*alpha,axis=1),
        "D3_mbb_ci_high": np.nanpercentile(estimates,100*(1-alpha),axis=1),
        "D3_mbb_block_length_time": lengths,
        "D3_mbb_effective_blocks": effective_blocks,
        "D3_mbb_num_boot": int(num_boot),
        "D3_mbb_confidence": float(confidence),
    }

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

    period = float(cfg.get("period", 2.0 * np.pi))
    lon = np.mod(lon, period)
    lat = np.mod(lat, period)
    x_min, x_max = 0.0, period
    y_min, y_max = 0.0, period

    x_edges = np.linspace(
        0.0,
        period,
        cfg["nblock_x"] + 1
    )

    y_edges = np.linspace(
        0.0,
        period,
        cfg["nblock_y"] + 1
    )

    nblock = (
        cfg["nblock_x"]
        * cfg["nblock_y"]
    )

    # --------------------------------------------------------
    # Accumulators
    # --------------------------------------------------------

    third_diag = init_third_order(
        nscale,
        reservoir_size=int(cfg.get("third_pdf_reservoir_size", 100_000)),
        random_seed=int(cfg.get("random_seed", 1)) + 12345,
    )

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

    do_mbb = bool(cfg.get("do_moving_block_bootstrap", False))
    if do_mbb:
        time_D2_sum = np.zeros((ntime, nscale), dtype=np.float64)
        time_D3_sum = np.zeros((ntime, nscale), dtype=np.float64)
        time_pair_count = np.zeros((ntime, nscale), dtype=np.int64)

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
        ii, jj, distance, dl, dt = gb.condensed_pairs_periodic(
            x,
            y,
            uu,
            vv,
            period=period
        )

        # Pair midpoint determines spatial block
        mid_x = gb.periodic_midpoint(
            x[ii], x[jj], period=period
        )

        mid_y = gb.periodic_midpoint(
            y[ii], y[jj], period=period
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

        update_third_order(third_diag, scale_bin, dl, dt)

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

        if do_mbb:
            time_pair_count[it] += np.bincount(
                scale_bin, minlength=nscale
            ).astype(np.int64)
            time_D2_sum[it] += np.bincount(
                scale_bin, weights=dl**2 + dt**2, minlength=nscale
            )
            time_D3_sum[it] += np.bincount(
                scale_bin, weights=dl**3 + dl*dt**2, minlength=nscale
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

    mbb_output = {}
    if do_mbb:
        # Turnover time tau_r=r/sqrt(D2); convert it to consecutive time records.
        tau_r = np.divide(
            dist_axis, np.sqrt(D2), out=np.full(nscale, np.nan),
            where=np.isfinite(D2) & (D2 > 0)
        )
        if cfg.get("mbb_block_length_time") is None:
            block_length_time = np.ceil(tau_r / float(cfg["dt"]))
            block_length_time[~np.isfinite(block_length_time)] = ntime
            block_length_time = np.maximum(block_length_time, 1).astype(int)
        else:
            supplied = np.asarray(cfg["mbb_block_length_time"], int).ravel()
            block_length_time = (
                np.full(nscale, int(supplied[0])) if supplied.size == 1 else supplied
            )
        mbb_output = moving_time_block_bootstrap(
            time_sum=time_D3_sum,
            time_count=time_pair_count,
            block_length=block_length_time,
            num_boot=int(cfg.get("mbb_num_boot", 1000)),
            confidence=float(cfg.get("mbb_confidence", 0.95)),
            random_seed=int(cfg.get("mbb_random_seed", cfg.get("random_seed", 1) + 271828)),
        )
        mbb_output.update({
            "D3_mbb_tau_r": tau_r,
            "D3_time_sum": time_D3_sum,
            "D2_time_sum": time_D2_sum,
            "pair_count_time": time_pair_count,
        })

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

    output["do_moving_block_bootstrap"] = do_mbb
    output.update(mbb_output)

    third_output = finish_third_order(third_diag)
    output.update(third_output)

    default_name = (
        f"HIT2d_pars_P{num_to_select}"
        f"T{tr[-1] + 1}_EulerianScalesHL"
        f"{'_MBB' if do_mbb else ''}.mat"
    )
    output_path = Path(
        cfg.get("output_path", input_dir / default_name)
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

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