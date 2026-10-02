import numpy as np

from tqdm import tqdm
from pathlib import Path
import scipy.io as sio
from third_order_diagnostics import init_third_order, update_third_order, finish_third_order






def moving_time_block_bootstrap(
    time_sum,
    time_count,
    block_length,
    num_boot=1000,
    confidence=0.95,
    random_seed=12345,
):
    """Circular moving-time-block bootstrap of count-weighted means."""
    time_sum = np.asarray(time_sum, dtype=float)
    time_count = np.asarray(time_count, dtype=float)
    if time_sum.shape != time_count.shape or time_sum.ndim != 2:
        raise ValueError("time_sum and time_count require identical (ntime,nscale) shapes")
    ntime, nscale = time_sum.shape
    lengths = np.asarray(block_length, dtype=int).ravel()
    if lengths.size == 1:
        lengths = np.full(nscale, int(lengths[0]), dtype=int)
    if lengths.size != nscale:
        raise ValueError("block_length must be scalar or contain one value per scale")
    rng = np.random.default_rng(random_seed)
    estimates = np.full((nscale, int(num_boot)), np.nan)
    effective_blocks = np.zeros(nscale, dtype=int)
    for ir in range(nscale):
        length = int(np.clip(lengths[ir], 1, ntime))
        effective_blocks[ir] = ntime // length
        if effective_blocks[ir] < 2 or np.sum(time_count[:, ir]) <= 0:
            continue
        ndraw = int(np.ceil(ntime / length))
        offsets = np.arange(length)
        for ib in range(int(num_boot)):
            starts = rng.integers(0, ntime, size=ndraw)
            indices = ((starts[:, None] + offsets[None, :]) % ntime).ravel()[:ntime]
            denominator = np.sum(time_count[indices, ir])
            if denominator > 0:
                estimates[ir, ib] = np.sum(time_sum[indices, ir]) / denominator
    alpha = 0.5 * (1.0 - float(confidence))
    return {
        "D3_mbb": estimates,
        "D3_mbb_mean": np.nanmean(estimates, axis=1),
        "D3_mbb_stderr": np.nanstd(estimates, axis=1, ddof=1),
        "D3_mbb_ci_low": np.nanpercentile(estimates, 100.0 * alpha, axis=1),
        "D3_mbb_ci_high": np.nanpercentile(estimates, 100.0 * (1.0-alpha), axis=1),
        "D3_mbb_block_length_time": lengths,
        "D3_mbb_effective_blocks": effective_blocks,
        "D3_mbb_num_boot": int(num_boot),
        "D3_mbb_confidence": float(confidence),
    }


# ============================================================
# Calculate normalized spatial variability
# ============================================================
def calculate_H_from_local(
    local_values,
    min_valid_blocks
):
    """
    Calculate

        H = std /
            (abs(mean) + rms)

    across spatial blocks.
    """

    local_values = np.asarray(
        local_values,
        dtype=float
    )

    nscale = local_values.shape[0]

    H = np.full(nscale, np.nan)
    mean_value = np.full(nscale, np.nan)
    std_value = np.full(nscale, np.nan)
    rms_value = np.full(nscale, np.nan)
    cancellation_ratio = np.full(nscale, np.nan)

    nvalid_blocks = np.zeros(
        nscale,
        dtype=int
    )

    for ir in range(nscale):

        values = local_values[ir]
        values = values[
            np.isfinite(values)
        ]

        nvalid_blocks[ir] = values.size

        if values.size < min_valid_blocks:
            continue

        mean_value[ir] = np.mean(values)

        std_value[ir] = np.std(
            values,
            ddof=0
        )

        rms_value[ir] = np.sqrt(
            np.mean(values**2)
        )

        denominator = (
            np.abs(mean_value[ir])
            + rms_value[ir]
        )

        if denominator > 0:
            H[ir] = (
                std_value[ir]
                / denominator
            )

        if rms_value[ir] > 0:
            cancellation_ratio[ir] = (
                np.abs(mean_value[ir])
                / rms_value[ir]
            )

    return {
        "H": H,
        "mean": mean_value,
        "std": std_value,
        "rms": rms_value,
        "cancellation_ratio": cancellation_ratio,
        "nvalid_blocks": nvalid_blocks,
    }


# ============================================================
# Eulerian grid-block structure functions and H metrics
# ============================================================
def calculate_eulerian_homogeneity(
    u,
    v,
    r_values,
    dx,
    dy=None,
    nblock_y=4,
    nblock_x=4,
    periodic=False,
    directions=None,
    min_samples_per_block=500,
    min_valid_blocks=None,
    time_chunk=50,
    require_same_block=False,
    reservoir_size=100_000,
    random_seed=12345,
    do_moving_block_bootstrap=False,
    snapshot_dt=None,
    mbb_num_boot=1000,
    mbb_confidence=0.95,
    mbb_block_length_time=None,
    mbb_random_seed=271828,
    show_progress=True,
    progress_name="Eulerian homogeneity"
):
    """
    Calculate Eulerian grid-block first-, second-, and
    third-order structure functions.

    Calculated components
    ---------------------
    First order:
        Dl   = <du_L>
        Dt   = <du_T>

    Second order:
        Dll  = <du_L^2>
        Dtt  = <du_T^2>
        D2   = Dll + Dtt

    Third order:
        Dlll = <du_L^3>
        Dltt = <du_L du_T^2>
        D3   = Dlll + Dltt

    The H metric is calculated across fixed spatial blocks:

        H = std(local D) /
            (abs(mean(local D)) + rms(local D))

    Parameters
    ----------
    u, v : ndarray
        Shape (ntime, ny, nx) or (ny, nx).

    r_values : array
        Requested separations, using the same units as dx and dy.

    dx, dy : float
        Grid spacing.

    periodic : bool
        True for HIT; False for Iceland.

    time_chunk : int
        Number of time records processed simultaneously.

    require_same_block : bool
        False:
            assign each increment to the block containing
            its reference point.

        True:
            retain only increments whose two endpoints are
            within the same spatial block.

    Returns
    -------
    result : dict
        Flat dictionary containing all local moments,
        H metrics, means, standard deviations, RMS values,
        cancellation ratios, and sampling information.
    """

    # --------------------------------------------------------
    # Input checks
    # --------------------------------------------------------
    u = np.asarray(u)
    v = np.asarray(v)

    if u.shape != v.shape:
        raise ValueError(
            "u and v must have the same shape."
        )

    if u.ndim == 2:
        u = u[None, :, :]
        v = v[None, :, :]

    if u.ndim != 3:
        raise ValueError(
            "u and v must have shape "
            "(ntime, ny, nx) or (ny, nx)."
        )

    if dy is None:
        dy = dx

    r_values = np.asarray(
        r_values,
        dtype=float
    ).ravel()

    if np.any(r_values <= 0):
        raise ValueError(
            "All r_values must be positive."
        )

    ntime, ny, nx = u.shape
    nscale = len(r_values)
    nblock = nblock_y * nblock_x

    third_diag = init_third_order(
        nscale,
        reservoir_size=int(reservoir_size),
        random_seed=int(random_seed),
    )

    if do_moving_block_bootstrap:
        time_D2_sum = np.zeros((ntime, nscale), dtype=np.float64)
        time_D3_sum = np.zeros((ntime, nscale), dtype=np.float64)
        time_pair_count = np.zeros((ntime, nscale), dtype=np.int64)

    if min_valid_blocks is None:
        min_valid_blocks = max(
            2,
            nblock // 2
        )

    if directions is None:
        directions = [
            (0, 1),
            (1, 1),
            (1, 0),
            (1, -1),
            (0, -1),
            (-1, -1),
            (-1, 0),
            (-1, 1),
        ]

    # Structure-function components
    component_names = [
        "Dl",
        "Dt",
        "Dll",
        "Dtt",
        "D2",
        "Dlll",
        "Dltt",
        "D3",
    ]

    print(
        f"{progress_name}\n"
        f"  velocity shape: {u.shape}\n"
        f"  requested scales: {nscale}\n"
        f"  directions: {len(directions)}\n"
        f"  spatial blocks: "
        f"{nblock_y} x {nblock_x}\n"
        f"  time chunk: {time_chunk}\n"
        f"  periodic: {periodic}\n"
        f"  require same block: {require_same_block}"
    )

    # --------------------------------------------------------
    # Fixed spatial-block map
    # --------------------------------------------------------
    y_block = np.minimum(
        np.arange(ny) * nblock_y // ny,
        nblock_y - 1
    )

    x_block = np.minimum(
        np.arange(nx) * nblock_x // nx,
        nblock_x - 1
    )

    block_id = (
        y_block[:, None] * nblock_x
        + x_block[None, :]
    ).astype(np.int64)

    # --------------------------------------------------------
    # Storage
    # --------------------------------------------------------
    local = {
        name: np.full(
            (nscale, nblock),
            np.nan,
            dtype=float
        )
        for name in component_names
    }

    sample_count = np.zeros(
        (nscale, nblock),
        dtype=np.int64
    )

    r_actual = np.full(
        nscale,
        np.nan,
        dtype=float
    )

    # --------------------------------------------------------
    # Nonperiodic slicing helper
    # --------------------------------------------------------
    def paired_slices(size, shift):

        if shift > 0:
            reference = slice(0, size - shift)
            displaced = slice(shift, size)

        elif shift < 0:
            reference = slice(-shift, size)
            displaced = slice(0, size + shift)

        else:
            reference = slice(0, size)
            displaced = slice(0, size)

        return reference, displaced

    # --------------------------------------------------------
    # Progress setup
    # --------------------------------------------------------
    nchunks = int(
        np.ceil(ntime / time_chunk)
    )

    total_steps = (
        nscale
        * len(directions)
        * nchunks
    )

    if show_progress:
        progress = tqdm(
            total=total_steps,
            desc=progress_name,
            unit="chunk"
        )
    else:
        progress = None

    # ========================================================
    # Separation loop
    # ========================================================
    for ir, target_r in enumerate(r_values):

        # Accumulated sum for each quantity and block
        sum_block = {
            name: np.zeros(
                nblock,
                dtype=np.float64
            )
            for name in component_names
        }

        count_block = np.zeros(
            nblock,
            dtype=np.int64
        )

        actual_separations = []

        # ====================================================
        # Direction loop
        # ====================================================
        for base_di, base_dj in directions:

            base_distance = np.sqrt(
                (base_di * dy)**2
                + (base_dj * dx)**2
            )

            if base_distance == 0:
                if progress is not None:
                    progress.update(nchunks)
                continue

            multiplier = max(
                1,
                int(np.round(
                    target_r / base_distance
                ))
            )

            di = int(base_di * multiplier)
            dj = int(base_dj * multiplier)

            if periodic:
                # Only minimum-image displacements are distinct.
                if abs(di) > ny // 2 or abs(dj) > nx // 2:
                    if progress is not None:
                        progress.update(nchunks)
                    continue
            else:
                if abs(di) >= ny or abs(dj) >= nx:
                    if progress is not None:
                        progress.update(nchunks)
                    continue

            separation = np.sqrt(
                (di * dy)**2
                + (dj * dx)**2
            )

            if separation == 0:
                if progress is not None:
                    progress.update(nchunks)
                continue

            actual_separations.append(
                separation
            )

            # Unit vectors
            ex = dj * dx / separation
            ey = di * dy / separation

            # Transverse direction:
            # e_T = (-ey, ex)
            tx = -ey
            ty = ex

            # ------------------------------------------------
            # Spatial layout for this direction
            # ------------------------------------------------
            if periodic:

                block_reference = block_id
                block_displaced = np.roll(
                    block_id, shift=(-di, -dj), axis=(0, 1)
                )

                # Assign the increment to its periodic midpoint block,
                # matching the particle-pair implementation.
                grid_y, grid_x = np.indices((ny, nx), dtype=float)
                midpoint_y = np.mod(grid_y + 0.5 * di, ny)
                midpoint_x = np.mod(grid_x + 0.5 * dj, nx)
                midpoint_by = np.minimum(
                    (midpoint_y * nblock_y // ny).astype(int), nblock_y - 1
                )
                midpoint_bx = np.minimum(
                    (midpoint_x * nblock_x // nx).astype(int), nblock_x - 1
                )
                increment_block = midpoint_bx * nblock_y + midpoint_by

            else:

                y0, y1 = paired_slices(
                    ny,
                    di
                )

                x0, x1 = paired_slices(
                    nx,
                    dj
                )

                block_reference = block_id[
                    y0,
                    x0
                ]

                block_displaced = block_id[
                    y1,
                    x1
                ]
                increment_block = block_reference

            if require_same_block:
                same_block_spatial = (
                    block_reference
                    == block_displaced
                )
            else:
                same_block_spatial = np.ones(
                    block_reference.shape,
                    dtype=bool
                )

            spatial_blocks = (
                increment_block.ravel()
            )

            spatial_same_block = (
                same_block_spatial.ravel()
            )

            # =================================================
            # Time-chunk loop
            # =================================================
            for t0 in range(
                0,
                ntime,
                time_chunk
            ):
                t1 = min(
                    t0 + time_chunk,
                    ntime
                )

                u_chunk = u[t0:t1]
                v_chunk = v[t0:t1]

                if periodic:

                    u_displaced = np.roll(
                        u_chunk,
                        shift=(-di, -dj),
                        axis=(1, 2)
                    )

                    v_displaced = np.roll(
                        v_chunk,
                        shift=(-di, -dj),
                        axis=(1, 2)
                    )

                    du = (
                        u_displaced
                        - u_chunk
                    )

                    dv = (
                        v_displaced
                        - v_chunk
                    )

                else:

                    u_reference = u_chunk[
                        :,
                        y0,
                        x0
                    ]

                    u_displaced = u_chunk[
                        :,
                        y1,
                        x1
                    ]

                    v_reference = v_chunk[
                        :,
                        y0,
                        x0
                    ]

                    v_displaced = v_chunk[
                        :,
                        y1,
                        x1
                    ]

                    du = (
                        u_displaced
                        - u_reference
                    )

                    dv = (
                        v_displaced
                        - v_reference
                    )

                # --------------------------------------------
                # Longitudinal and transverse increments
                # --------------------------------------------
                du_L = (
                    du * ex
                    + dv * ey
                )

                du_T = (
                    du * tx
                    + dv * ty
                )

                if do_moving_block_bootstrap:
                    time_valid = (
                        np.isfinite(du_L)
                        & np.isfinite(du_T)
                        & same_block_spatial[None, ...]
                    )
                    for jt in range(t1 - t0):
                        good_t = time_valid[jt]
                        if np.any(good_t):
                            Lt = du_L[jt][good_t]
                            Tt = du_T[jt][good_t]
                            time_pair_count[t0 + jt, ir] += Lt.size
                            time_D2_sum[t0 + jt, ir] += np.sum(Lt**2 + Tt**2)
                            time_D3_sum[t0 + jt, ir] += np.sum(Lt**3 + Lt*Tt**2)

                # Flatten once
                du_L_flat = du_L.ravel()
                du_T_flat = du_T.ravel()

                chunk_length = t1 - t0

                blocks_flat = np.tile(
                    spatial_blocks,
                    chunk_length
                )

                same_block_flat = np.tile(
                    spatial_same_block,
                    chunk_length
                )

                valid = (
                    np.isfinite(du_L_flat)
                    & np.isfinite(du_T_flat)
                    & same_block_flat
                )

                if np.any(valid):

                    L = du_L_flat[valid]
                    T = du_T_flat[valid]
                    B = blocks_flat[valid]

                    update_third_order(
                        third_diag,
                        np.full(L.size, ir, dtype=np.int64),
                        L,
                        T,
                    )

                    # Powers reused by several components
                    L2 = L * L
                    T2 = T * T
                    L3 = L2 * L
                    LTT = L * T2

                    # Pair count
                    count_block += np.bincount(
                        B,
                        minlength=nblock
                    ).astype(np.int64)

                    # First order
                    sum_block["Dl"] += np.bincount(
                        B,
                        weights=L,
                        minlength=nblock
                    )

                    sum_block["Dt"] += np.bincount(
                        B,
                        weights=T,
                        minlength=nblock
                    )

                    # Second order
                    sum_block["Dll"] += np.bincount(
                        B,
                        weights=L2,
                        minlength=nblock
                    )

                    sum_block["Dtt"] += np.bincount(
                        B,
                        weights=T2,
                        minlength=nblock
                    )

                    sum_block["D2"] += np.bincount(
                        B,
                        weights=L2 + T2,
                        minlength=nblock
                    )

                    # Third order
                    sum_block["Dlll"] += np.bincount(
                        B,
                        weights=L3,
                        minlength=nblock
                    )

                    sum_block["Dltt"] += np.bincount(
                        B,
                        weights=LTT,
                        minlength=nblock
                    )

                    sum_block["D3"] += np.bincount(
                        B,
                        weights=L3 + LTT,
                        minlength=nblock
                    )

                    del L, T, B
                    del L2, T2, L3, LTT

                del du, dv
                del du_L, du_T
                del du_L_flat, du_T_flat
                del blocks_flat, same_block_flat
                del valid

                if periodic:
                    del u_displaced
                    del v_displaced
                else:
                    del u_reference
                    del u_displaced
                    del v_reference
                    del v_displaced

                if progress is not None:
                    progress.set_postfix_str(
                        f"r={target_r:.3g}, "
                        f"dir=({base_di},{base_dj})"
                    )
                    progress.update(1)

        # ----------------------------------------------------
        # Convert sums to block means
        # ----------------------------------------------------
        valid_block = (
            count_block
            >= min_samples_per_block
        )

        for name in component_names:

            local[name][
                ir,
                valid_block
            ] = (
                sum_block[name][valid_block]
                / count_block[valid_block]
            )

        sample_count[ir] = count_block

        if actual_separations:
            r_actual[ir] = np.mean(
                actual_separations
            )

    if progress is not None:
        progress.close()

    # ========================================================
    # Calculate H for each quantity
    # ========================================================
    diagnostic = {}

    for name in component_names:

        diagnostic[name] = (
            calculate_H_from_local(
                local[name],
                min_valid_blocks
            )
        )

    # ========================================================
    # Flat output dictionary
    # ========================================================
    result = {
        "r_requested": r_values,
        "r_actual": r_actual,
        "sample_count": sample_count,
        "nblock_y": nblock_y,
        "nblock_x": nblock_x,
        "min_samples_per_block":
            min_samples_per_block,
        "min_valid_blocks":
            min_valid_blocks,
        "periodic": periodic,
        "require_same_block":
            require_same_block,
        "time_chunk": time_chunk,
    }

    # Save all local quantities and diagnostics
    for name in component_names:

        result[f"local_{name}"] = (
            local[name]
        )

        result[f"H_{name}"] = (
            diagnostic[name]["H"]
        )

        result[f"mean_{name}"] = (
            diagnostic[name]["mean"]
        )

        result[f"std_{name}"] = (
            diagnostic[name]["std"]
        )

        result[f"rms_{name}"] = (
            diagnostic[name]["rms"]
        )

        result[
            f"cancellation_ratio_{name}"
        ] = diagnostic[name][
            "cancellation_ratio"
        ]

        result[
            f"nvalid_blocks_{name}"
        ] = diagnostic[name][
            "nvalid_blocks"
        ]

    # Convenient aliases matching manuscript notation
    result["H1L"] = result["H_Dl"]
    result["H1T"] = result["H_Dt"]

    result["H2L"] = result["H_Dll"]
    result["H2T"] = result["H_Dtt"]
    result["H2"] = result["H_D2"]

    result["H3L"] = result["H_Dlll"]
    result["H3LTT"] = result["H_Dltt"]
    result["H3"] = result["H_D3"]

    print(
        f"\nCompleted: {progress_name}"
    )

    for name in [
        "Dl",
        "Dt",
        "Dll",
        "Dtt",
        "D2",
        "Dlll",
        "Dltt",
        "D3",
    ]:
        nvalid = np.count_nonzero(
            np.isfinite(
                result[f"H_{name}"]
            )
        )

        print(
            f"  H_{name}: "
            f"{nvalid}/{nscale} valid scales"
        )

    # Particle-compatible aliases used by the six-column plot.
    result["dist_axis"] = result["r_actual"]
    result["pair_count"] = result["sample_count"]
    result["count_total"] = np.sum(result["sample_count"], axis=1)
    # Pair-count-weighted global moments, matching the particle code.
    for name in component_names:
        local_value = result[f"local_{name}"]
        numerator = np.nansum(local_value * result["sample_count"], axis=1)
        denominator = np.sum(
            np.where(np.isfinite(local_value), result["sample_count"], 0), axis=1
        )
        result[name] = np.divide(
            numerator, denominator, out=np.full(nscale, np.nan), where=denominator > 0
        )
    result["do_moving_block_bootstrap"] = bool(do_moving_block_bootstrap)
    if do_moving_block_bootstrap:
        if mbb_block_length_time is None:
            if snapshot_dt is None or float(snapshot_dt) <= 0:
                raise ValueError(
                    "snapshot_dt must be positive when automatic MBB block length is used"
                )
            tau_r = np.divide(
                result["dist_axis"], np.sqrt(result["D2"]),
                out=np.full(nscale, np.nan),
                where=np.isfinite(result["D2"]) & (result["D2"] > 0),
            )
            lengths = np.ceil(tau_r / float(snapshot_dt))
            lengths[~np.isfinite(lengths)] = ntime
            lengths = np.maximum(lengths, 1).astype(int)
        else:
            tau_r = np.divide(
                result["dist_axis"], np.sqrt(result["D2"]),
                out=np.full(nscale, np.nan),
                where=np.isfinite(result["D2"]) & (result["D2"] > 0),
            )
            supplied = np.asarray(mbb_block_length_time, dtype=int).ravel()
            lengths = np.full(nscale, int(supplied[0])) if supplied.size == 1 else supplied
        mbb = moving_time_block_bootstrap(
            time_D3_sum, time_pair_count, lengths,
            num_boot=mbb_num_boot, confidence=mbb_confidence,
            random_seed=mbb_random_seed,
        )
        result.update(mbb)
        result.update({
            "D3_mbb_tau_r": tau_r,
            "D3_time_sum": time_D3_sum,
            "D2_time_sum": time_D2_sum,
            "pair_count_time": time_pair_count,
            "snapshot_dt": float(snapshot_dt) if snapshot_dt is not None else np.nan,
        })

    result.update(finish_third_order(third_diag))

    return result


# ============================================================
# Save result as MATLAB file
# ============================================================
def save_homogeneity_mat(
    result,
    output_path,
    case_name,
    dx,
    dy,
    length_unit
):
    """
    Save all structure-function components and H metrics
    as a MATLAB .mat file.
    """

    output_path = Path(output_path)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    mat_data = {}

    for key, value in result.items():

        if isinstance(
            value,
            (np.ndarray, list, tuple)
        ):
            mat_data[key] = np.asarray(value)

        elif isinstance(
            value,
            (int, float, bool, np.number)
        ):
            mat_data[key] = np.asarray(value)

        elif isinstance(value, str):
            mat_data[key] = value

    mat_data["case_name"] = case_name
    mat_data["dx"] = np.asarray(dx)
    mat_data["dy"] = np.asarray(dy)
    mat_data["length_unit"] = length_unit

    if result["require_same_block"]:
        mat_data["block_definition"] = (
            "both endpoints required to lie "
            "within the same spatial block"
        )
    else:
        mat_data["block_definition"] = (
            "increment assigned to the block "
            "containing its reference point"
        )

    mat_data["definitions"] = (
        "Dl=<duL>; Dt=<duT>; "
        "Dll=<duL^2>; Dtt=<duT^2>; "
        "D2=Dll+Dtt; "
        "Dlll=<duL^3>; "
        "Dltt=<duL*duT^2>; "
        "D3=Dlll+Dltt"
    )

    sio.savemat(
        output_path,
        mat_data,
        do_compression=True
    )

    print(
        f"\nSaved MATLAB file:\n"
        f"{output_path}"
    )


def load_hit_velocity_fields(grid_dir, timerange_matlab, file_offset=4000, show_progress=True):
    """Reconstruct fixed-grid HIT velocities from spectral MAT files."""
    grid_dir = Path(grid_dir)
    times = np.asarray(timerange_matlab, dtype=int).ravel()
    if times.size == 0:
        raise ValueError("timerange_matlab cannot be empty.")
    parameters = sio.loadmat(grid_dir / "HIT2D_Parameters.mat")
    diag, mm, kk = parameters["Diag"], parameters["mm"], parameters["kk"]
    numbers = int(file_offset) + times - 1
    u = np.empty((numbers.size,) + diag.shape, dtype=np.float64)
    v = np.empty_like(u)
    iterator = tqdm(numbers, desc="Load HIT Eulerian fields", unit="field", disable=not show_progress)
    for it, number in enumerate(iterator):
        hq = sio.loadmat(grid_dir / f"HIT2D_t_{number}.mat")["hq"]
        hpsi = diag * hq
        u[it] = np.real(np.fft.ifft2(-1j * mm * hpsi))
        v[it] = np.real(np.fft.ifft2(1j * kk * hpsi))
    return u, v, numbers


def run_hit_eulerian(cfg):
    """Load HIT fields, calculate diagnostics, optionally save, and return result/path."""
    u, v, numbers = load_hit_velocity_fields(
        cfg["grid_dir"], cfg["timerange_matlab"],
        cfg.get("file_offset", 4000), cfg.get("show_progress", True)
    )
    result = calculate_eulerian_homogeneity(
        u=u, v=v, r_values=cfg["r_requested"],
        dx=float(cfg.get("dx", 2.0 * np.pi / u.shape[-1])), dy=cfg.get("dy"),
        nblock_y=int(cfg.get("nblock_y", 4)), nblock_x=int(cfg.get("nblock_x", 4)),
        periodic=True, directions=cfg.get("directions"),
        min_samples_per_block=int(cfg.get("min_samples_per_block", 500)),
        min_valid_blocks=int(cfg.get("min_valid_blocks", 8)),
        time_chunk=int(cfg.get("time_chunk", 10)),
        require_same_block=bool(cfg.get("require_same_block", False)),
        reservoir_size=int(cfg.get("third_pdf_reservoir_size", 100_000)),
        random_seed=int(cfg.get("random_seed", 12345)),
        do_moving_block_bootstrap=bool(cfg.get("do_moving_block_bootstrap", False)),
        snapshot_dt=cfg.get("snapshot_dt"),
        mbb_num_boot=int(cfg.get("mbb_num_boot", 1000)),
        mbb_confidence=float(cfg.get("mbb_confidence", 0.95)),
        mbb_block_length_time=cfg.get("mbb_block_length_time"),
        mbb_random_seed=int(cfg.get("mbb_random_seed", 271828)),
        show_progress=bool(cfg.get("show_progress", True)),
        progress_name="HIT Eulerian reference",
    )
    result["timerange"] = np.asarray(cfg["timerange_matlab"], dtype=int)
    result["source_file_numbers"] = numbers
    result["Case"] = "HIT2d Eulerian"
    output_path = cfg.get("output_path")
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        sio.savemat(output_path, result, do_compression=True, oned_as="row")
        print(f"Saved Eulerian result:\n{output_path}")
    return result, output_path



