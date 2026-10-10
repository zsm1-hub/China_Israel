"""Online third-order diagnostics with paired uniform reservoirs."""
import numpy as np


def init_third_order(nscale, reservoir_size=100_000, random_seed=12345):
    nscale = int(nscale)
    return {
        "count": np.zeros(nscale, dtype=np.int64),
        "positive_count": np.zeros(nscale, dtype=np.int64),
        "negative_count": np.zeros(nscale, dtype=np.int64),
        "sum": np.zeros(nscale), "abs_sum": np.zeros(nscale), "sq_sum": np.zeros(nscale),
        "reservoir_size": int(reservoir_size),
        "dl_values": [np.empty(0) for _ in range(nscale)],
        "dtr_values": [np.empty(0) for _ in range(nscale)],
        "reservoir_keys": [np.empty(0) for _ in range(nscale)],
        "rng": np.random.default_rng(random_seed),
    }


def update_third_order(d, scale_bin, dl, dt):
    """Update full D3 statistics and a pair-aligned (dl, dtr) reservoir."""
    dl = np.asarray(dl, float).ravel()
    dt = np.asarray(dt, float).ravel()
    b = np.asarray(scale_bin, int).ravel()
    if not (dl.size == dt.size == b.size):
        raise ValueError("scale_bin, dl, and dt must have equal lengths")
    x = dl**3 + dl*dt**2
    good = np.isfinite(x) & np.isfinite(dl) & np.isfinite(dt) & (b >= 0) & (b < d["count"].size)
    for k in np.unique(b[good]):
        q = good & (b == k)
        longitudinal, transverse, values = dl[q], dt[q], x[q]
        d["count"][k] += values.size
        d["positive_count"][k] += np.count_nonzero(values > 0)
        d["negative_count"][k] += np.count_nonzero(values < 0)
        d["sum"][k] += values.sum()
        d["abs_sum"][k] += np.abs(values).sum()
        d["sq_sum"][k] += np.square(values).sum()
        keys = d["rng"].random(values.size)
        all_dl = np.r_[d["dl_values"][k], longitudinal]
        all_dtr = np.r_[d["dtr_values"][k], transverse]
        all_keys = np.r_[d["reservoir_keys"][k], keys]
        cap = d["reservoir_size"]
        if all_keys.size > cap:
            keep = np.argpartition(all_keys, -cap)[-cap:]
            all_dl, all_dtr, all_keys = all_dl[keep], all_dtr[keep], all_keys[keep]
        d["dl_values"][k], d["dtr_values"][k], d["reservoir_keys"][k] = all_dl, all_dtr, all_keys


def finish_third_order(d):
    """Return all-sample statistics and fixed-size paired reservoirs."""
    n = d["count"].astype(float)
    mean = np.divide(d["sum"], n, out=np.full_like(n, np.nan), where=n > 0)
    second = np.divide(d["sq_sum"], n, out=np.full_like(n, np.nan), where=n > 0)
    std = np.sqrt(np.maximum(second - mean**2, 0.0))
    se = np.divide(std, np.sqrt(n), out=np.full_like(n, np.nan), where=n > 0)
    nscale, cap = d["count"].size, d["reservoir_size"]
    dl_matrix = np.full((nscale, cap), np.nan)
    dtr_matrix = np.full((nscale, cap), np.nan)
    count = np.zeros(nscale, dtype=np.int64)
    for k in range(nscale):
        number = min(d["dl_values"][k].size, d["dtr_values"][k].size, cap)
        if number:
            dl_matrix[k, :number] = d["dl_values"][k][:number]
            dtr_matrix[k, :number] = d["dtr_values"][k][:number]
            count[k] = number
    dl3 = dl_matrix**3
    dldtr2 = dl_matrix*dtr_matrix**2
    d3 = dl3 + dldtr2
    return {
        "third_count": d["count"],
        "third_positive_count": d["positive_count"],
        "third_negative_count": d["negative_count"],
        "third_positive_fraction": np.divide(d["positive_count"], n, out=np.full_like(n, np.nan), where=n > 0),
        "third_mean": mean, "third_std": std, "third_standard_error": se,
        "third_relative_standard_error": np.divide(se, np.abs(mean), out=np.full_like(n, np.nan), where=np.abs(mean) > 0),
        "third_cancellation_ratio": np.divide(np.abs(d["sum"]), d["abs_sum"], out=np.full_like(n, np.nan), where=d["abs_sum"] > 0),
        "dl_reservoir": dl_matrix, "dtr_reservoir": dtr_matrix,
        "dl3_reservoir": dl3, "dldtr2_reservoir": dldtr2, "d3_reservoir": d3,
        "dl_reservoir_count": count, "dtr_reservoir_count": count.copy(), "d3_reservoir_count": count.copy(),
        "third_pdf_reservoir_size": cap,
    }
