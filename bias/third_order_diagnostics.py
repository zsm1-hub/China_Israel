"""Online third-order diagnostics with uniform reservoirs of delta-u_L."""
import numpy as np


def init_third_order(nscale, reservoir_size=100_000, random_seed=12345):
    return {
        "count": np.zeros(nscale, dtype=np.int64),
        "positive_count": np.zeros(nscale, dtype=np.int64),
        "negative_count": np.zeros(nscale, dtype=np.int64),
        "sum": np.zeros(nscale), "abs_sum": np.zeros(nscale), "sq_sum": np.zeros(nscale),
        "reservoir_size": int(reservoir_size),
        "dl_values": [np.empty(0) for _ in range(nscale)],
        "dl_keys": [np.empty(0) for _ in range(nscale)],
        "rng": np.random.default_rng(random_seed),
    }


def update_third_order(d, scale_bin, dl, dt):
    dl = np.asarray(dl, float); dt = np.asarray(dt, float); b = np.asarray(scale_bin, int)
    x = dl**3 + dl*dt**2
    good = np.isfinite(x) & np.isfinite(dl) & (b >= 0) & (b < d["count"].size)
    for k in np.unique(b[good]):
        mask = good & (b == k); v = x[mask]; longitudinal = dl[mask]
        d["count"][k] += v.size
        d["positive_count"][k] += np.count_nonzero(v > 0)
        d["negative_count"][k] += np.count_nonzero(v < 0)
        d["sum"][k] += v.sum(); d["abs_sum"][k] += np.abs(v).sum(); d["sq_sum"][k] += np.square(v).sum()
        keys = d["rng"].random(longitudinal.size)
        values = np.r_[d["dl_values"][k], longitudinal]
        all_keys = np.r_[d["dl_keys"][k], keys]
        cap = d["reservoir_size"]
        if values.size > cap:
            keep = np.argpartition(all_keys, -cap)[-cap:]
            values, all_keys = values[keep], all_keys[keep]
        d["dl_values"][k], d["dl_keys"][k] = values, all_keys


def finish_third_order(d):
    n = d["count"].astype(float)
    mean = np.divide(d["sum"], n, out=np.full_like(n, np.nan), where=n > 0)
    second = np.divide(d["sq_sum"], n, out=np.full_like(n, np.nan), where=n > 0)
    std = np.sqrt(np.maximum(second - mean**2, 0.0))
    se = np.divide(std, np.sqrt(n), out=np.full_like(n, np.nan), where=n > 1)
    count = np.array([len(x) for x in d["dl_values"]], dtype=np.int64)
    matrix = np.full((n.size, max(1, count.max(initial=0))), np.nan)
    for k, values in enumerate(d["dl_values"]): matrix[k, :values.size] = values
    return {
        "third_count": d["count"],
        "third_positive_count": d["positive_count"],
        "third_negative_count": d["negative_count"],
        "third_positive_fraction": np.divide(d["positive_count"], n, out=np.full_like(n, np.nan), where=n > 0),
        "third_mean": mean, "third_std": std, "third_standard_error": se,
        "third_relative_standard_error": np.divide(se, np.abs(mean), out=np.full_like(n, np.nan), where=np.abs(mean) > 0),
        "third_cancellation_ratio": np.divide(np.abs(d["sum"]), d["abs_sum"], out=np.full_like(n, np.nan), where=d["abs_sum"] > 0),
        "dl_reservoir": matrix, "dl_reservoir_count": count,
        "third_pdf_reservoir_size": d["reservoir_size"],
    }
