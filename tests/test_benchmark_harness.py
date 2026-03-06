import copy
import numpy as np
import pytest

from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data_with_additive


def _make_config_with_dt(base_config: dict, dt: int) -> dict:
    cfg = copy.deepcopy(base_config)
    total_T = int(cfg["n_timesteps"] * cfg.get("dt", 1))
    cfg["dt"] = int(dt)
    cfg["n_timesteps"] = int(total_T // dt)
    if cfg["n_timesteps"] < 5:
        cfg["n_timesteps"] = 5
    cfg["solver_t_eval"] = "timesteps"
    return cfg


def _make_data_release(base_data: dict, model: str, params: dict) -> dict:
    d = copy.deepcopy(base_data)
    d["additive_release"] = {"model": model, "params": dict(params)}
    return d


def _total_additive(out) -> np.ndarray:
    return out.A_part.sum(axis=0) + out.A_aq


def _released_fraction(out) -> np.ndarray:
    A0 = float(out.A_part[:, 0].sum() + out.A_aq[0])
    return out.A_aq / max(A0, 1e-30)


def _t50(t: np.ndarray, rel_frac: np.ndarray) -> float:
    if np.all(rel_frac < 0.5):
        return float("nan")
    idx = np.where(rel_frac >= 0.5)[0][0]
    if idx == 0:
        return float(t[0])
    t0, t1 = t[idx - 1], t[idx]
    y0, y1 = rel_frac[idx - 1], rel_frac[idx]
    if y1 == y0:
        return float(t1)
    return float(t0 + (0.5 - y0) * (t1 - t0) / (y1 - y0))


def _peak_flux_time(t: np.ndarray, A_aq: np.ndarray) -> float:
    dt = np.diff(t)
    dA = np.diff(A_aq)
    flux = dA / np.maximum(dt, 1e-30)
    if flux.size == 0:
        return float("nan")
    return float(t[1:][np.argmax(flux)])


def _tols_from_bi(Bi: float):
    """
    Option B: tolerances reflect the *approximate* analytical model regimes.
    """
    if Bi <= 1.0:
        # external-transfer dominated: should match pretty well
        return dict(ts=0.05, t50=0.10, tpeak=0.20)
    if Bi < 100.0:
        # bridge approximation: allow a bit more error
        return dict(ts=0.10, t50=0.15, tpeak=0.25)
    # Bi >= 100: analytical uses Dirichlet series (perfect sink), so expect
    # larger discrepancy
    return dict(ts=0.25, t50=0.50, tpeak=0.60)


@pytest.mark.parametrize(
    "D_p,K_pw,stress_kfrag",
    [
        (1e-18, 1e3, 0.005),
        (1e-16, 1e4, 0.01),
        (1e-14, 1e5, 0.02),
    ],
)
def test_benchmark_analytical_vs_numerical_agree_and_conserve(D_p, K_pw, stress_kfrag):
    cfg = copy.deepcopy(minimal_config)

    data_base = copy.deepcopy(minimal_data_with_additive)
    data_base["k_frag"] = float(stress_kfrag)

    N = cfg["n_size_classes"]
    data_base["initial_additive_concs"] = [1.0] * N

    # keep consistent with your analytical Bi definition
    D_w = 1e-9
    common = {"D_p": D_p, "D_w": D_w, "K_pw": K_pw}

    # Higher resolution in "verification mode" (still affordable)
    data_a = _make_data_release(data_base, "analytical", {**common, "n_terms": 200})
    out_a = FragmentMNP(cfg, data_a).run()

    data_n = _make_data_release(data_base, "numerical", {**common, "n_r": 200, "n_substeps": 80})
    out_n = FragmentMNP(cfg, data_n).run()

    # 1) mass conservation
    A_tot_a = _total_additive(out_a)
    A_tot_n = _total_additive(out_n)
    A0 = float(A_tot_a[0])

    assert np.max(np.abs(A_tot_a - A0) / max(A0, 1e-30)) < 1e-10
    assert np.max(np.abs(A_tot_n - A0) / max(A0, 1e-30)) < 1e-10

    # 2) agreement (robust error metric) with regime-dependent tolerance
    ra = _released_fraction(out_a)
    rn = _released_fraction(out_n)

    den = np.maximum(np.maximum(np.abs(ra), np.abs(rn)), 1e-12)
    sym_rel_err = np.abs(ra - rn) / den

    # Ignore near-zero regime where relative error is meaningless
    mask = (ra > 1e-6) | (rn > 1e-6)
    max_err = float(np.nanmax(sym_rel_err[mask])) if np.any(mask) else 0.0

    Bi = (K_pw * D_w) / D_p
    tols = _tols_from_bi(Bi)

    assert max_err < tols["ts"]

    # 3) key summary metrics (also regime-dependent)
    t50_a, t50_n = _t50(out_a.t, ra), _t50(out_n.t, rn)
    if np.isfinite(t50_a) and np.isfinite(t50_n):
        assert abs(t50_a - t50_n) / max(t50_a, 1e-30) < tols["t50"]

    tp_a, tp_n = _peak_flux_time(out_a.t, out_a.A_aq), _peak_flux_time(out_n.t, out_n.A_aq)
    if np.isfinite(tp_a) and np.isfinite(tp_n):
        assert abs(tp_a - tp_n) / max(tp_a, 1e-30) < tols["tpeak"]


def test_numerical_converges_with_dt_refinement():
    cfg0 = copy.deepcopy(minimal_config)
    cfg0["dt"] = 200
    cfg0["n_timesteps"] = 50  # total duration = 10000 s
    cfg0["solver_t_eval"] = "timesteps"

    data = copy.deepcopy(minimal_data_with_additive)
    data["k_frag"] = 0.02
    N = cfg0["n_size_classes"]
    data["initial_additive_concs"] = [1.0] * N
    data["additive_release"] = {
        "model": "numerical",
        "params": {
            "D_p": 1e-16,
            "D_w": 1e-9,
            "K_pw": 1e4,
            "n_r": 200,
            "n_substeps": 80,
        },
    }

    dts = [200, 100, 50]
    finals = []
    for dt in dts:
        cfg = _make_config_with_dt(cfg0, dt=dt)
        out = FragmentMNP(cfg, data).run()
        finals.append(float(_released_fraction(out)[-1]))

    # Check for numerical blow-ups first
    assert np.all(np.isfinite(finals)), f"Non-finite final released fractions: {finals}"

    err1 = abs(finals[0] - finals[1])
    err2 = abs(finals[1] - finals[2])

    # Expect refinement to reduce the change, allowing for a small noise floor
    noise_floor = 1e-10
    assert (err2 + noise_floor) < (err1 + noise_floor)