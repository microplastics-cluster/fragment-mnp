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


def _make_config_release(base_config: dict, model: str, solver: dict | None = None) -> dict:
    cfg = copy.deepcopy(base_config)
    cfg["additive_release"] = {
        "model": model,
        "solver": {} if solver is None else dict(solver),
    }
    return cfg


def _make_data_release(base_data: dict, physical_params: dict) -> dict:
    d = copy.deepcopy(base_data)
    d["additive_release"] = dict(physical_params)
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


def _peak_flux_value(t: np.ndarray, A_aq: np.ndarray) -> float:
    dt = np.diff(t)
    dA = np.diff(A_aq)
    flux = dA / np.maximum(dt, 1e-30)
    if flux.size == 0:
        return float("nan")
    return float(np.max(flux))


def _classify_bi_regime(Bi: float) -> str:
    if Bi <= 1.0:
        return "low Bi: external transfer limited"
    if Bi < 100.0:
        return "intermediate Bi"
    return "high Bi: internal diffusion limited"


def _tols_from_bi(Bi: float):
    """
    Tolerances reflect the approximate analytical-model regime.
    """
    if Bi <= 1.0:
        return dict(ts=0.05, t50=0.10, tpeak=0.20)
    if Bi < 100.0:
        return dict(ts=0.10, t50=0.15, tpeak=0.25)
    return dict(ts=0.25, t50=0.50, tpeak=0.60)


def _run_case(D_p: float, K_pw: float, stress_kfrag: float, D_w: float = 1e-9):
    cfg = copy.deepcopy(minimal_config)

    data_base = copy.deepcopy(minimal_data_with_additive)
    data_base["k_frag"] = float(stress_kfrag)
    data_base["initial_additive_concs"] = [1.0] * cfg["n_size_classes"]

    common = {"D_p": D_p, "D_w": D_w, "K_pw": K_pw}

    cfg_a = _make_config_release(cfg, "analytical", {"n_terms": 200})
    data_a = _make_data_release(data_base, common)
    out_a = FragmentMNP(cfg_a, data_a).run()

    cfg_n = _make_config_release(
        cfg,
        "numerical",
        {
            "n_r": 200,
            "n_substeps": 80,
            "theta": 1.0,
        },
    )
    data_n = _make_data_release(data_base, common)
    out_n = FragmentMNP(cfg_n, data_n).run()

    A_tot_a = _total_additive(out_a)
    A_tot_n = _total_additive(out_n)
    A0 = float(A_tot_a[0])

    ra = _released_fraction(out_a)
    rn = _released_fraction(out_n)
    den = np.maximum(np.maximum(np.abs(ra), np.abs(rn)), 1e-12)
    sym_rel_err = np.abs(ra - rn) / den
    mask = (ra > 1e-6) | (rn > 1e-6)
    max_err = float(np.nanmax(sym_rel_err[mask])) if np.any(mask) else 0.0
    rmse = float(np.sqrt(np.mean((ra - rn) ** 2)))

    Bi = (K_pw * D_w) / D_p

    return {
        "Bi": float(Bi),
        "regime": _classify_bi_regime(Bi),
        "max_sym_rel_err": max_err,
        "rmse_released_fraction": rmse,
        "t50_analytical": _t50(out_a.t, ra),
        "t50_numerical": _t50(out_n.t, rn),
        "t_peak_analytical": _peak_flux_time(out_a.t, out_a.A_aq),
        "t_peak_numerical": _peak_flux_time(out_n.t, out_n.A_aq),
        "peak_flux_analytical": _peak_flux_value(out_a.t, out_a.A_aq),
        "peak_flux_numerical": _peak_flux_value(out_n.t, out_n.A_aq),
        "released_final_analytical": float(ra[-1]),
        "released_final_numerical": float(rn[-1]),
        "mass_balance_relerr_analytical": float(np.max(np.abs(A_tot_a - A0) / max(A0, 1e-30))),
        "mass_balance_relerr_numerical": float(np.max(np.abs(A_tot_n - A0) / max(A0, 1e-30))),
    }


REGIME_CASES = [
    pytest.param(
        1e-12, 1.0, 0.005,
        id="low-bi_external-transfer-limited"
    ),
    pytest.param(
        1e-12, 1e2, 0.01,
        id="intermediate-bi"
    ),
    pytest.param(
        1e-16, 1e10, 0.02,
        id="high-bi_internal-diffusion-limited"
    ),
]


@pytest.mark.parametrize("D_p,K_pw,stress_kfrag", REGIME_CASES)
def test_benchmark_analytical_vs_numerical_by_release_regime(D_p, K_pw, stress_kfrag):
    metrics = _run_case(D_p=D_p, K_pw=K_pw, stress_kfrag=stress_kfrag)

    Bi = metrics["Bi"]
    regime = metrics["regime"]
    tols = _tols_from_bi(Bi)

    assert regime == _classify_bi_regime(Bi)

    assert metrics["mass_balance_relerr_analytical"] < 1e-10
    assert metrics["mass_balance_relerr_numerical"] < 1e-10

    assert metrics["max_sym_rel_err"] < tols["ts"]

    t50_a = metrics["t50_analytical"]
    t50_n = metrics["t50_numerical"]
    if np.isfinite(t50_a) and np.isfinite(t50_n):
        assert abs(t50_a - t50_n) / max(t50_a, 1e-30) < tols["t50"]

    tp_a = metrics["t_peak_analytical"]
    tp_n = metrics["t_peak_numerical"]
    if np.isfinite(tp_a) and np.isfinite(tp_n):
        assert abs(tp_a - tp_n) / max(tp_a, 1e-30) < tols["tpeak"]


def test_release_regime_labels_are_assigned_as_expected():
    assert _classify_bi_regime(0.1) == "low Bi: external transfer limited"
    assert _classify_bi_regime(10.0) == "intermediate Bi"
    assert _classify_bi_regime(1000.0) == "high Bi: internal diffusion limited"


def test_numerical_converges_with_dt_refinement():
    cfg0 = copy.deepcopy(minimal_config)
    cfg0["dt"] = 200
    cfg0["n_timesteps"] = 50
    cfg0["solver_t_eval"] = "timesteps"
    cfg0["additive_release"] = {
        "model": "numerical",
        "solver": {
            "n_r": 200,
            "n_substeps": 80,
            "theta": 1.0,
        },
    }

    data = copy.deepcopy(minimal_data_with_additive)
    data["k_frag"] = 0.02
    data["initial_additive_concs"] = [1.0] * cfg0["n_size_classes"]
    data["additive_release"] = {
        "D_p": 1e-16,
        "D_w": 1e-9,
        "K_pw": 1e4,
    }

    dts = [200, 100, 50]
    finals = []
    for dt in dts:
        cfg = _make_config_with_dt(cfg0, dt=dt)
        out = FragmentMNP(cfg, data).run()
        finals.append(float(_released_fraction(out)[-1]))

    assert np.all(np.isfinite(finals)), f"Non-finite final released fractions: {finals}"

    err1 = abs(finals[0] - finals[1])
    err2 = abs(finals[1] - finals[2])
    noise_floor = 1e-10
    assert (err2 + noise_floor) < (err1 + noise_floor)
