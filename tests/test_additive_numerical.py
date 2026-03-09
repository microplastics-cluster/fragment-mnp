import numpy as np
import pytest

from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data_with_additive


def _series_reference(Fo: float, n_terms: int = 600) -> float:
    """
    Classical diffusion-in-sphere (Dirichlet boundary) reference:
        F = sum_{n=1}^inf (6/(n^2*pi^2)) exp(-n^2*pi^2*Fo)

    Used as a high-Bi / strong-sink benchmark.
    """
    n = np.arange(1, n_terms + 1, dtype=float)
    lam2 = (n * np.pi) ** 2
    terms = (6.0 / lam2) * np.exp(-lam2 * Fo)
    return float(np.sum(terms))


def test_numerical_km_zero_no_release():
    """
    If k_m = 0, there is no flux to water, so remaining fraction must be 1.
    """
    F = FragmentMNP._numerical_additive_remaining_fraction(
        t=3600.0,
        r=1e-4,
        D_p=1e-16,
        K_pw=1e4,
        k_m=0.0,
        n_r=60,
        n_substeps=10,
        theta=1.0
    )
    assert np.isclose(F, 1.0, rtol=0.0, atol=1e-12)


def test_numerical_high_km_matches_dirichlet_series():
    """
    With very large k_m (strong external sink), the Robin BC approaches
    the Dirichlet boundary condition C_s ~ 0.

    In that limit, the remaining fraction should match the classical
    diffusion-in-sphere series reasonably well.
    """
    r = 1e-4
    D_p = 1e-16
    t = 1e6  # seconds (~11.6 days)

    Fo = D_p * t / (r * r)
    F_ref = _series_reference(Fo=Fo, n_terms=2000)

    # huge k_m makes boundary almost Dirichlet
    F_num = FragmentMNP._numerical_additive_remaining_fraction(
        t=t,
        r=r,
        D_p=D_p,
        K_pw=1.0,      # keep simple; high k_m dominates anyway
        k_m=1e6,       # very large [m/s] => near-perfect sink
        n_r=120,
        n_substeps=40,
        theta=1.0
    )

    # Numerical diffusion + discretisation means we allow a modest tolerance.
    # Tighten later if you increase n_r/n_substeps.
    assert np.isclose(F_num, F_ref, rtol=2e-2, atol=2e-3)


def test_full_model_additive_mass_conserved_numerical():
    """
    End-to-end test: when additive tracking is enabled, total additive mass
    should be conserved between particulate and aqueous pools.
    """
    data = dict(minimal_data_with_additive)

    # Ensure we are actually testing the numerical branch
    config = dict(minimal_config)
    config["additive_release"] = {
        "model": "numerical",
        "solver": {
            "n_r": 60,
            "n_substeps": 10,
            "theta": 1.0,
        }
    }

    data["additive_release"] = dict(data["additive_release"])
    data["additive_release"].update({
        "D_w": data["additive_release"].get("D_w", 1e-9),
    })

    out = FragmentMNP(config, data).run()

    assert out.A_part is not None
    assert out.A_aq is not None

    total0 = float(np.sum(out.A_part[:, 0]) + out.A_aq[0])
    totalT = float(np.sum(out.A_part[:, -1]) + out.A_aq[-1])

    assert np.isclose(totalT, total0, rtol=1e-12, atol=1e-10)