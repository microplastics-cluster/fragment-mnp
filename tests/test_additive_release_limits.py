import numpy as np
import pytest
from fragmentmnp.fragmentmnp import FragmentMNP


def _series_reference(Fo: float, n_terms: int = 300) -> float:
    """High-accuracy reference for the Bi->∞ sphere series."""
    n = np.arange(1, n_terms + 1, dtype=float)
    lam2 = (n * np.pi) ** 2
    return float(np.sum((6.0 / lam2) * np.exp(-lam2 * Fo)))


def test_bi_to_zero_limit_release_near_zero():
    """
    Bi -> 0 should imply essentially no release in a timestep:
      F ≈ exp(-Bi*Fo) ≈ 1  => release ≈ 0
    """
    r = 1e-4       # 100 µm radius
    t = 3600.0     # 1 hour

    # Make Bi extremely small by making K_pw*D_w << D_p
    D_p = 1e-12
    D_w = 1e-15
    K_pw = 1e-6
    Bi = (K_pw * D_w) / D_p
    assert Bi < 1e-6

    F = FragmentMNP._analytical_additive_remaining_fraction(
        t=t, r=r, D_p=D_p, D_w=D_w, K_pw=K_pw, n_terms=50
    )

    # Remaining fraction should be extremely close to 1
    assert np.isclose(F, 1.0, rtol=0.0, atol=1e-8)


def test_bi_to_infinity_limit_matches_series():
    """
    Bi -> ∞ should select the classical diffusion-in-sphere series.

    IMPORTANT:
    Use a Fourier number (Fo) where the truncated series converges quickly.
    Very small Fo requires extremely many terms and makes the test unreliable.
    """
    # Use a smaller radius to increase Fo to a convergent regime
    r = 1e-6       # 1 µm radius
    t = 3600.0     # 1 hour

    # Make Bi huge while keeping D_p large enough to avoid tiny Fo
    D_p = 1e-16
    D_w = 1e-9
    K_pw = 1e10
    Bi = (K_pw * D_w) / D_p
    assert Bi > 1e6

    Fo = D_p * t / (r * r)
    assert Fo > 1e-3  # ensure we're not in the slow-convergence regime

    F_impl = FragmentMNP._analytical_additive_remaining_fraction(
        t=t, r=r, D_p=D_p, D_w=D_w, K_pw=K_pw, n_terms=50
    )
    F_ref = _series_reference(Fo=Fo, n_terms=300)

    assert np.isclose(F_impl, F_ref, rtol=1e-6, atol=1e-10)
