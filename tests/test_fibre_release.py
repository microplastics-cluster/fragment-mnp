import copy
import numpy as np
from scipy.special import jn_zeros

from fragmentmnp import FragmentMNP
from fragmentmnp.examples import (
    minimal_config_with_fibre,
    minimal_data_with_fibre_additive_analytical,
    minimal_data_with_fibre_additive_numerical,
)


def _cylinder_dirichlet_reference(Fo: float, n_terms: int = 1000) -> float:
    alpha = jn_zeros(0, n_terms)
    return float(np.sum((4.0 / alpha**2) * np.exp(-(alpha**2) * Fo)))


def test_cylindrical_analytical_high_bi_matches_dirichlet_bessel_series():
    R = 10e-6
    D_p = 1e-16
    t = 1e6
    Fo = D_p * t / R**2
    ref = _cylinder_dirichlet_reference(Fo, n_terms=1000)

    got = FragmentMNP._analytical_cylindrical_remaining_fraction(
        t=t,
        r=R,
        D_p=D_p,
        D_w=1e-9,
        K_pw=1e10,
        n_terms=120,
    )
    assert np.isclose(got, ref, rtol=2e-5, atol=2e-6)


def test_cylindrical_numerical_km_zero_has_no_release():
    F = FragmentMNP._numerical_cylindrical_remaining_fraction(
        t=3600.0,
        r=10e-6,
        D_p=1e-16,
        K_pw=1e4,
        k_m=0.0,
        n_r=80,
        n_substeps=20,
        theta=1.0,
    )
    assert np.isclose(F, 1.0, rtol=0.0, atol=1e-12)


def test_cylindrical_numerical_high_transfer_matches_dirichlet_reference():
    R = 10e-6
    D_p = 1e-16
    t = 1e5
    Fo = D_p * t / R**2
    ref = _cylinder_dirichlet_reference(Fo, n_terms=1000)

    got = FragmentMNP._numerical_cylindrical_remaining_fraction(
        t=t,
        r=R,
        D_p=D_p,
        K_pw=1.0,
        k_m=1e6,
        n_r=180,
        n_substeps=80,
        theta=1.0,
    )
    assert np.isclose(got, ref, rtol=2e-2, atol=2e-3)


def test_fibre_analytical_additive_mass_is_conserved_end_to_end():
    out = FragmentMNP(
        minimal_config_with_fibre,
        minimal_data_with_fibre_additive_analytical,
    ).run()
    total0 = out.c_chem_part_total[0, :, 0].sum() + out.c_chem_medium_total[0, 0]
    totalT = out.c_chem_part_total[0, :, -1].sum() + out.c_chem_medium_total[0, -1]
    assert np.isclose(totalT, total0, rtol=1e-10, atol=1e-12)


def test_fibre_numerical_additive_mass_is_conserved_end_to_end():
    cfg = copy.deepcopy(minimal_config_with_fibre)
    cfg['n_timesteps'] = 30
    data = copy.deepcopy(minimal_data_with_fibre_additive_numerical)
    data['additives'][0]['pools'][0]['release']['solver'].update({
        'n_r': 50,
        'n_substeps': 10,
        'theta': 1.0,
    })
    out = FragmentMNP(cfg, data).run()
    total0 = out.c_chem_part_total[0, :, 0].sum() + out.c_chem_medium_total[0, 0]
    totalT = out.c_chem_part_total[0, :, -1].sum() + out.c_chem_medium_total[0, -1]
    assert np.isclose(totalT, total0, rtol=1e-10, atol=1e-12)


def test_constant_diameter_fibre_release_radius_is_independent_of_length():
    fmnp = FragmentMNP(
        minimal_config_with_fibre,
        minimal_data_with_fibre_additive_analytical,
    )
    assert np.allclose(fmnp.release_radii, fmnp.release_radii[0])
