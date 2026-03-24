import numpy as np

from fragmentmnp import FragmentMNP
from fragmentmnp.examples import (
    minimal_config,
    minimal_config_with_additive,
    minimal_data_with_additive,
    minimal_data_with_multi_additives,
    minimal_data_with_fate,
)


def test_additive_outputs_exist_when_enabled():
    out = FragmentMNP(minimal_config_with_additive, minimal_data_with_additive).run()
    assert out.A_part is not None
    assert out.A_aq is not None


def test_additive_mass_conserved_total_particulate_plus_aqueous():
    """
    With analytical release enabled, particulate additive decreases but
    aqueous additive increases. Total additive mass must be conserved.
    """
    out = FragmentMNP(minimal_config_with_additive, minimal_data_with_additive).run()

    total0 = np.sum(minimal_data_with_additive["initial_additive_concs"])
    totalT = out.A_part[:, -1].sum() + out.A_aq[-1]

    assert np.isclose(total0, totalT, rtol=1e-10, atol=1e-12)

def test_multi_additive_outputs_exist():
    out = FragmentMNP(minimal_config, minimal_data_with_multi_additives).run()

    assert out.c_chem_part_species is not None
    assert out.c_chem_medium_species is not None
    assert out.c_chem_part_total is not None
    assert out.c_chem_medium_total is not None


def test_multi_additive_shapes():
    out = FragmentMNP(minimal_config, minimal_data_with_multi_additives).run()

    # AO168 has 2 pools, UV328 has 1 pool
    assert out.c_chem_part_species.shape[0] == 3
    assert out.c_chem_medium_species.shape[0] == 3
    assert out.c_chem_part_total.shape[0] == 2
    assert out.c_chem_medium_total.shape[0] == 2


def test_multi_additive_mass_conserved_per_additive():
    out = FragmentMNP(minimal_config, minimal_data_with_multi_additives).run()

    for a_idx in range(out.c_chem_part_total.shape[0]):
        total0 = out.c_chem_part_total[a_idx, :, 0].sum() + out.c_chem_medium_total[a_idx, 0]
        totalT = out.c_chem_part_total[a_idx, :, -1].sum() + out.c_chem_medium_total[a_idx, -1]
        assert np.isclose(total0, totalT, rtol=1e-10, atol=1e-12)


def test_species_sum_matches_additive_total():
    out = FragmentMNP(minimal_config, minimal_data_with_multi_additives).run()

    # species 0 and 1 belong to AO168
    np.testing.assert_allclose(
        out.c_chem_part_species[0] + out.c_chem_part_species[1],
        out.c_chem_part_total[0]
    )
    np.testing.assert_allclose(
        out.c_chem_medium_species[0] + out.c_chem_medium_species[1],
        out.c_chem_medium_total[0]
    )

def test_zero_fate_rates_preserve_total_mass_conservation():
    out = FragmentMNP(minimal_config, minimal_data_with_multi_additives).run()

    total0 = out.c_chem_part_total[0, :, 0].sum() + out.c_chem_medium_total[0, 0]
    totalT = out.c_chem_part_total[0, :, -1].sum() + out.c_chem_medium_total[0, -1]
    assert np.isclose(total0, totalT, rtol=1e-10, atol=1e-12)


def test_pool_transfer_moves_mass_to_target_pool():
    out = FragmentMNP(minimal_config, minimal_data_with_fate).run()

    pool1_initial = out.c_chem_part_species[0, :, 0].sum()
    pool2_initial = out.c_chem_part_species[1, :, 0].sum()
    pool2_later = out.c_chem_part_species[1, :, -1].sum()

    assert np.isclose(pool2_initial, 0.0)
    assert pool1_initial > 0.0
    assert pool2_later > 0.0


def test_k_deg_and_k_loss_reduce_modeled_total_mass():
    out = FragmentMNP(minimal_config, minimal_data_with_fate).run()

    total0 = out.c_chem_part_total[0, :, 0].sum() + out.c_chem_medium_total[0, 0]
    totalT = out.c_chem_part_total[0, :, -1].sum() + out.c_chem_medium_total[0, -1]

    assert totalT < total0
