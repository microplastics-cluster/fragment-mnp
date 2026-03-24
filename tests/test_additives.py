import numpy as np

from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_config_with_additive, minimal_data_with_additive, minimal_data_with_multi_additives, minimal_data_with_phase2


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

def test_phase2_medium_pools_exist_and_transform():
    out = FragmentMNP(minimal_config, minimal_data_with_phase2).run()
    assert out.c_medium_pool_species is not None
    assert 'Additive A:dissolved_parent' in out.medium_pool_names
    assert 'Additive A:transformed_product' in out.medium_pool_names
    i_parent = out.get_medium_pool_index('Additive A:dissolved_parent')
    i_prod = out.get_medium_pool_index('Additive A:transformed_product')
    assert out.c_medium_pool_species[i_parent, -1] >= 0.0
    assert out.c_medium_pool_species[i_prod, -1] > 0.0

def test_phase2_medium_totals_match_named_medium_pool_sum():
    out = FragmentMNP(minimal_config, minimal_data_with_phase2).run()
    summed = out.c_medium_pool_species.sum(axis=0)
    assert np.allclose(summed, out.c_chem_medium_total[0])
