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

def _make_inheritance_case(mode, **kwargs):
    data = {
        'initial_concs': [42.0] * 7,
        'density': 1380,
        'k_frag': 0.05,
        'k_min': 0.0,
        'additives': [
            {
                'name': 'Additive A',
                'pools': [
                    {
                        'name': 'Pool 1',
                        'initial_concs': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                        'release': {
                            'model': 'analytical',
                            'params': {
                                'D_p': 1e-12,
                                'D_w': 1e-15,
                                'K_pw': 1e-6,
                            }
                        },
                        'inheritance': {
                            'mode': mode,
                            **kwargs
                        }
                    }
                ],
                'medium_pools': [
                    {'name': 'medium', 'initial_mass': 0.0, 'fate': {}}
                ]
            }
        ]
    }
    return data


def test_inheritance_switch_preserves_total_additive_mass():
    for inh in [
        {"mode": "proportional"},
        {"mode": "size_biased", "beta": -1.0},
        {"mode": "surface_enriched", "gamma": 1.0},
    ]:
        data = {
            'initial_concs': [42.0] * 7,
            'density': 1380,
            'k_frag': 0.05,
            'k_min': 0.0,
            'additives': [
                {
                    'name': 'Additive A',
                    'pools': [
                        {
                            'name': 'Pool 1',
                            'initial_concs': [1.0] * 7,
                            'release': {
                                'model': 'analytical',
                                'params': {'D_p': 1e-16, 'D_w': 1e-9, 'K_pw': 1e4}
                            },
                            'inheritance': inh
                        }
                    ]
                }
            ]
        }

        out = FragmentMNP(minimal_config, data).run()
        total0 = out.c_chem_part_total[0, :, 0].sum() + out.c_chem_medium_total[0, 0]
        totalT = out.c_chem_part_total[0, :, -1].sum() + out.c_chem_medium_total[0, -1]
        assert np.isclose(total0, totalT, rtol=1e-10, atol=1e-12)


def test_proportional_and_size_biased_modes_give_different_size_distributions():
    out_prop = FragmentMNP(
        minimal_config,
        _make_inheritance_case("proportional")
    ).run()

    out_bias = FragmentMNP(
        minimal_config,
        _make_inheritance_case("size_biased", beta=-2.0)
    ).run()

    part_prop = out_prop.c_chem_part_total[0, :, -1]
    part_bias = out_bias.c_chem_part_total[0, :, -1]

    total0_prop = out_prop.c_chem_part_total[0, :, 0].sum() + out_prop.c_chem_medium_total[0, 0]
    totalT_prop = out_prop.c_chem_part_total[0, :, -1].sum() + out_prop.c_chem_medium_total[0, -1]

    total0_bias = out_bias.c_chem_part_total[0, :, 0].sum() + out_bias.c_chem_medium_total[0, 0]
    totalT_bias = out_bias.c_chem_part_total[0, :, -1].sum() + out_bias.c_chem_medium_total[0, -1]

    assert np.isclose(totalT_prop, total0_prop, rtol=1e-10, atol=1e-12)
    assert np.isclose(totalT_bias, total0_bias, rtol=1e-10, atol=1e-12)

    assert not np.allclose(part_prop, part_bias)


def test_size_biased_negative_beta_enriches_smaller_sizes_relative_to_proportional():
    out_prop = FragmentMNP(
        minimal_config,
        _make_inheritance_case("proportional")
    ).run()

    out_bias = FragmentMNP(
        minimal_config,
        _make_inheritance_case("size_biased", beta=-2.0)
    ).run()

    part_prop = out_prop.c_chem_part_total[0, :, -1]
    part_bias = out_bias.c_chem_part_total[0, :, -1]

    # smallest size class should hold relatively more additive under negative beta
    assert part_bias[0] > part_prop[0]

def test_size_dependent_kdeg_removes_more_from_smaller_classes():
    data = {
        'initial_concs': [42.0] * 7,
        'density': 1380,
        'k_frag': 0.0,
        'k_min': 0.0,
        'additives': [
            {
                'name': 'Additive A',
                'pools': [
                    {
                        'name': 'Pool 1',
                        'initial_concs': [1.0] * 7,
                        'release': {
                            'model': 'analytical',
                            'params': {
                                'D_p': 1e-12,
                                'D_w': 1e-15,
                                'K_pw': 1e-6,
                            }
                        },
                        'fate': {
                            'k_deg': [0.1, 0.08, 0.06, 0.04, 0.02, 0.01, 0.0]
                        }
                    }
                ]
            }
        ]
    }

    out = FragmentMNP(minimal_config, data).run()
    part_final = out.c_chem_part_total[0, :, -1]

    # Smaller classes should be depleted more strongly
    assert part_final[0] < part_final[-1]

def test_size_dependent_transfer_moves_more_from_smaller_classes():
    data = {
        'initial_concs': [42.0] * 7,
        'density': 1380,
        'k_frag': 0.0,
        'k_min': 0.0,
        'additives': [
            {
                'name': 'Additive A',
                'pools': [
                    {
                        'name': 'Pool 1',
                        'initial_concs': [1.0] * 7,
                        'release': {
                            'model': 'analytical',
                            'params': {
                                'D_p': 1e-12,
                                'D_w': 1e-15,
                                'K_pw': 1e-6,
                            }
                        },
                        'fate': {
                            'transfers': [
                                {
                                    'to': 'Pool 2',
                                    'k': [0.1, 0.08, 0.06, 0.04, 0.02, 0.01, 0.0]
                                }
                            ]
                        }
                    },
                    {
                        'name': 'Pool 2',
                        'initial_concs': [0.0] * 7,
                        'release': {
                            'model': 'analytical',
                            'params': {
                                'D_p': 1e-12,
                                'D_w': 1e-15,
                                'K_pw': 1e-6,
                            }
                        }
                    }
                ]
            }
        ]
    }

    out = FragmentMNP(minimal_config, data).run()

    i1 = out.get_species_index('Additive A:Pool 1')
    i2 = out.get_species_index('Additive A:Pool 2')

    pool2_final = out.c_chem_part_species[i2, :, -1]

    # Smaller classes should have received more transferred mass
    assert pool2_final[0] > pool2_final[-1]