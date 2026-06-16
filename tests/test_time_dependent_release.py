import copy
import numpy as np

from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config
from fragmentmnp.validation import validate_config, validate_data


def _base_time_dependent_data(D_p):
    return {
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
                            'solver': {'n_terms': 200},
                            'params': {
                                'D_p': D_p,
                                'D_w': 1e-9,
                                'K_pw': 1e10,
                            },
                        },
                    }
                ],
                'medium_pools': [
                    {'name': 'medium', 'initial_mass': 0.0, 'fate': {}},
                ],
            }
        ],
    }


def test_time_dependent_release_parameter_validates():
    config = copy.deepcopy(minimal_config)
    data = _base_time_dependent_data(
        {
            'times': [0.0, 25.0, 50.0, 100.0],
            'values': [1e-20, 1e-18, 1e-16, 1e-15],
        }
    )

    validated = validate_data(data, validate_config(config))
    D_p = validated['additives'][0]['pools'][0]['release']['params']['D_p']
    assert isinstance(D_p, dict)
    assert D_p['times'][-1] == 100.0


def test_time_dependent_Dp_increases_final_release():
    config = copy.deepcopy(minimal_config)
    config['n_timesteps'] = 100
    config['dt'] = 1

    # Always-low diffusivity case.
    data_low = _base_time_dependent_data(1e-20)
    out_low = FragmentMNP(config, data_low).run()

    # Weathering case: D_p increases with model time.
    data_td = _base_time_dependent_data(
        {
            'times': [0.0, 25.0, 50.0, 100.0],
            'values': [1e-20, 1e-18, 1e-16, 1e-15],
        }
    )
    out_td = FragmentMNP(config, data_td).run()

    assert out_td.c_chem_medium_total[0, -1] > out_low.c_chem_medium_total[0, -1]


def test_time_dependent_size_resolved_release_parameter_runs_and_conserves_mass():
    config = copy.deepcopy(minimal_config)
    config['n_timesteps'] = 60
    config['dt'] = 1
    n = config['n_size_classes']

    values = np.vstack([
        np.full(n, 1e-20),
        np.logspace(-18, -15, n),
        np.logspace(-17, -14, n),
    ])

    data = _base_time_dependent_data(
        {
            'times': [0.0, 30.0, 60.0],
            'values': values.tolist(),  # shape = (time, size)
        }
    )

    out = FragmentMNP(config, data).run()
    total0 = out.c_chem_part_total[0, :, 0].sum() + out.c_chem_medium_total[0, 0]
    totalT = out.c_chem_part_total[0, :, -1].sum() + out.c_chem_medium_total[0, -1]
    assert np.isclose(totalT, total0, rtol=1e-10, atol=1e-12)


def test_parametric_linear_release_parameter_runs():
    config = copy.deepcopy(minimal_config)
    data = _base_time_dependent_data(
        {
            'type': 'linear',
            'initial': 1e-20,
            'rate': 1e-1,
            'minimum': 1e-20,
            'maximum': 1e-15,
        }
    )

    out = FragmentMNP(config, data).run()
    assert out.c_chem_medium_total[0, -1] >= 0.0
