import copy
import numpy as np

from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data_with_components
from fragmentmnp.validation import validate_config, validate_data


def test_component_input_validates_and_defaults():
    config = copy.deepcopy(minimal_config)
    data = copy.deepcopy(minimal_data_with_components)
    # Exercise component-only input by removing aggregate fields that can be derived.
    data.pop('initial_concs')
    data.pop('density')
    data.pop('k_frag')

    validated = validate_data(data, validate_config(config))
    assert validated['components'] is not None
    assert len(validated['components']) == 3
    summed = np.sum([c['initial_concs'] for c in validated['components']], axis=0)
    np.testing.assert_allclose(validated['initial_concs'], summed)


def test_component_outputs_are_available_and_aggregate_correctly():
    out = FragmentMNP(minimal_config, minimal_data_with_components).run()

    assert out.c_component is not None
    assert out.c_component.shape[0] == 3
    assert out.c_component.shape[1] == minimal_config['n_size_classes']
    np.testing.assert_allclose(out.c, out.c_component.sum(axis=0))
    np.testing.assert_allclose(out.c_diss, out.c_diss_component.sum(axis=0))
    np.testing.assert_allclose(out.c_min, out.c_min_component.sum(axis=0))


def test_component_mass_is_conserved_without_dissolution():
    out = FragmentMNP(minimal_config, minimal_data_with_components).run()
    for j in range(out.c_component.shape[0]):
        total = out.c_component[j].sum(axis=0) + out.c_diss_component[j] + out.c_min_component[j]
        assert np.allclose(total, total[0], rtol=5e-2, atol=1e-10)


def test_layer_thickness_is_tracked():
    out = FragmentMNP(minimal_config, minimal_data_with_components).run()
    assert out.layer_thickness_component is not None
    assert out.layer_thickness_component.shape[0] == 3
    ts = out.get_component_timeseries('PE_outer_layer')
    assert 'layer_thickness' in ts
    assert ts['layer_thickness'][0] > 0.0


def test_component_specific_fragmentation_rates_change_outputs():
    base = copy.deepcopy(minimal_data_with_components)
    out_component = FragmentMNP(minimal_config, base).run()

    same_rate = copy.deepcopy(minimal_data_with_components)
    for comp in same_rate['components']:
        comp['k_frag'] = 0.010
    out_same = FragmentMNP(minimal_config, same_rate).run()

    assert not np.allclose(
        out_component.get_component_timeseries('EVOH_barrier_layer')['particulate_by_size'],
        out_same.get_component_timeseries('EVOH_barrier_layer')['particulate_by_size'],
    )
