"""
Integration tests for the full model
"""
import numpy as np
from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data


def test_model_init():
    """
    Test the model is correctly initialised by checking a few
    attributes after creation
    """
    fmnp = FragmentMNP(minimal_config, minimal_data)
    assert (
        fmnp.n_size_classes == minimal_config['n_size_classes'] and
        fmnp.n_timesteps == minimal_config['n_timesteps']
    )


def test_model_run():
    """
    Test the model run by check the outputs
    """
    output = FragmentMNP(minimal_config, minimal_data).run()
    assert (
        np.array_equal(output.t, np.arange(minimal_config['n_timesteps'])) and
        output.n.sum() == 29400.0
    )


def test_fsd():
    """
    Test the fragment size distribution is calculated correctly
    when an even split amongst daughter size classes is assumed
    """
    N = minimal_config['n_size_classes']
    fsd = np.zeros((N, N))
    for i in np.arange(N):
        fsd[i, :] = 1 / i if i != 0 else 0
    fsd = np.tril(fsd, k=-1)
    fmnp = FragmentMNP(minimal_config, minimal_data)
    # Do the assertion
    np.testing.assert_array_equal(fmnp.fsd, fsd)


def test_f_surface_area():
    """
    Test the surface area scaling factor is calculated correctly
    """
    # With default parameters, f_surface_area should equal this
    f_surface_area = np.logspace(3, -3, 7)
    fmnp = FragmentMNP(minimal_config, minimal_data)
    model_f_surface_area = fmnp._f_surface_area(fmnp.psd)
    # Do the assertion
    np.testing.assert_allclose(model_f_surface_area, f_surface_area)


def test_k_diss_scaling_method_equivalence():
    """
    Test that setting gamma=0 when the scaling method is surface_area
    results in the same output as constant scaling method
    """
    config_sa = minimal_config.copy()
    config_sa['k_diss_scaling_method'] = 'surface_area'
    data_sa = minimal_data.copy()
    data_sa['k_diss'] = 0.0001
    data_sa['k_diss_gamma'] = 0
    data_c = minimal_data.copy()
    data_c['k_diss'] = 0.0001
    output_sa = FragmentMNP(config_sa, data_sa).run()
    output_c = FragmentMNP(minimal_config, data_c).run()
    # Check the output is the same
    assert (
        np.array_equal(output_sa.n, output_c.n) and
        np.array_equal(output_sa.c_diss, output_c.c_diss)
    )
