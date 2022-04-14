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
