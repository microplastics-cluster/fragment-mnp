"""
Unit tests for model output
"""
import numpy as np
from fragmentmnp.output import FMNPOutput


# Generate some arbitrary data
t = np.arange(0, 100)
size_bins = np.array([1e-6, 1e-3])
c = np.array([
    np.linspace(0,42,100),
    np.linspace(0,10,100)
])
n = c / (4.0 / 3.0) * np.pi * size_bins[:, np.newaxis] ** 3
c_diss = np.linspace(0, 1e6, 100)
n_diss = np.linspace(0, 1e-2, 100)
# Create output based on this
example_output = FMNPOutput(t, c, n, c_diss, n_diss)


def test_saving_data():
    """
    Test that saving some arbitrary data to the output
    object returns the same arbitrary data
    """
    assert (
        np.array_equal(t, example_output.t) and
        np.array_equal(c, example_output.c) and
        np.array_equal(n, example_output.n) and
        np.array_equal(c_diss, example_output.c_diss) and
        np.array_equal(n_diss, example_output.n_diss)
    )


def test_n_timesteps_size_classes():
    """
    Test that the number of timesteps and size classes
    is set correctly 
    """
    assert (
        example_output.n_timesteps == 100 and
        example_output.n_size_classes == 2
    )
