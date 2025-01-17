"""
Unit tests for model output
"""
import numpy as np

from _mock_output import mock_output, t, c, n, c_diss, c_min


def test_saving_data():
    """
    Test that saving some arbitrary data to the output
    object returns the same arbitrary data
    """
    assert (
        np.array_equal(t, mock_output.t) and
        np.array_equal(c, mock_output.c) and
        np.array_equal(n, mock_output.n) and
        np.array_equal(c_diss, mock_output.c_diss) and
        np.array_equal(c_min, mock_output.c_min)
    )


def test_n_timesteps_size_classes():
    """
    Test that the number of timesteps and size classes
    is set correctly
    """
    assert (
        mock_output.n_timesteps == 100 and
        mock_output.n_size_classes == 2
    )
