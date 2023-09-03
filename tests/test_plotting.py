"""
Unit tests for model output plotting
"""
import numpy as np
from _mock_output import mock_output


def test_specifying_size_classes():
    """
    Test that providing a list of size classes via
    the `size_classes_to_plot` option works as
    expected
    """
    fig, ax = mock_output.plot(size_classes_to_plot=[0, 1],
                               units={'mass': 'mg', 'length': 'µm', 'time': 'days', 'volume': 'l'})
    n_lines = len(ax.get_lines())
    assert n_lines == 2
