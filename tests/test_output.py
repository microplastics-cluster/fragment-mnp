"""
Unit tests for model output, including additive summary helpers.
"""
import numpy as np

from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data_with_multi_additives
from _mock_output import mock_output, t, c, n, c_diss, c_min


def test_saving_data():
    assert (
        np.array_equal(t, mock_output.t) and
        np.array_equal(c, mock_output.c) and
        np.array_equal(n, mock_output.n) and
        np.array_equal(c_diss, mock_output.c_diss) and
        np.array_equal(c_min, mock_output.c_min)
    )


def test_n_timesteps_size_classes():
    assert (
        mock_output.n_timesteps == 100 and
        mock_output.n_size_classes == 2
    )


def test_summary_helpers_multi_additive():
    out = FragmentMNP(minimal_config, minimal_data_with_multi_additives).run()

    summary = out.summary_records(level='additive')
    assert len(summary) == len(out.additive_names)
    assert all('fraction_released' in row for row in summary)
    assert all('t50' in row for row in summary)
    assert all('conservation_residual' in row for row in summary)

    size_rows = out.size_class_contribution_records(level='additive')
    assert len(size_rows) == len(out.additive_names) * out.n_size_classes


def test_additive_lookup_helpers():
    out = FragmentMNP(minimal_config, minimal_data_with_multi_additives).run()

    idx = out.get_additive_index('Additive A')
    assert idx == 0

    species_idx = out.get_species_index('Additive A:Pool 1')
    assert species_idx == 0

    ts = out.get_additive_timeseries('Additive A')
    assert 'particulate_by_size' in ts
    assert 'aqueous' in ts
    assert ts['particulate_by_size'].shape[1] == out.n_timesteps
