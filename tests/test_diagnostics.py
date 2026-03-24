from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data_with_phase2
from fragmentmnp.diagnostics import plot_additive_fate_dashboard


def test_diagnostic_plot_runs():
    out = FragmentMNP(minimal_config, minimal_data_with_phase2).run()
    fig, ax = plot_additive_fate_dashboard(out, additive=0)
    assert fig is not None
    assert ax is not None