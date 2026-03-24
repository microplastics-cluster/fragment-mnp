from fragmentmnp.examples import minimal_config, minimal_data_with_phase2
from fragmentmnp.sensitivity import run_sensitivity_grid


def test_sensitivity_runner_returns_dataframe():
    df = run_sensitivity_grid(
        minimal_config,
        minimal_data_with_phase2,
        {
            "data.k_frag": [0.01, 0.02],
            "data.additives.0.pools.0.release.params.D_p": [1e-17, 1e-16],
        },
    )
    assert len(df) == 4
    assert "fraction_released" in df.columns
    assert "t50_release" in df.columns