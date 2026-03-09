import numpy as np

from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config_with_additive, minimal_data_with_additive


def test_additive_outputs_exist_when_enabled():
    out = FragmentMNP(minimal_config_with_additive, minimal_data_with_additive).run()
    assert out.A_part is not None
    assert out.A_aq is not None


def test_additive_mass_conserved_total_particulate_plus_aqueous():
    """
    With analytical release enabled, particulate additive decreases but
    aqueous additive increases. Total additive mass must be conserved.
    """
    out = FragmentMNP(minimal_config_with_additive, minimal_data_with_additive).run()

    total0 = np.sum(minimal_data_with_additive["initial_additive_concs"])
    totalT = out.A_part[:, -1].sum() + out.A_aq[-1]

    assert np.isclose(total0, totalT, rtol=1e-10, atol=1e-12)