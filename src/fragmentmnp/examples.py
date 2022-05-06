"""
Example config and data (:mod:`fragmentmnp.examples`)
=====================================================

Provides example config and data dictionaries for use
in the FRAGMENT-MNP model.

Dictionaries
------------
- ``full_config``: Config with all variables complete
- ``minimal_config``: Config with only variables that
  don't have defaults
- ``full_data``: Data with all variables complete
- ``mininal_data``: Data with only variables that don't
  have defaults
"""

# Example model config specifying all options
full_config = {
    'n_size_classes': 7,
    'particle_size_range': [-9, -3],
    'n_timesteps': 100,
    'dt': 1,
    'allow_loss': False
}

# Example model config specifying only required options
minimal_config = {
    'n_size_classes': 7,
    'particle_size_range': [-9, -3],
    'n_timesteps': 100,
}

# Example model data
full_data = {
    'initial_concs': [42.0] * 7,
    'k_frag': 0.01,
    'theta_1': 0.0
}

# Example model data specifying only required data
minimal_data = {
    'initial_concs': [42.0] * 7,
    'k_frag': 0.01,
}
