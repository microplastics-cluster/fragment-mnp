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
import numpy as np

# Example model config specifying all options
full_config = {
    'n_size_classes': 7,
    'particle_size_range': [-9, -3],
    'n_timesteps': 100,
    'dt': 1,
    'solver_method': 'RK45',
    'solver_rtol': 1e-3,
    'solver_atol': 1e-6,
    'solver_max_step': np.inf,
    'solver_t_eval': 'timesteps'
}

# Example model config specifying only required options
minimal_config = {
    'n_size_classes': 7,
    'particle_size_range': [-9, -3],
    'n_timesteps': 100,
}

# Example rate constant distribution parameters
_k_dist_params = {}
for x in ['t', 's']:
    _k_dist_params[f'A_{x}'] = 1.0
    _k_dist_params[f'alpha_{x}'] = 0.0
    _k_dist_params[f'B_{x}'] = 1.0
    _k_dist_params[f'beta_{x}'] = 0.0
    _k_dist_params[f'C_{x}'] = None
    _k_dist_params[f'gamma_{x}'] = 1.0
    _k_dist_params[f'D_{x}'] = None
    _k_dist_params[f'delta1_{x}'] = 1.0
    _k_dist_params[f'delta2_{x}'] = None

# Example model data
full_data = {
    'initial_concs': [42.0] * 7,
    'initial_concs_diss': 0.0,
    'density': 1380,              # PET density [kg/m3]
    'k_frag': {'k_f': 0.01,
               'k_0': 0.0,
               'is_compound': True,
               **_k_dist_params},
    'k_diss': {'k_f': 0.0,
               'k_0': 0.0,
               'is_compound': True,
               **_k_dist_params},
    'fsd_beta': 0.0
}

# Example model data specifying only required data
minimal_data = {
    'initial_concs': [42.0] * 7,
    'density': 1380,              # PET density [kg/m3]
    'k_frag': 0.01,
}
