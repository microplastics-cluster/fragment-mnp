"""
Example config and data (:mod:`fragmentmnp.examples`)
=====================================================

Provides example config and data dictionaries for use
in the FRAGMENT-MNP model.
"""
import numpy as np

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
"""Example model config with all available variables."""

minimal_config = {
    'n_size_classes': 7,
    'particle_size_range': [-9, -3],
    'n_timesteps': 100,
}
"""Example model config with only required variables.
Other variables will take their default values."""


# Example rate constant distribution parameters
def _k_dist_params(dims):
    k_dist_params = {}
    for x in dims:
        k_dist_params[f'A_{x}'] = 1.0
        k_dist_params[f'alpha_{x}'] = 0.0
        k_dist_params[f'B_{x}'] = 1.0
        k_dist_params[f'beta_{x}'] = 0.0
        k_dist_params[f'C_{x}'] = None
        k_dist_params[f'gamma_{x}'] = 1.0
        k_dist_params[f'D_{x}'] = None
        k_dist_params[f'delta1_{x}'] = 1.0
        k_dist_params[f'delta2_{x}'] = None
    return k_dist_params


full_data = {
    'initial_concs': [42.0] * 7,
    'initial_concs_diss': 0.0,
    'density': 1380,              # PET density [kg/m3]
    'k_frag': {'k_f': 0.01,
               'k_0': 0.0,
               'is_compound': True,
               **_k_dist_params(['t', 's'])},
    # No dissolution or mineralisation in example data
    'k_diss': {'k_f': 0.0,
               'k_0': 0.0,
               'is_compound': True,
               **_k_dist_params(['t', 's'])},
    'k_min': {'k_f': 0.0,
              'k_0': 0.0,
              'is_compound': True,
              **_k_dist_params(['t'])},
    'fsd_beta': 0.0
}
"""Example model data with all available variables."""

minimal_data = {
    'initial_concs': [42.0] * 7,
    'density': 1380,              # PET density [kg/m3]
    'k_frag': 0.01,
    'k_min': 0.0  # Added k_min to minimal data
}
"""Example model data with only required variables.
Other variables will take their default values."""
