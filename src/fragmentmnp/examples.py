"""
Example config and data (:mod:`fragmentmnp.examples`)
=====================================================

This module provides example config and data dictionaries for use
with FRAGMENT-MNP.

We keep the original examples intact, and add a *new* example that
activates additive tracking and analytical release coupling.
"""
import numpy as np


# -------------------------
# Example model configuration
# -------------------------

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

minimal_config_with_additive = {
    **minimal_config,
    'additive_release': {
        'model': 'analytical',
        'solver': {
            'n_terms': 50
        }
    }
}
"""Example model config with only required variables.
Other variables will take their default values."""


# -------------------------
# Example rate constant distribution parameters
# -------------------------

def _k_dist_params(dims):
    """
    Helper function to create default regression parameter dict entries
    for the FRAGMENT-MNP rate-constant distribution builder.
    """
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
    'k_frag': 0.01
}
"""Example model data with only required variables.
Other variables will take their default values."""


# --------------------------------------------------------------------
# Legacy single-additive example
# --------------------------------------------------------------------
minimal_data_with_additive = {
    'initial_concs': [42.0] * 7,
    'density': 1380,
    'k_frag': 0.01,
    'k_min': 0.0,

    # Additive mass concentration in each size class at t=0
    # (same binning as initial_concs)
    'initial_additive_concs': [1.0] * 7,

    # Additive release model configuration.
    # This is where analytical and/or numerical solution parameters live.
    #
    # IMPORTANT:
    # - The actual analytical and/or numerical formula is implemented in
    #   FragmentMNP
    # - Here the parameter values are stored

    'additive_release': {
        'model': 'analytical',
        'params': {
            'D_p': 1e-16,    # diffusion in polymer (m2/s)
            'D_w': 1e-9,     # diffusion in water (m2/s)
            'K_pw': 1e4,     # polymer-water partition coefficient (-)
        }
    }
}
"""Example data that activates single additive bookkeeping + analytical release."""

# --------------------------------------------------------------------
# New canonical multi-additive / multi-pool example
# --------------------------------------------------------------------
minimal_data_with_multi_additives = {
    'initial_concs': [42.0] * 7,
    'density': 1380,
    'k_frag': 0.01,
    'k_min': 0.0,

    'additives': [
        {
            'name': 'Additive A',

            'pools': [
                {
                    'name': 'Pool 1',

                    'initial_concs': [0.3] * 7,

                    'release': {
                        'model': 'analytical',
                        'solver': {'n_terms': 50},
                        'params': {
                            'D_p': 1e-16,
                            'D_w': 1e-9,
                            'K_pw': 1e4,
                        }
                    },

                    'inheritance': {
                        'mode': 'proportional'
                    }
                },

                {
                    'name': 'Pool 2',

                    'initial_concs': [0.7] * 7,

                    'release': {
                        'model': 'numerical',

                        'solver': {
                            'n_r': 60,
                            'n_substeps': 20,
                            'theta': 1.0,
                        },

                        'params': {
                            'D_p': 1e-18,
                            'K_pw': 1e5,
                            'k_m': 1e-8,
                        }
                    },

                    'inheritance': {
                        'mode': 'size_biased',
                        'beta': -1.0
                    }
                }
            ]
        },

        {
            'name': 'Additive B',

            'pools': [
                {
                    'name': 'Pool 1',

                    'initial_concs': [0.5] * 7,

                    'release': {
                        'model': 'analytical',
                        'solver': {'n_terms': 50},
                        'params': {
                            'D_p': 1e-15,
                            'D_w': 1e-9,
                            'K_pw': 1e3,
                        }
                    },

                    'inheritance': {
                        'mode': 'surface_enriched',
                        'gamma': 1.0
                    }
                }
            ]
        }
    ]
}
"""Canonical multi-additive / multi-pool example."""


# --------------------------------------------------------------------
# Phase 2 example with named medium pools and transformed-product generation
# --------------------------------------------------------------------
minimal_data_with_phase2 = {
    'initial_concs': [42.0] * 7,
    'density': 1380,
    'k_frag': 0.01,
    'k_min': 0.0,
    'additives': [
        {
            'name': 'Additive A',
            'pools': [
                {
                    'name': 'fast_domain',
                    'initial_concs': [0.8] * 7,
                    'release': {
                        'model': 'analytical',
                        'target': 'Additive A:dissolved_parent',
                        'solver': {'n_terms': 50},
                        'params': {'D_p': 1e-16, 'D_w': 1e-9, 'K_pw': 1e4},
                    },
                    'inheritance': {
                        'mode': 'proportional'
                    },
                    'fate': {
                        'k_deg': [1e-5, 2e-5, 4e-5, 8e-5, 1.6e-4, 3.2e-4, 6.4e-4],
                        'transfers': [
                            {'to': 'slow_domain', 'k': [1e-4, 1e-4, 2e-4, 3e-4, 5e-4, 8e-4, 1e-3]}
                        ]
                    },
                },
                {
                    'name': 'slow_domain',
                    'initial_concs': [0.2] * 7,
                    'release': {
                        'model': 'analytical',
                        'target': 'Additive A:dissolved_parent',
                        'solver': {'n_terms': 50},
                        'params': {'D_p': 1e-18, 'D_w': 1e-9, 'K_pw': 1e5},
                    },
                    'inheritance': {
                        'mode': 'size_biased',
                        'beta': -1.0
                    },
                    'fate': {'transfers': [{'to': 'fast_domain', 'k': 2e-4}]},
                },
            ],
            'medium_pools': [
                {'name': 'dissolved_parent', 'initial_mass': 0.0, 'fate': {'transfers': [{'to': 'transformed_product', 'k': 5e-3}]}},
                {'name': 'transformed_product', 'initial_mass': 0.0, 'fate': {'k_loss': 1e-4}},
                {'name': 'sorbed', 'initial_mass': 0.0, 'fate': {}},
            ],
        }
    ],
}

# Optional fragmentation inheritance rule for additive redistribution
# during polymer fragmentation:
# - proportional
# - size_biased (beta)
# - surface_enriched (gamma)
"""Example data showing alternative fragmentation inheritance modes."""
minimal_data_with_inheritance = {
    'initial_concs': [42.0] * 7,
    'density': 1380,
    'k_frag': 0.01,
    'k_min': 0.0,
    'additives': [
        {
            'name': 'Additive A',
            'pools': [
                {
                    'name': 'reference_pool',
                    'initial_concs': [1.0] * 7,
                    'release': {
                        'model': 'analytical',
                        'solver': {'n_terms': 50},
                        'params': {
                            'D_p': 1e-16,
                            'D_w': 1e-9,
                            'K_pw': 1e4,
                        }
                    },
                    'inheritance': {
                        'mode': 'proportional'
                    }
                },
                {
                    'name': 'small_fragment_enriched_pool',
                    'initial_concs': [1.0] * 7,
                    'release': {
                        'model': 'analytical',
                        'solver': {'n_terms': 50},
                        'params': {
                            'D_p': 1e-16,
                            'D_w': 1e-9,
                            'K_pw': 1e4,
                        }
                    },
                    'inheritance': {
                        'mode': 'size_biased',
                        'beta': -1.0
                    }
                },
                {
                    'name': 'surface_enriched_pool',
                    'initial_concs': [1.0] * 7,
                    'release': {
                        'model': 'analytical',
                        'solver': {'n_terms': 50},
                        'params': {
                            'D_p': 1e-16,
                            'D_w': 1e-9,
                            'K_pw': 1e4,
                        }
                    },
                    'inheritance': {
                        'mode': 'surface_enriched',
                        'gamma': 1.0
                    }
                }
            ]
        }
    ]
}

# --------------------------------------------------------------------
# Example with time-dependent additive release parameters
# --------------------------------------------------------------------
minimal_data_with_time_dependent_release = {
    'initial_concs': [42.0] * 7,
    'density': 1380,
    'k_frag': 0.01,
    'k_min': 0.0,
    'additives': [
        {
            'name': 'Weathering-sensitive additive',
            'pools': [
                {
                    'name': 'matrix_pool',
                    'initial_concs': [1.0] * 7,
                    'release': {
                        'model': 'analytical',
                        'solver': {'n_terms': 100},
                        'params': {
                            # D_p increases with weathering/ageing time.
                            # Times use the same units as the model time grid.
                            'D_p': {
                                'times': [0.0, 25.0, 50.0, 100.0],
                                'values': [1e-20, 1e-18, 1e-16, 1e-15],
                            },
                            'D_w': 1e-9,
                            'K_pw': 1e10,
                        }
                    },
                    'inheritance': {'mode': 'proportional'},
                }
            ],
            'medium_pools': [
                {'name': 'medium', 'initial_mass': 0.0, 'fate': {}},
            ],
        }
    ],
}
"""Example where polymer diffusivity D_p evolves over model time."""

# --------------------------------------------------------------------
# Component / multilayer packaging example
# --------------------------------------------------------------------
minimal_data_with_components = {
    # Legacy aggregate fields are optional when components are provided, but
    # included here to make the aggregate/component consistency explicit.
    'initial_concs': [40.0, 25.0, 10.0, 0.0, 0.0, 0.0, 0.0],
    'density': 1012.5,
    'k_frag': 0.0,
    'k_diss': 0.0,
    'k_min': 0.0,
    'fsd_beta': 0.0,
    'components': [
        {
            'name': 'PE_outer_layer',
            'initial_concs': [25.0, 15.0, 5.0, 0.0, 0.0, 0.0, 0.0],
            'density': 930.0,
            'k_frag': 0.015,
            'k_diss': 0.0,
            'k_min': 0.0,
            'fsd_beta': -0.2,
            'layer_thickness': 40e-6,
        },
        {
            'name': 'EVOH_barrier_layer',
            'initial_concs': [5.0, 3.0, 2.0, 0.0, 0.0, 0.0, 0.0],
            'density': 1190.0,
            'k_frag': 0.004,
            'k_diss': 0.0,
            'k_min': 0.0,
            'fsd_beta': 0.0,
            'layer_thickness': 5e-6,
        },
        {
            'name': 'PP_inner_layer',
            'initial_concs': [10.0, 7.0, 3.0, 0.0, 0.0, 0.0, 0.0],
            'density': 900.0,
            'k_frag': 0.010,
            'k_diss': 0.0,
            'k_min': 0.0,
            'fsd_beta': -0.1,
            'layer_thickness': 25e-6,
        },
    ],
}
"""Example data for multilayer packaging with component-specific density,
fragmentation rates and tracked layer thickness."""
