"""
Example config and data for use in the FRAGMENT-MNP model
"""

# Example model configuration
config = {
    'n_size_classes': 7,
    'particle_size_range': [-9, -3],
    'n_timesteps': 100,
    'dt': 1
}

data = {
    'initial_concs': [42.0] * 7,
    'k_frag': 0.01
}
