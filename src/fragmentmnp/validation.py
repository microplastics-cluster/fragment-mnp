"""
Validation of config and data (:mod:`fragmentmnp.validation`)
=============================================================

Provides config and input data validation for the FRAGMENT-MNP model
"""
from schema import Schema, Or, And, Optional

# The schema that particle size ranges should follow
particle_size_range_schema = And(Or((int, float), [int, float]),
                                 And(lambda d: len(d) == 2,
                                     error='particle_size_range must ' +
                                           'be a length-2 iterable'))


# The schema that the config dict should follow
config_schema = Schema({
    # There should be <= 100 size classes
    'n_size_classes': And(int, lambda d: d <= 100),
    # Size range should be a length-2 iterable of type int or float
    Optional('particle_size_range'): particle_size_range_schema,
    # Size classes should be a list of ints of floats
    Optional('particle_size_classes'): And([Or(int, float)]),
    # Timesteps should be an integer
    'n_timesteps': int,
    # Length of timesteps should be an integer (unit of seconds)
    Optional('dt', default=1): int
})


# The schema that the data dict should follow
data_schema = Schema({
    # Initial concs must be a list and >= 0
    'initial_concs': [And(Or(int, float), lambda x: x >= 0.0)],
    # Density must either be a float/int and greater than 0
    'density': And(Or(int, float), lambda x: x >= 0.0),
    # k_frag must either be a float/int, or a list of flaots/ints,
    # and greater than 0
    'k_frag': And(Or(int, float), lambda x: x >= 0.0),
    # theta1 (surface energy empirical parameter) must be a float or int
    Optional('theta_1', default=0.0): Or(int, float),
    # k_diss (dissolution) must be int or float
    Optional('k_diss', default=0.0): Or(int, float)
})


def validate_config(config: dict) -> dict:
    """
    Validate the given config dict against required schema

    Parameters
    ----------
    config : dict
        Model config options

    Returns
    -------
    dict
        Validated config dict
    """
    # Run the validation, and return the validated dict
    # if it passes
    validated = Schema(config_schema).validate(config)
    return validated


def validate_data(data: dict) -> dict:
    """
    Validate the given data dict against required schema

    Parameters
    ----------
    data : dict
        Model input data

    Returns
    -------
    dict
        Validated data dict
    """
    # Run the validation, and return the validated dict
    # if it passes
    validated = Schema(data_schema).validate(data)
    return validated
