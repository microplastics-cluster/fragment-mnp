"""
Validation of config and data (:mod:`fragmentmnp.validation`)
=============================================================

Provides config and input data validation for the FRAGMENT-MNP model
"""
from schema import Schema, Or, And, Optional
from ._errors import FMNPIncorrectDistributionLength

# The schema that particle size ranges should follow
particle_size_range_schema = And(Or((int, float), [int, float]),
                                 And(lambda d: len(d) == 2,
                                     error='particle_size_range must ' +
                                           'be a length-2 iterable'))


def _is_positive_array(arr):
    """
    Check if arr is iterable and all elements
    are positive
    """
    is_array = True
    try:
        for x in arr:
            if x < 0.0:
                is_array = False
    except TypeError:
        is_array = False
    return is_array


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
    Optional('dt', default=1): int,
    # How should the average k_diss provided be scaled across
    # size classes?
    Optional('k_diss_scaling_method',
             default='constant'): lambda x: x in ('constant', 'surface_area')
})


# The schema that the data dict should follow
data_schema = Schema({
    # Initial concs must be a list and >= 0
    'initial_concs': _is_positive_array,
    # Density must either be a float/int and greater than 0
    'density': And(Or(int, float), lambda x: x >= 0.0),
    # k_frag must either be a float/int, or a list of floats/ints,
    # and greater than 0
    'k_frag': Or(And(Or(int, float), lambda x: x >= 0.0), _is_positive_array),
    # theta1 (surface energy empirical parameter) must be a float or int
    Optional('theta_1', default=0.0): Or(int, float),
    # k_diss (dissolution) must be int or float
    Optional('k_diss', default=0.0): Or(And(Or(int, float),
                                            lambda x: x >= 0.0),
                                            _is_positive_array),
    # k_diss_gamma is an empirical param that linearly scales
    # the affect of surface area on dissolution rates
    Optional('k_diss_gamma', default=1.0): Or(int, float)
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


def validate_data(data: dict, config: dict) -> dict:
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

    # If k_frag isn't a scalar (i.e. we're not calculating
    # a k_frag distribution interally), then check the
    # provided distribution is the correct length
    if not isinstance(validated['k_frag'], (int, float)):
        if len(data['k_frag']) != config['n_size_classes']:
            raise FMNPIncorrectDistributionLength(
                'k_frag distribution provided in input data ' +
                'is not the same length as particle size distribution. ' +
                f'Expecting {config["n_size_classes"]}-length array. ' +
                f'Received {len(data["k_frag"])}-length array.'
            )

    # If k_diss isn't a scalar (i.e. we're not calculating
    # a k_diss distribution interally), then check the
    # provided distribution is the correct length
    if not isinstance(validated['k_diss'], (int, float)):
        if len(data['k_diss']) != config['n_size_classes']:
            raise FMNPIncorrectDistributionLength(
                'k_diss distribution provided in input data ' +
                'is not the same length as particle size distribution. ' +
                f'Expecting {config["n_size_classes"]}-length array. ' +
                f'Received {len(data["k_diss"])}-length array.'
            )

    # Check initial conc distribution is the correct length
    if len(data['initial_concs']) != config['n_size_classes']:
        raise FMNPIncorrectDistributionLength(
            'initial_concs distribution provided in input data ' +
            'is not the same length as particle size distribution. ' +
            f'Expecting {config["n_size_classes"]}-length array. ' +
            f'Received {len(data["initial_concs"])}-length array.'
        )

    # TODO extra validation here, e.g. check lengths are n_size_classes
    return validated

