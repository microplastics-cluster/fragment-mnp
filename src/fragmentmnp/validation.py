"""
Validation of config and data (:mod:`fragmentmnp.validation`)
=============================================================

Provides config and input data validation for the FRAGMENT-MNP model
"""
import numpy as np
from schema import And, Optional, Or, Schema

from ._errors import FMNPIncorrectDistributionLength


def _is_positive_array(arr):
    """Check if arr is iterable and all elements are positive"""
    is_array = True
    try:
        # Check it's iterable
        _ = iter(arr)
        # Check if any element is < 0
        nparr = np.array(arr)
        if np.any(nparr < 0.0):
            is_array = False
    except TypeError:
        is_array = False
    return is_array


def _is_array(arr):
    """Check if arr is iterable"""
    is_array = True
    try:
        # Check it's iterable
        _ = iter(arr)
    except TypeError:
        is_array = False
    return is_array


# The schema that particle size ranges should follow
particle_size_range_schema = And(Or((int, float), [int, float]),
                                 And(lambda d: len(d) == 2,
                                     error='particle_size_range must ' +
                                           'be a length-2 iterable'))


# The schema that rate constant 2D (time and surface area)
# distributions, like k_frag and k_diss, should follow. Either
# a scalar is given (and it is treated as constant), or
# a dict is given with the params required to calculate the 2D
# distribution. Basic checks here ensure k values given are
# greater than zero, and auditing during the distribution
# calculate makes sure no values in the calculated distribution
# are less than zero
k_dist_2d_schema = Or(
    Or(And(int, lambda x: x >= 0.0),
       And(float, lambda x: x >= 0.0)),
    {
        'k_f': And(Or(int, float), lambda x: x >= 0.0),
        Optional('k_0', default=0.0): Or(int, float),
        Optional('is_compound', default=True): bool,
        **{
            Optional(f'{name}_{x}'): Or(int, float)
            for x in ['t', 's']
            for name in ['alpha', 'B', 'beta', 'gamma',
                         'delta1']
        },
        **{
            Optional(f'{name}_{x}'): Or(int, float, None)
            for x in ['t', 's']
            for name in ['C', 'D', 'delta2']
        },
        **{
            Optional(f'A_{x}'): Or(int, float, _is_array)
            for x in ['t', 's']
        },
    }
)


# The schema that the config dict should follow
config_schema = Schema({
    # There should be <= 100 size classes
    'n_size_classes': And(int, lambda d: d <= 100),
    # Size range should be a length-2 iterable of type int or float
    Optional('particle_size_range'): particle_size_range_schema,
    # Size classes should be a list of ints of floats
    Optional('particle_size_classes'): _is_positive_array,
    # Timesteps should be an integer
    'n_timesteps': int,
    # Length of timesteps should be an integer (unit of seconds)
    Optional('dt', default=1): int,
    # What ODE solver method should be used? Should be one of
    # those available in scipy.solve_ivp:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
    Optional('solver_method', default='LSODA'): str,
    # Error tolerances for the ODE solver
    Optional('solver_atol', default=1e-6): Or(float, [float]),
    Optional('solver_rtol', default=1e-3): float,
    # Max step size for the ODE solver
    Optional('solver_max_step', default=np.inf): float,
    Optional('solver_t_eval', default='timesteps'): Or(_is_positive_array,
                                                       'timesteps',
                                                       None)
})


# The schema that the data dict should follow
data_schema = Schema({
    # Initial concs must be a list and >= 0
    'initial_concs': _is_positive_array,
    # Initial dissolved fraction concentration
    Optional('initial_concs_diss', default=0.0): Or(float, int),
    # Density must either be a float/int and greater than 0
    'density': And(Or(int, float), lambda x: x >= 0.0),
    # k_frag must either be a float/int, or a dict containing
    # and average value and parameters to create distribution from.
    # Both default to zero
    'k_frag': k_dist_2d_schema,
    Optional('k_diss', default=0.0): k_dist_2d_schema,
    # fsd_beta is an empirical param that scales the depedence
    # of the fragment size distribution on particle diameter d
    # accordingly to d^beta. beta=0 means an equal split
    Optional('fsd_beta', default=0.0): Or(int, float)
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
