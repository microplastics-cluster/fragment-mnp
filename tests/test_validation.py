"""
Unit tests for the config and data validation
"""
import numpy as np
from schema import SchemaError
from fragmentmnp import FragmentMNP
from fragmentmnp.validation import validate_config, validate_data
import fragmentmnp.examples
from fragmentmnp._errors import FMNPIncorrectDistributionLength, FMNPDistributionValueError


# Get some valid config from the examples module
valid_config = fragmentmnp.examples.full_config
valid_minimal_config = fragmentmnp.examples.minimal_config
valid_data = fragmentmnp.examples.full_data
valid_minimal_data = fragmentmnp.examples.minimal_data


def test_valid_config():
    """
    Test for validating a correct config dict
    """
    # Validated config will return the config dict if
    # it passes
    validated = validate_config(valid_config)
    assert validated == valid_config


def test_valid_minimal_config():
    """
    Test the minimal config example passes with
    defaults filled in
    """
    validated = validate_config(valid_minimal_config)
    # dt should have been defaulted to 1
    assert validated['dt'] == 1


def test_invalid_config_missing_key():
    """
    Test for validating an incorrect config dict with
    a missing key
    """
    # Take a copy of the valid config and remove the n_size_classes key
    invalid_config = valid_config.copy()
    del invalid_config['n_size_classes']
    # If incorrect_config in invalid, this will raise an error,
    # so catch this and say we've passed the test if this happens
    try:
        validate_config(invalid_config)
        assert False
    except SchemaError:
        # If an exception has been raised, this test should pass
        assert True


def test_invalid_config_bad_datatype():
    """
    Test for validating an incorrect config dict with
    incorrect data type
    """
    invalid_config = valid_config | {'particle_size_range': 'foo'}
    # If incorrect_config in invalid, this will raise an error,
    # so catch this and say we've passed the test if this happens
    try:
        validate_config(invalid_config)
        assert False
    except SchemaError:
        # If an exception has been raised, this test should pass
        assert True


def test_invalid_config_too_many_size_classes():
    """
    Test for validating an incorrect config dict with
    too many size classes (over 1000)
    """
    invalid_config = valid_config | {'n_size_classes': 101}
    # If incorrect_config in invalid, this will raise an error,
    # so catch this and say we've passed the test if this happens
    try:
        validate_config(invalid_config)
        assert False
    except SchemaError:
        # If an exception has been raised, this test should pass
        assert True


def test_invalid_config_size_range_not_length_2():
    """
    Test for validating an incorrect config dict with
    particle_size_range not as a length-2 iterable
    """
    invalid_config = valid_config | {'particle_size_range': [0, 1, 2]}
    # If incorrect_config in invalid, this will raise an error,
    # so catch this and say we've passed the test if this happens
    try:
        validate_config(invalid_config)
        assert False
    except SchemaError:
        # If an exception has been raised, this test should pass
        assert True


def test_valid_data():
    """
    Test for validating a correct data dict
    """
    validated = validate_data(valid_data, valid_config)
    assert validated == valid_data


def test_valid_minimal_data():
    """
    Test the minimal config example passes with
    defaults filled in
    """
    validated = validate_data(valid_minimal_data, valid_config)
    # theta_1 should have been defaulted to 0
    assert validated['theta_1'] == 0.0



def test_invalid_k_frag_distribution_length():
    """
    Test that inputting a k_frag distribution that
    isn't the same length as the number of size classes
    results in an error
    """
    invalid_data = valid_minimal_data.copy()
    invalid_data['k_frag'] = [1]
    try:
        validate_data(invalid_data, valid_config)
        assert False
    except FMNPIncorrectDistributionLength:
        assert True


def test_invalid_initial_concs_distribution_length():
    """
    Test that inputting an initial_concs distribution that
    isn't the same length as the number of size classes
    results in an error
    """
    invalid_data = valid_minimal_data.copy()
    invalid_data['initial_concs'] = [1]
    try:
        validate_data(invalid_data, valid_config)
        assert False
    except FMNPIncorrectDistributionLength:
        assert True


def test_negative_initial_concs():
    """
    Test that inputting an initial_concs distribution that
    has a negative value results in an error
    """
    invalid_data = valid_minimal_data.copy()
    invalid_data['initial_concs'][0] = -1
    try:
        validate_data(invalid_data, valid_config)
        assert False
    except SchemaError:
        assert True


def test_k_distribution_negative_values():
    """
    Test that specifying k distribution parameters that
    result in negative values returns an error
    """
    invalid_data = valid_minimal_data.copy()
    # Set a negative baseline correction with a constant
    # distribution to make the entire k_frag array negative
    invalid_data['k_frag'] = {'k_f': 0.0, 'k_1': -1.0}
    try:
        _ = FragmentMNP(valid_minimal_config, invalid_data)
        assert False
    except FMNPDistributionValueError:
        assert True


def test_noniterable_array():
    """
    Test that inputting an initial_concs distribution that
    isn't iterable results in an error
    """
    invalid_data = valid_minimal_data.copy()
    invalid_data['initial_concs'] = 1
    try:
        validate_data(invalid_data, valid_config)
        assert False
    except SchemaError:
        assert True


def test_atol_array_and_scalar():
    """
    Test that atol can be input as an array or a scalar
    """
    valid_config_scalar = valid_minimal_config.copy()
    valid_config_array = valid_minimal_config.copy()
    valid_config_scalar['solver_atol'] = 1e-6
    valid_config_array['solver_atol'] = [1e-6] * 7
    try:
        validate_config(valid_config_scalar)
        validate_config(valid_config_array)
        assert True
    except SchemaError:
        assert False


def test_t_eval():
    """
    Test that solver_t_eval can be input as an array
    """
    config_arr = valid_minimal_config.copy()
    config_list = valid_minimal_config.copy()
    config_str = valid_minimal_config.copy()
    config_none = valid_minimal_config.copy()
    # Arbitrarily spaced timestep evaluation points
    arr = np.arange(0, config_arr['n_timesteps'], 42)
    config_arr['solver_t_eval'] = arr
    config_list['solver_t_eval'] = list(arr)
    config_str['solver_t_eval'] = 'integer'
    config_none['solver_t_eval'] = None

    try:
        validate_config(config_arr)
        validate_config(config_str)
        validate_config(config_none)
        assert True
    except SchemaError:
        assert False
