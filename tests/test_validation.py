"""
Tests for the config and data validation
"""
from fragmentmnp.validation import validate_config
from schema import SchemaError


def test_correct_config_validation():
    """
    Test for validating a correct config dict
    """
    correct_config = {
        'n_size_classes': 5,
        'particle_size_range': [0, 10],
        'n_timesteps': 100
    }
    # Validated config will return the config dict if
    # it passes
    validated = validate_config(correct_config)
    assert validated == correct_config


def test_incorrect_config_validation_missing_key():
    """
    Test for validating an incorrect config dict with
    a missing key
    """
    incorrect_config = {}
    # If incorrect_config in invalid, this will raise an error,
    # so catch this and say we've passed the test if this happens
    try:
        validate_config(incorrect_config)
        assert False
    except SchemaError:
        # If an exception has been raised, this test should pass
        assert True


def test_incorrect_config_validation_bad_datatype():
    """
    Test for validating an incorrect config dict with
    incorrect data type
    """
    incorrect_config = {
        'n_size_classes': 5,
        'particle_size_range': 'foo',
        'n_timesteps': 100
    }
    # If incorrect_config in invalid, this will raise an error,
    # so catch this and say we've passed the test if this happens
    try:
        validate_config(incorrect_config)
        assert False
    except SchemaError:
        # If an exception has been raised, this test should pass
        assert True


def test_incorrect_config_too_many_size_classes():
    """
    Test for validating an incorrect config dict with
    too many size classes (over 1000)
    """
    incorrect_config = {
        'n_size_classes': 101,
        'n_timesteps': 100
    }
    # If incorrect_config in invalid, this will raise an error,
    # so catch this and say we've passed the test if this happens
    try:
        validate_config(incorrect_config)
        assert False
    except SchemaError:
        # If an exception has been raised, this test should pass
        assert True


def test_incorrect_config_size_range_not_length_2():
    """
    Test for validating an incorrect config dict with
    particle_size_range not as a length-2 iterable
    """
    incorrect_config = {
        'n_size_classes': 5,
        'particle_size_range': [0, 1, 2],
        'n_timesteps': 100
    }
    # If incorrect_config in invalid, this will raise an error,
    # so catch this and say we've passed the test if this happens
    try:
        validate_config(incorrect_config)
        assert False
    except SchemaError:
        # If an exception has been raised, this test should pass
        assert True
