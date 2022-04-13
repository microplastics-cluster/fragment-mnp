"""
Tests for the config and data validation
"""
from fragmentmnp.validation import validate_config
import fragmentmnp.examples
from schema import SchemaError


# Get some valid config from the examples module
valid_config = fragmentmnp.examples.config


def test_correct_config_validation():
    """
    Test for validating a correct config dict
    """
    # Validated config will return the config dict if
    # it passes
    validated = validate_config(valid_config)
    assert validated == valid_config


def test_incorrect_config_validation_missing_key():
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


def test_incorrect_config_validation_bad_datatype():
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


def test_incorrect_config_too_many_size_classes():
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


def test_incorrect_config_size_range_not_length_2():
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
