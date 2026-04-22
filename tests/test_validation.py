"""
Unit tests for the config and data validation
"""
import copy
import pytest
import numpy as np
from schema import SchemaError
from fragmentmnp import FragmentMNP
from fragmentmnp.validation import validate_config, validate_data
import fragmentmnp.examples
from fragmentmnp._errors import FMNPIncorrectDistributionLength, FMNPDistributionValueError
from fragmentmnp.examples import minimal_config, minimal_data


# Get some valid config from the examples module
valid_config = fragmentmnp.examples.full_config
valid_minimal_config = fragmentmnp.examples.minimal_config
valid_data = fragmentmnp.examples.full_data
valid_minimal_data = fragmentmnp.examples.minimal_data


def test_valid_config():
    """
    Test for validating a correct config dict
    """
    validated = validate_config(valid_config)

    for key, value in valid_config.items():
        assert validated[key] == value

    assert "additive_release" in validated
    assert validated["additive_release"] is None


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
    Test for validating a correct data dict.

    NOTE:
    validate_data() returns a dict with defaults filled in. As the model evolves,
    new Optional(..., default=...) entries may appear in the validated dict.

    Therefore, we test that:
      1) All keys provided by the user are preserved exactly.
      2) New additive-related default keys exist and default to None.
    """
    validated = validate_data(valid_data, valid_config)

    # 1) Make sure all original user-provided keys are unchanged
    for k, v in valid_data.items():
        assert k in validated
        assert validated[k] == v

    # 2) Check new collaboration defaults are present
    assert "initial_additive_concs" in validated
    assert validated["initial_additive_concs"] is None

def test_valid_minimal_config_additive_defaults():
        config = copy.deepcopy(valid_minimal_config)
        config["additive_release"] = {
            "model": "analytical"
        }

        validated = validate_config(config)

        assert validated["additive_release"]["model"] == "analytical"
        assert validated["additive_release"]["solver"]["n_terms"] == 50
        assert validated["additive_release"]["solver"]["n_r"] == 60
        assert validated["additive_release"]["solver"]["n_substeps"] == 20
        assert validated["additive_release"]["solver"]["theta"] == 1.0


def test_valid_minimal_data():
    """
    Test the minimal config example passes with defaults filled in
    """
    validated = validate_data(valid_minimal_data, valid_config)
    assert validated['k_diss'] == 0.0

    # Optional: additive defaults should also be present
    assert validated["initial_additive_concs"] is None
    assert validated["additive_release"] is None


def test_invalid_initial_concs_distribution_length():
    """
    Test that inputting an initial_concs distribution that
    isn't the same length as the number of size classes
    results in an error
    """
    invalid_data = copy.deepcopy(valid_minimal_data)
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
    invalid_data = copy.deepcopy(valid_minimal_data)
    invalid_data['initial_concs'][0] = -1
    try:
        validate_data(invalid_data, valid_config)
        assert False
    except SchemaError:
        assert True


def test_k_dist_negative_values():
    """
    Test that specifying k distribution parameters that
    result in negative values returns an error
    """
    invalid_data = copy.deepcopy(valid_minimal_data)
    # Set a negative baseline correction with a constant
    # distribution to make the entire k_frag array negative
    invalid_data['k_frag'] = {'k_f': 0.0, 'k_0': -1.0}
    print(invalid_data)
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
    invalid_data = copy.deepcopy(valid_minimal_data)
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
    valid_config_scalar = copy.deepcopy(valid_minimal_config)
    valid_config_array = copy.deepcopy(valid_minimal_config)
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
    config_arr = copy.deepcopy(valid_minimal_config)
    config_list = copy.deepcopy(valid_minimal_config)
    config_str = copy.deepcopy(valid_minimal_config)
    config_none = copy.deepcopy(valid_minimal_config)
    # Arbitrarily spaced timestep evaluation points
    arr = np.arange(0, config_arr['n_timesteps'], 42)
    config_arr['solver_t_eval'] = arr
    config_list['solver_t_eval'] = list(arr)
    config_str['solver_t_eval'] = 'timesteps'
    config_none['solver_t_eval'] = None

    try:
        validate_config(config_arr)
        validate_config(config_str)
        validate_config(config_none)
        assert True
    except SchemaError:
        assert False

def test_additive_release_analytical_requires_physical_params():
    config = copy.deepcopy(valid_minimal_config)
    config["additive_release"] = {
        "model": "analytical",
        "solver": {"n_terms": 25}
    }

    data = copy.deepcopy(valid_minimal_data)
    data["initial_additive_concs"] = [1.0] * config["n_size_classes"]
    data["additive_release"] = {
        "D_p": 1e-16,
        "K_pw": 1e4,
        # missing D_w on purpose
    }

    try:
        validate_data(data, validate_config(config))
        assert False
    except SchemaError:
        assert True

def test_additive_release_numerical_requires_km_or_Dw():
    config = copy.deepcopy(valid_minimal_config)
    config["additive_release"] = {
        "model": "numerical",
        "solver": {"n_r": 40, "n_substeps": 10, "theta": 1.0}
    }

    data = copy.deepcopy(valid_minimal_data)
    data["initial_additive_concs"] = [1.0] * config["n_size_classes"]
    data["additive_release"] = {
        "D_p": 1e-16,
        "K_pw": 1e4,
        # neither k_m nor D_w provided
    }

    try:
        validate_data(data, validate_config(config))
        assert False
    except SchemaError:
        assert True

def test_config_requires_particle_size_definition():
    config = {
        'n_size_classes': 7,
        'n_timesteps': 100,
    }

    with pytest.raises(SchemaError):
        validate_config(config)

def test_config_rejects_invalid_additive_model():
    config = copy.deepcopy(valid_minimal_config)
    config['additive_release'] = {
        'model': 'foo',
        'solver': {}
    }

    with pytest.raises(SchemaError):
        validate_config(config)

def test_valid_multi_additives():
    config = copy.deepcopy(minimal_config)
    data = copy.deepcopy(minimal_data)

    data["additives"] = [
        {
            "name": "AO168",
            "pools": [
                {
                    "name": "fast",
                    "initial_concs": [0.2] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "solver": {"n_terms": 25},
                        "params": {
                            "D_p": 1e-16,
                            "D_w": 1e-9,
                            "K_pw": 1e4,
                        },
                    },
                },
                {
                    "name": "slow",
                    "initial_concs": [0.8] * config["n_size_classes"],
                    "release": {
                        "model": "numerical",
                        "solver": {
                            "n_r": 40,
                            "n_substeps": 10,
                            "theta": 1.0,
                        },
                        "params": {
                            "D_p": 1e-18,
                            "K_pw": 1e5,
                            "k_m": 1e-8,
                        },
                    },
                },
            ],
        }
    ]

    validated = validate_data(data, validate_config(config))
    assert validated["additives"] is not None
    assert len(validated["additives"]) == 1
    assert len(validated["additives"][0]["pools"]) == 2


def test_invalid_multi_additives_missing_pool_release():
    config = copy.deepcopy(minimal_config)
    data = copy.deepcopy(minimal_data)

    data["additives"] = [
        {
            "name": "AO168",
            "pools": [
                {
                    "name": "fast",
                    "initial_concs": [0.2] * config["n_size_classes"],
                }
            ],
        }
    ]

    with pytest.raises(SchemaError):
        validate_data(data, validate_config(config))


def test_invalid_multi_additives_wrong_initial_conc_length():
    config = copy.deepcopy(minimal_config)
    data = copy.deepcopy(minimal_data)

    data["additives"] = [
        {
            "name": "AO168",
            "pools": [
                {
                    "name": "fast",
                    "initial_concs": [0.2],   # wrong length
                    "release": {
                        "model": "analytical",
                        "params": {
                            "D_p": 1e-16,
                            "D_w": 1e-9,
                            "K_pw": 1e4,
                        },
                    },
                }
            ],
        }
    ]

    with pytest.raises(FMNPIncorrectDistributionLength):
        validate_data(data, validate_config(config))


def test_valid_multi_additives_with_fate_and_transfers():
    config = copy.deepcopy(valid_minimal_config)
    data = copy.deepcopy(valid_minimal_data)
    data["additives"] = [
        {
            "name": "Additive A",
            "pools": [
                {
                    "name": "Pool 1",
                    "initial_concs": [1.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    },
                    "fate": {
                        "k_deg": 0.01,
                        "k_loss": 0.02,
                        "transfers": [{"to": "Pool 2", "k": 0.03}]
                    }
                },
                {
                    "name": "Pool 2",
                    "initial_concs": [0.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    }
                }
            ]
        }
    ]

    validated = validate_data(data, config)
    fate = validated["additives"][0]["pools"][0]["fate"]
    assert fate["k_deg"] == 0.01
    assert fate["k_loss"] == 0.02
    assert fate["transfers"][0]["to"] == "Additive A:Pool 2"


def test_invalid_transfer_target_raises():
    config = copy.deepcopy(valid_minimal_config)
    data = copy.deepcopy(valid_minimal_data)
    data["additives"] = [
        {
            "name": "Additive A",
            "pools": [
                {
                    "name": "Pool 1",
                    "initial_concs": [1.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    },
                    "fate": {
                        "transfers": [{"to": "Missing Pool", "k": 0.03}]
                    }
                }
            ]
        }
    ]

    with pytest.raises(SchemaError):
        validate_data(data, config)

def test_valid_inheritance_block_defaults_to_proportional():
    config = copy.deepcopy(valid_minimal_config)
    data = copy.deepcopy(valid_minimal_data)

    data["additives"] = [
        {
            "name": "Additive A",
            "pools": [
                {
                    "name": "Pool 1",
                    "initial_concs": [1.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    }
                }
            ]
        }
    ]

    validated = validate_data(data, config)
    inh = validated["additives"][0]["pools"][0]["inheritance"]
    assert inh["mode"] == "proportional"


def test_valid_size_biased_inheritance_block():
    config = copy.deepcopy(valid_minimal_config)
    data = copy.deepcopy(valid_minimal_data)

    data["additives"] = [
        {
            "name": "Additive A",
            "pools": [
                {
                    "name": "Pool 1",
                    "initial_concs": [1.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    },
                    "inheritance": {
                        "mode": "size_biased",
                        "beta": -1.0
                    }
                }
            ]
        }
    ]

    validated = validate_data(data, config)
    inh = validated["additives"][0]["pools"][0]["inheritance"]
    assert inh["mode"] == "size_biased"
    assert inh["beta"] == -1.0


def test_invalid_inheritance_mode_raises():
    config = copy.deepcopy(valid_minimal_config)
    data = copy.deepcopy(valid_minimal_data)

    data["additives"] = [
        {
            "name": "Additive A",
            "pools": [
                {
                    "name": "Pool 1",
                    "initial_concs": [1.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    },
                    "inheritance": {
                        "mode": "not_a_mode"
                    }
                }
            ]
        }
    ]

    with pytest.raises(SchemaError):
        validate_data(data, config)


def test_valid_particulate_fate_allows_size_dependent_rates():
    config = copy.deepcopy(valid_minimal_config)
    data = copy.deepcopy(valid_minimal_data)

    data["additives"] = [
        {
            "name": "Additive A",
            "pools": [
                {
                    "name": "Pool 1",
                    "initial_concs": [1.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    },
                    "fate": {
                        "k_deg": [1e-5] * config["n_size_classes"],
                        "k_loss": [0.0] * config["n_size_classes"],
                    }
                },
                {
                    "name": "Pool 2",
                    "initial_concs": [0.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    }
                }
            ]
        }
    ]

    validated = validate_data(data, config)
    assert isinstance(validated["additives"][0]["pools"][0]["fate"]["k_deg"], list)

def test_invalid_size_dependent_fate_wrong_length_raises():
    config = copy.deepcopy(valid_minimal_config)
    data = copy.deepcopy(valid_minimal_data)

    data["additives"] = [
        {
            "name": "Additive A",
            "pools": [
                {
                    "name": "Pool 1",
                    "initial_concs": [1.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    },
                    "fate": {
                        "k_deg": [1e-5, 2e-5]
                    }
                }
            ]
        }
    ]

    with pytest.raises(SchemaError):
        validate_data(data, config)

def test_invalid_negative_fate_rate_raises():
    config = copy.deepcopy(valid_minimal_config)
    data = copy.deepcopy(valid_minimal_data)
    data["additives"] = [
        {
            "name": "Additive A",
            "pools": [
                {
                    "name": "Pool 1",
                    "initial_concs": [1.0] * config["n_size_classes"],
                    "release": {
                        "model": "analytical",
                        "params": {"D_p": 1e-16, "D_w": 1e-9, "K_pw": 1e4}
                    },
                    "fate": {
                        "k_deg": -0.01
                    }
                }
            ]
        }
    ]

    with pytest.raises(SchemaError):
        validate_data(data, config)
