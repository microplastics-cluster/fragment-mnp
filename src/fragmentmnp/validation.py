"""
Validation of config and data (:mod:`fragmentmnp.validation`)
=============================================================

Provides config and input data validation for the FRAGMENT-MNP model.

This file is extended to optionally validate:
- initial_chemical_concs (array, length n_size_classes, >=0)
- chemical_release settings split across:
    config: model + solver controls
    data:   physical/material parameters
"""
import numpy as np
from schema import And, Optional, Or, Schema, SchemaError

from ._errors import FMNPIncorrectDistributionLength


def _is_positive_array(arr):
    """Check if arr is iterable and all elements are positive."""
    is_array = True
    try:
        _ = iter(arr)
        nparr = np.array(arr)
        if np.any(nparr < 0.0):
            is_array = False
    except TypeError:
        is_array = False
    return is_array


def _is_array(arr):
    """Check if arr is iterable."""
    is_array = True
    try:
        _ = iter(arr)
    except TypeError:
        is_array = False
    return is_array


particle_size_range_schema = And(
    Or((int, float), [int, float]),
    And(
        lambda d: len(d) == 2,
        error='particle_size_range must be a length-2 iterable'
    )
)


def k_dist_schema(dims):
    k_dist_schema_ = Or(
        Or(
            And(int, lambda x: x >= 0.0),
            And(float, lambda x: x >= 0.0)
        ),
        {
            'k_f': And(Or(int, float), lambda x: x >= 0.0),
            Optional('k_0', default=0.0): Or(int, float),
            Optional('is_compound', default=True): bool,
            **{
                Optional(f'{name}_{x}'): Or(int, float)
                for x in dims
                for name in ['alpha', 'B', 'beta', 'gamma', 'delta1']
            },
            **{
                Optional(f'{name}_{x}'): Or(int, float, None)
                for x in dims
                for name in ['C', 'D', 'delta2']
            },
            **{
                Optional(f'A_{x}'): Or(int, float, _is_array)
                for x in dims
            },
        }
    )
    return k_dist_schema_


k_dist_2d_schema = k_dist_schema(['t', 's'])
k_dist_t_schema = k_dist_schema(['t'])


def _normalize_chemical_release_input(data: dict, config: dict) -> tuple[dict, dict]:
    """
    Backward compatibility for old additive/chemical release schemas.

    Old:
        data['additive_release'] = {
            'model': 'analytical'|'numerical',
            'params': {...}
        }

    New:
        config['chemical_release'] = {
            'model': ...,
            'solver': {...}
        }
        data['chemical_release'] = {
            'D_p': ..., 'D_w': ..., 'K_pw': ..., 'k_m': ...
        }
    """
    data = dict(data)
    config = dict(config)

    # Old -> new key aliasing
    if 'initial_additive_concs' in data and 'initial_chemical_concs' not in data:
        data['initial_chemical_concs'] = data.pop('initial_additive_concs')

    if 'additive_release' in config and 'chemical_release' not in config:
        config['chemical_release'] = config.pop('additive_release')

    if 'additive_release' in data and 'chemical_release' not in data:
        data['chemical_release'] = data.pop('additive_release')

    add_data = data.get('chemical_release', None)

    # Backward compatibility for old all-in-data schema:
    # data['additive_release'/'chemical_release'] = {'model': ..., 'params': {...}}
    if isinstance(add_data, dict) and ('model' in add_data or 'params' in add_data):
        old_model = str(add_data.get('model', 'analytical')).lower()
        old_params = dict(add_data.get('params', {}))

        # Split solver vs physical params
        solver = {}
        for key in ['n_r', 'n_substeps', 'theta', 'n_terms']:
            if key in old_params:
                solver[key] = old_params.pop(key)

        # Only set config chemical_release if not already explicitly set
        if config.get('chemical_release', None) is None:
            config['chemical_release'] = {
                'model': old_model,
                'solver': solver,
            }

        data['chemical_release'] = old_params

    return data, config


def chemical_release_config_schema():
    """
    Config-side chemical release settings:
      - model selection
      - numerical solver controls
    """
    return {
        Optional('model', default='analytical'): And(
            str,
            lambda x: x.lower() in ['analytical', 'numerical']
        ),
        Optional('solver', default={}): {
            Optional('n_r', default=60): And(int, lambda x: x > 2),
            Optional('n_substeps', default=20): And(int, lambda x: x >= 1),
            Optional('theta', default=1.0): And(
                Or(int, float),
                lambda x: 0.0 <= x <= 1.0
            ),
            Optional('n_terms', default=50): And(int, lambda x: x >= 1),
        }
    }


def chemical_release_data_schema():
    """
    Data-side chemical release parameters:
      - physical/material parameters only

    Cross-field requirements depending on model are checked in validate_data().
    """
    return {
        Optional('D_p'): And(Or(int, float), lambda x: x >= 0.0),
        Optional('D_w'): And(Or(int, float), lambda x: x >= 0.0),
        Optional('K_pw'): And(Or(int, float), lambda x: x > 0.0),
        Optional('k_m'): And(Or(int, float), lambda x: x >= 0.0),
    }


config_schema = Schema({
    'n_size_classes': And(int, lambda d: d <= 100),
    Optional('particle_size_range'): particle_size_range_schema,
    Optional('particle_size_classes'): _is_positive_array,
    'n_timesteps': int,
    Optional('dt', default=1): int,
    Optional('solver_method', default='LSODA'): str,
    Optional('solver_atol', default=1e-6): Or(float, [float]),
    Optional('solver_rtol', default=1e-3): float,
    Optional('solver_max_step', default=np.inf): float,
    Optional('solver_t_eval', default='timesteps'): Or(
        _is_positive_array,
        'timesteps',
        None
    ),
    Optional('chemical_release', default=None): Or(
        None,
        chemical_release_config_schema()
    ),
})


data_schema = Schema({
    'initial_concs': _is_positive_array,
    Optional('initial_concs_diss', default=0.0): Or(float, int),
    'density': And(Or(int, float), lambda x: x >= 0.0),
    'k_frag': k_dist_2d_schema,
    Optional('k_diss', default=0.0): k_dist_2d_schema,
    Optional('k_min', default=0.0): k_dist_t_schema,
    Optional('fsd_beta', default=0.0): Or(int, float),

    Optional('initial_chemical_concs', default=None): Or(None, _is_positive_array),

    Optional('chemical_release', default=None): Or(
        None,
        chemical_release_data_schema()
    ),
})


def _validate_config_cross_checks(config: dict) -> dict:
    """
    Additional config validation that depends on combinations of fields,
    beyond what the schema package can express cleanly.
    """
    has_ps_classes = 'particle_size_classes' in config
    has_ps_range = 'particle_size_range' in config

    if not (has_ps_classes or has_ps_range):
        raise SchemaError(
            "Model config must contain either 'particle_size_classes' "
            "or 'particle_size_range'."
        )

    chem_cfg = config.get('chemical_release', None)
    if chem_cfg is not None:
        model = str(chem_cfg.get('model', 'analytical')).lower()
        if model not in ['analytical', 'numerical']:
            raise SchemaError(
                "config.chemical_release.model must be one of "
                "['analytical', 'numerical']."
            )

    return config


def _fill_chemical_release_defaults(config: dict) -> dict:
    """
    Fill nested defaults for chemical_release.solver, because nested schema
    defaults are not always fully materialized when the parent dict is present
    but child keys are omitted.
    """
    chem_cfg = config.get('chemical_release', None)
    if chem_cfg is None:
        return config

    solver = dict(chem_cfg.get('solver', {}))
    solver.setdefault('n_r', 60)
    solver.setdefault('n_substeps', 20)
    solver.setdefault('theta', 1.0)
    solver.setdefault('n_terms', 50)

    chem_cfg['solver'] = solver
    config['chemical_release'] = chem_cfg
    return config


def _validate_chemical_release_cross_checks(data: dict, config: dict) -> None:
    """
    Cross-validation of chemical_release config/data split.
    """
    chem_cfg = config.get('chemical_release', None)
    chem_data = data.get('chemical_release', None)

    if chem_cfg is None:
        return

    model = str(chem_cfg.get('model', 'analytical')).lower()

    if chem_data is None:
        raise SchemaError(
            f"config.chemical_release is set for model '{model}', "
            "but data.chemical_release is missing."
        )

    if model == 'analytical':
        required = ['D_p', 'D_w', 'K_pw']
        missing = [k for k in required if k not in chem_data]
        if missing:
            raise SchemaError(
                "data.chemical_release missing required keys for analytical model: "
                f"{missing}"
            )

    elif model == 'numerical':
        required = ['D_p', 'K_pw']
        missing = [k for k in required if k not in chem_data]
        if missing:
            raise SchemaError(
                "data.chemical_release missing required keys for numerical model: "
                f"{missing}"
            )

        if ('k_m' not in chem_data) and ('D_w' not in chem_data):
            raise SchemaError(
                "data.chemical_release for numerical model must contain either "
                "'k_m' or 'D_w'."
            )


def validate_config(config: dict) -> dict:
    # Normalize old config-side names before schema validation
    config = dict(config)

    if 'additive_release' in config and 'chemical_release' not in config:
        config['chemical_release'] = config.pop('additive_release')

    validated = Schema(config_schema).validate(config)
    validated = _validate_config_cross_checks(validated)
    validated = _fill_chemical_release_defaults(validated)

    # Backward-compatible alias in validated output
    if 'chemical_release' in validated and 'additive_release' not in validated:
        validated['additive_release'] = validated['chemical_release']

    return validated


def validate_data(data: dict, config: dict) -> dict:
    data, config = _normalize_chemical_release_input(data, config)
    validated = Schema(data_schema).validate(data)

    if len(validated['initial_concs']) != config['n_size_classes']:
        raise FMNPIncorrectDistributionLength(
            'initial_concs distribution provided in input data '
            'is not the same length as particle size distribution. '
            f'Expecting {config["n_size_classes"]}-length array. '
            f'Received {len(validated["initial_concs"])}-length array.'
        )

    if validated.get('initial_chemical_concs') is not None:
        if len(validated['initial_chemical_concs']) != config['n_size_classes']:
            raise FMNPIncorrectDistributionLength(
                'initial_chemical_concs distribution provided in input data '
                'is not the same length as particle size distribution. '
                f'Expecting {config["n_size_classes"]}-length array. '
                f'Received {len(validated["initial_chemical_concs"])}-length array.'
            )

    _validate_chemical_release_cross_checks(validated, config)

    # Backward-compatible aliases in validated output
    if 'initial_chemical_concs' in validated and 'initial_additive_concs' not in validated:
        validated['initial_additive_concs'] = validated['initial_chemical_concs']

    if 'chemical_release' in validated and 'additive_release' not in validated:
        validated['additive_release'] = validated['chemical_release']

    return validated