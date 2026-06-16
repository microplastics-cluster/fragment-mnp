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


def _is_release_parameter_like(x):
    """
    Accept release parameters as constants, size vectors, or time-dependent
    dictionaries. Full cross-checking is done after we know the model type and
    number of size classes.
    """
    if isinstance(x, (int, float)):
        return True
    if isinstance(x, dict):
        return True
    try:
        arr = np.asarray(x, dtype=float)
    except Exception:
        return False
    return arr.ndim in (1, 2)


def chemical_release_data_schema():
    """
    Data-side chemical release parameters.

    Each physical parameter can now be:
      - a scalar constant, e.g. 1e-16
      - a size-class vector, e.g. [1e-16, ...]
      - a time-series dict, e.g. {"times": [...], "values": [...]}
      - a simple model dict, e.g. {"type": "linear", "initial": 1e-16, "rate": 1e-3}
      - a FRAGMENT-MNP distribution dict with k_f/k_0 and t/s regression terms

    Cross-field requirements depending on model are checked in validate_data().
    """
    return {
        Optional('D_p'): _is_release_parameter_like,
        Optional('D_w'): _is_release_parameter_like,
        Optional('K_pw'): _is_release_parameter_like,
        Optional('k_m'): _is_release_parameter_like,
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
    Optional('additives', default=None): Or(None, list),
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


def _validate_time_series_release_parameter(name: str,
                                            value: dict,
                                            n_size_classes: int,
                                            strictly_positive: bool) -> None:
    """Validate {times, values} release-parameter dictionaries."""
    if 'times' not in value or 'values' not in value:
        raise SchemaError(
            f"Release parameter '{name}' time-series dict must contain 'times' and 'values'."
        )
    times = np.asarray(value['times'], dtype=float)
    values = np.asarray(value['values'], dtype=float)

    if times.ndim != 1 or times.size < 2:
        raise SchemaError(f"Release parameter '{name}' times must be a 1D array with at least two values.")
    if np.any(~np.isfinite(times)) or np.any(np.diff(times) <= 0.0):
        raise SchemaError(f"Release parameter '{name}' times must be finite and strictly increasing.")
    if np.any(~np.isfinite(values)):
        raise SchemaError(f"Release parameter '{name}' values must be finite.")

    if values.ndim == 1:
        if values.size != times.size:
            raise SchemaError(
                f"Release parameter '{name}' values must have the same length as times."
            )
    elif values.ndim == 2:
        valid_time_by_size = values.shape == (times.size, n_size_classes)
        valid_size_by_time = values.shape == (n_size_classes, times.size)
        if not (valid_time_by_size or valid_size_by_time):
            raise SchemaError(
                f"Release parameter '{name}' 2D values must have shape "
                f"(len(times), n_size_classes) or (n_size_classes, len(times))."
            )
    else:
        raise SchemaError(f"Release parameter '{name}' values must be 1D or 2D.")

    if strictly_positive:
        ok = np.all(values > 0.0)
        msg = f"Release parameter '{name}' values must be positive."
    else:
        ok = np.all(values >= 0.0)
        msg = f"Release parameter '{name}' values must be non-negative."
    if not ok:
        raise SchemaError(msg)


def _validate_release_parameter(name: str,
                                value,
                                n_size_classes: int,
                                strictly_positive: bool = True) -> None:
    """
    Validate one physical release parameter while preserving the user's input.

    Accepted forms:
      scalar constant
      length-n_size_classes vector
      {times, values} time series, where values can be scalar over size or size-resolved
      {type/model: constant|linear|exponential|logistic, ...}
      k-distribution style dict with k_f/k_0 and optional t/s terms
    """
    def _check_values(arr, what):
        arr = np.asarray(arr, dtype=float)
        if np.any(~np.isfinite(arr)):
            raise SchemaError(f"Release parameter '{name}' {what} must be finite.")
        if strictly_positive:
            if np.any(arr <= 0.0):
                raise SchemaError(f"Release parameter '{name}' {what} must be positive.")
        else:
            if np.any(arr < 0.0):
                raise SchemaError(f"Release parameter '{name}' {what} must be non-negative.")

    if isinstance(value, (int, float)):
        _check_values([value], 'value')
        return

    if isinstance(value, dict):
        # Time-series interpolation form.
        if 'times' in value or 'values' in value:
            _validate_time_series_release_parameter(
                name=name,
                value=value,
                n_size_classes=n_size_classes,
                strictly_positive=strictly_positive,
            )
            return

        # Existing FRAGMENT-MNP distribution style, reusing k_f/k_0 with t/s regressions.
        if 'k_f' in value:
            _check_values([value.get('k_f', 0.0)], 'k_f')
            if 'k_0' in value:
                _check_values([value.get('k_0')], 'k_0')
            return

        # Simple parametric ageing/weathering models.
        model = str(value.get('type', value.get('model', 'constant'))).lower()
        valid_models = {'constant', 'linear', 'exponential', 'logistic'}
        if model not in valid_models:
            raise SchemaError(
                f"Release parameter '{name}' model/type must be one of {sorted(valid_models)}."
            )

        base = value.get('initial', value.get('base', value.get('value', None)))
        if base is None:
            raise SchemaError(
                f"Release parameter '{name}' model dict must contain one of 'initial', 'base', or 'value'."
            )
        _check_values(base, 'base/initial/value')

        if 'size_factor' in value:
            sf = np.asarray(value['size_factor'], dtype=float)
            if sf.ndim != 1 or sf.size != n_size_classes:
                raise SchemaError(
                    f"Release parameter '{name}' size_factor must be a length-{n_size_classes} array."
                )
            _check_values(sf, 'size_factor')
        return

    try:
        arr = np.asarray(value, dtype=float)
    except Exception as exc:
        raise SchemaError(f"Release parameter '{name}' has unsupported type.") from exc

    if arr.ndim == 1:
        if arr.size != n_size_classes:
            raise SchemaError(
                f"Release parameter '{name}' vector must have length {n_size_classes}; "
                "use {'times': ..., 'values': ...} for time series."
            )
        _check_values(arr, 'vector')
        return

    if arr.ndim == 2:
        # This is accepted as a compact (time, size) matrix for advanced users.
        if arr.shape[1] != n_size_classes and arr.shape[0] != n_size_classes:
            raise SchemaError(
                f"Release parameter '{name}' 2D array must have one dimension equal to n_size_classes."
            )
        _check_values(arr, 'matrix')
        return

    raise SchemaError(f"Release parameter '{name}' must be scalar, vector, matrix, or dict.")


def _validate_release_params_for_model(params: dict, model: str, n_size_classes: int) -> None:
    """Validate required release parameters for analytical and numerical models."""
    model = str(model).lower()
    if model == 'analytical':
        required = ['D_p', 'D_w', 'K_pw']
    elif model == 'numerical':
        required = ['D_p', 'K_pw']
    else:
        raise SchemaError(f"Unknown release model '{model}'.")

    missing = [k for k in required if k not in params]
    if missing:
        label = 'Analytical' if model == 'analytical' else 'Numerical'
        raise SchemaError(f"{label} release is missing required params: {missing}")

    if model == 'numerical' and ('k_m' not in params) and ('D_w' not in params):
        raise SchemaError("Numerical release must contain either 'k_m' or 'D_w'.")

    # D_p, D_w, and K_pw must be positive when used. k_m may be zero.
    for key in ['D_p', 'D_w', 'K_pw']:
        if key in params:
            _validate_release_parameter(key, params[key], n_size_classes, strictly_positive=True)
    if 'k_m' in params:
        _validate_release_parameter('k_m', params['k_m'], n_size_classes, strictly_positive=False)


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

    _validate_release_params_for_model(
        params=chem_data,
        model=model,
        n_size_classes=config['n_size_classes'],
    )


def _default_release_solver(model: str) -> dict:
    model = str(model).lower()
    if model == 'analytical':
        return {'n_terms': 50}
    if model == 'numerical':
        return {'n_r': 60, 'n_substeps': 20, 'theta': 1.0}
    raise SchemaError(f"Unknown release model '{model}'.")


def _validate_release_block(release: dict, n_size_classes: int) -> dict:
    if not isinstance(release, dict):
        raise SchemaError("Each pool.release must be a dict.")

    model = str(release.get('model', 'analytical')).lower()
    if model not in ['analytical', 'numerical']:
        raise SchemaError("pool.release.model must be 'analytical' or 'numerical'.")

    solver = dict(_default_release_solver(model))
    solver.update(release.get('solver', {}))
    params = dict(release.get('params', {}))

    _validate_release_params_for_model(
        params=params,
        model=model,
        n_size_classes=n_size_classes,
    )

    return {
        'model': model,
        'solver': solver,
        'params': params
    }

def _validate_inheritance_block(inheritance: dict | None) -> dict:
    """
    Validate particulate fragmentation inheritance settings.

    Supported modes
    ---------------
    proportional
        Daughter fragments inherit additive in proportion to the parent
        additive/polymer concentration. This is the current behaviour.

    size_biased
        Daughter fragments inherit additive using a weighting proportional
        to daughter diameter**beta. beta < 0 enriches smaller daughters,
        beta > 0 enriches larger daughters.

    surface_enriched
        Daughter fragments inherit additive using daughter surface area
        weighting (surface area ** gamma). gamma=1 corresponds to direct
        surface-area weighting.

    Notes
    -----
    These rules redistribute the additive mass released by fragmentation
    across daughter size classes, while conserving additive mass.
    """
    if inheritance is None:
        inheritance = {}

    if not isinstance(inheritance, dict):
        raise SchemaError("pool.inheritance must be a dict if provided.")

    mode = str(inheritance.get("mode", "proportional")).lower()
    valid_modes = ["proportional", "size_biased", "surface_enriched"]
    if mode not in valid_modes:
        raise SchemaError(
            f"pool.inheritance.mode must be one of {valid_modes}."
        )

    out = {"mode": mode}

    if mode == "size_biased":
        beta = float(inheritance.get("beta", 0.0))
        out["beta"] = beta

    if mode == "surface_enriched":
        gamma = float(inheritance.get("gamma", 1.0))
        if gamma < 0.0:
            raise SchemaError("pool.inheritance.gamma must be non-negative.")
        out["gamma"] = gamma

    return out

def _is_nonnegative_scalar_or_array(x, n_expected=None):
    """
    Accept either:
    - non-negative scalar
    - iterable of non-negative values

    If n_expected is provided and x is iterable, require len(x) == n_expected.
    """
    if isinstance(x, (int, float)):
        return float(x) >= 0.0

    try:
        arr = np.array(x, dtype=float)
    except Exception:
        return False

    if arr.ndim != 1:
        return False
    if np.any(arr < 0.0):
        return False
    if n_expected is not None and len(arr) != n_expected:
        return False
    return True

def _validate_additives_structure(additives: list, n_size_classes: int) -> list:
    if not isinstance(additives, list) or len(additives) == 0:
        raise SchemaError("data.additives must be a non-empty list.")

    out = []
    for additive in additives:
        if not isinstance(additive, dict):
            raise SchemaError("Each additive entry must be a dict.")
        if 'name' not in additive:
            raise SchemaError("Each additive must contain 'name'.")
        if 'pools' not in additive:
            raise SchemaError(f"Additive '{additive['name']}' must contain 'pools'.")

        additive_name = additive['name']
        pools = additive['pools']
        if not isinstance(pools, list) or len(pools) == 0:
            raise SchemaError(
                f"Additive '{additive_name}' must contain a non-empty pools list."
            )

        medium_pools = list(additive.get('medium_pools', []))
        if len(medium_pools) == 0:
            medium_pools = [{'name': 'medium', 'initial_mass': 0.0, 'fate': {}}]

        valid_particulate_targets = {f"{additive_name}:{pool['name']}" for pool in pools}
        valid_medium_targets = {f"{additive_name}:{pool['name']}" for pool in medium_pools}

        pools_out = []
        for pool in pools:
            if not isinstance(pool, dict):
                raise SchemaError(
                    f"Each pool in additive '{additive_name}' must be a dict."
                )
            if 'name' not in pool:
                raise SchemaError(
                    f"Each pool in additive '{additive_name}' must contain 'name'."
                )
            if 'initial_concs' not in pool:
                raise SchemaError(
                    f"Pool '{pool.get('name', '?')}' in additive '{additive_name}' must contain 'initial_concs'."
                )
            if not _is_positive_array(pool['initial_concs']):
                raise SchemaError(f"Pool '{pool['name']}' initial_concs must be a non-negative array.")
            if len(pool['initial_concs']) != n_size_classes:
                raise FMNPIncorrectDistributionLength(
                    f"Pool '{pool['name']}' in additive '{additive_name}' has "
                    f"{len(pool['initial_concs'])} initial concentrations; expected {n_size_classes}."
                )
            if 'release' not in pool:
                raise SchemaError(
                    f"Pool '{pool['name']}' in additive '{additive_name}' must contain 'release'."
                )

            release = _validate_release_block(pool['release'], n_size_classes=n_size_classes)
            release_target = _normalize_transfer_target(
                pool.get('release', {}).get('target', f"{additive_name}:medium"), additive_name
            )
            if release_target not in valid_medium_targets:
                raise SchemaError(
                    f"Unknown release target '{release_target}' for particulate pool '{additive_name}:{pool['name']}'."
                )
            release['target'] = release_target

            pools_out.append({
                'name': pool['name'],
                'initial_concs': list(pool['initial_concs']),
                'release': release,
                'inheritance': _validate_inheritance_block(
                    pool.get('inheritance', {})
                ),
                'fate': _validate_fate_block(
                    pool.get('fate', {}),
                    additive_name=additive_name,
                    pool_name=pool['name'],
                    valid_targets=valid_particulate_targets,
                    n_size_classes=n_size_classes,
                    allow_size_dependent=True
                )
            })

        medium_pools_out = []
        for pool in medium_pools:
            if not isinstance(pool, dict):
                raise SchemaError(f"Each medium_pool in additive '{additive_name}' must be a dict.")
            if 'name' not in pool:
                raise SchemaError(f"Each medium_pool in additive '{additive_name}' must contain 'name'.")
            initial_mass = float(pool.get('initial_mass', 0.0))
            if initial_mass < 0.0:
                raise SchemaError(f"medium_pool initial_mass must be non-negative for '{additive_name}:{pool['name']}'.")
            medium_pools_out.append({
                'name': pool['name'],
                'initial_mass': initial_mass,
                'fate': _validate_fate_block(
                    pool.get('fate', {}),
                    additive_name=additive_name,
                    pool_name=pool['name'],
                    valid_targets=valid_medium_targets,
                    n_size_classes=None,
                    allow_size_dependent=False
                )
            })

        out.append({
            'name': additive_name,
            'pools': pools_out,
            'medium_pools': medium_pools_out,
        })

    return out



def _normalize_transfer_target(raw_target: str, additive_name: str) -> str:
    """Return fully-qualified target name."""
    target = str(raw_target)
    return target if ':' in target else f"{additive_name}:{target}"


def _validate_fate_block(fate: dict,
                         additive_name: str,
                         pool_name: str,
                         valid_targets: set[str],
                         n_size_classes: int | None = None,
                         allow_size_dependent: bool = False) -> dict:
    if not isinstance(fate, dict):
        raise SchemaError(
            f"pool.fate for '{additive_name}:{pool_name}' must be a dict."
        )

    def _validate_rate(name, value):
        if allow_size_dependent:
            if not _is_nonnegative_scalar_or_array(value, n_expected=n_size_classes):
                raise SchemaError(
                    f"{name} must be a non-negative scalar or a length-{n_size_classes} "
                    f"non-negative array for '{additive_name}:{pool_name}'."
                )
        else:
            value = float(value)
            if value < 0.0:
                raise SchemaError(
                    f"{name} must be non-negative for '{additive_name}:{pool_name}'."
                )
        return value

    k_deg = _validate_rate('k_deg', fate.get('k_deg', 0.0))
    k_loss = _validate_rate('k_loss', fate.get('k_loss', 0.0))

    transfers_out = []
    for tr in fate.get('transfers', []):
        if not isinstance(tr, dict):
            raise SchemaError(
                f"Each transfer in '{additive_name}:{pool_name}' fate must be a dict."
            )
        if 'to' not in tr or 'k' not in tr:
            raise SchemaError(
                f"Each transfer in '{additive_name}:{pool_name}' fate must contain 'to' and 'k'."
            )

        k_val = tr['k']
        if allow_size_dependent:
            if not _is_nonnegative_scalar_or_array(k_val, n_expected=n_size_classes):
                raise SchemaError(
                    f"Transfer rate k must be a non-negative scalar or a length-{n_size_classes} "
                    f"non-negative array for '{additive_name}:{pool_name}'."
                )
        else:
            k_val = float(k_val)
            if k_val < 0.0:
                raise SchemaError(
                    f"Transfer rate k must be non-negative for '{additive_name}:{pool_name}'."
                )

        target = _normalize_transfer_target(tr['to'], additive_name)
        this_name = f"{additive_name}:{pool_name}"
        if target == this_name:
            raise SchemaError(f"Self-transfer is not allowed for '{this_name}'.")
        if target not in valid_targets:
            raise SchemaError(
                f"Unknown transfer target '{target}' for '{this_name}'."
            )

        transfers_out.append({'to': target, 'k': k_val})

    return {'k_deg': k_deg, 'k_loss': k_loss, 'transfers': transfers_out}

def _normalize_to_additives(data: dict, config: dict) -> dict:
    """
    Normalize old single-additive input to the new canonical
    data['additives'] structure.
    """
    data = dict(data)
    config = dict(config)

    if data.get('additives', None) is not None:
        return data

    init_A = data.get('initial_chemical_concs', data.get('initial_additive_concs', None))
    chem_cfg = config.get('chemical_release', config.get('additive_release', None))
    chem_data = data.get('chemical_release', data.get('additive_release', None))

    if init_A is None or chem_cfg is None or chem_data is None:
        return data

    data['additives'] = [{
        'name': 'additive_0',
        'pools': [{
            'name': 'pool_0',
            'initial_concs': list(init_A),
            'release': {
                'model': str(chem_cfg.get('model', 'analytical')).lower(),
                'solver': dict(chem_cfg.get('solver', {})),
                'params': dict(chem_data),
                'target': 'additive_0:medium',
            }
        }],
        'medium_pools': [{
            'name': 'medium',
            'initial_mass': 0.0,
            'fate': {},
        }]
    }]
    return data


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
    data = _normalize_to_additives(data, config)
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

    if validated.get('additives') is not None:
        validated['additives'] = _validate_additives_structure(
            validated['additives'],
            config['n_size_classes']
        )
    else:
        _validate_chemical_release_cross_checks(validated, config)

    # Backward-compatible aliases in validated output
    if 'initial_chemical_concs' in validated and 'initial_additive_concs' not in validated:
        validated['initial_additive_concs'] = validated['initial_chemical_concs']

    if 'chemical_release' in validated and 'additive_release' not in validated:
        validated['additive_release'] = validated['chemical_release']

    return validated