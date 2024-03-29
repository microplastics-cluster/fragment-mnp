# Model configuration

When initialising the {class}`fragmentmnp.FragmentMNP` model class, a config dict must be provided. Examples are given in the {mod}`fragmentmnp.examples` module:

```python
>>> from fragmentmnp.examples import minimal_config, full_config
>>> minimal_config
{'n_size_classes': 7, 'particle_size_range': [-9, -3], 'n_timesteps': 100}
>>> full_config
{'n_size_classes': 7, 'particle_size_range': [-9, -3], 'n_timesteps': 100, 'dt': 1,
'k_diss_scaling_method': 'constant'}
```

The `minimal_config` contains only required variables, whilst `full_config` includes variables that have defaults. Here we take a look at the schema for this config dict.

`n_size_classes`
: *Required, int, less than or equal to 100.*
: The number of bins the model should use for particle size distributions.

`particle_size_classes`
: *Optional, list of int or floats of length `n_size_classes`, units: m.*
: The mean diameter of each particle size class. If `particle_size_classes` is not specified, then `particle_size_range` will be used to calculate this diameters instead. Note that whilst `particle_size_range` requires exponents of 10 to be provided, `particle_size_classes` requires absolute values.

`particle_size_range`
: *Optional, list/tuple of int or floats of length 2, units: m.*
: The particle size range that should be included in the model in the format [min, max]. The values given are used as the exponent of 10, so, for example, `[-9, -3]` results in a size range of 10e-9 m to 10e-3 m (1 nm to 1 mm). The size class bins themselves will be deduced by evenly splitting the range by `n_size_classes` on a log scale. If neither `particle_size_range` nor `particle_size_classes` is specified, you will get an error.

`n_timesteps`
: *Required, int.*
: The number of timesteps to run the model for. Each timestep is of length `dt`, which can be provided as a parameter and defaults to 1 second.

`dt`
: *Optional, int, units: s, default: 1*.
: The length of each timestep. Defaults to 1 second.

`k_diss_scaling_method`
: *Optional, str equal to "constant" or "surface_area", default: "constant".*
: The method by which the dissolution rate `k_diss` is scaled across size classes. If `constant`, then `k_diss` is the same across all size classes. If `surface_area`, then `k_diss` scales as $s^\gamma$, where $s$ is the surface area to volume ratio of the polymer particles, and $\gamma$ is an empirical scaling parameter set by the `k_diss_gamma` variable in the [model input data dict](input-data). If `k_diss` is given as a distribution in the input data, then `k_diss_scaling_method` is ignored and this distribution is used directly instead.

(config:solver_method)=
`solver_method`
: *Optional, str equal to a valid method, default: "RK45".*
: The method that is used to numerically solve the model differential equation. This method is passed directly to [`scipy.intergrate.solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) and must be a valid method available in SciPy ("RK45", "RK23", "DOP853", "Radau", "BDF" or "LSODA"). By default, RK45 is used. See the [SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) for more details. If you are having numerical stability issues, try "Radau", "BDF" or "LSODA", and read [](./advanced-usage/numerical-instability.ipynb).

(config:solver_rtol_atol)=
`solver_rtol`, `solver_atol`
: *Optional, float or list of floats for atol, defaults: rtol=1e-3, atol=1e-6*
: Relative and absolute tolerances passed to [scipy.integrate.solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html). The solver keeps the local error estimates less than `solver_atol + solver_rtol * abs(y)`. For `solver_atol`, a list of length `n_size_classes` can be given, such that a different absolute tolerance is used for each size class. Setting these tolerances might be particularly useful for stiff problems that cause numerical instability.

(config:solver_max_step)=
`solver_max_step`
: *Optional, float, default: `np.inf`.*
: Set a maximum step size for the ODE solver. This might be particularly useful for stiff problems that cause numerical instability.

(config:solver_t_eval)=
`solver_t_eval`
: *Optional, str equal to "integer", None, or list of integers or floats, default: "integer".*
: Time points at which to store the computed solution, passed to [scipy.integrate.solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html). If the string "integer" is provided, the model stores iterative integer timesteps with a step of 1 from 0 to `n_timesteps`, i.e. `solver_t_eval = np.arange(n_timesteps)`. If `None`, the points selected by the solver are used (which might be controlled by the `solver_max_step` option). If a list of integers or floats is provided, these are used directly as the time points. An error will be thrown if these lie outside of the model time span.