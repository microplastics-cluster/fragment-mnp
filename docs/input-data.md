# Input data

Input data must be passed to the {class}`fragmentmnp.FragmentMNP` model class. Examples are given in the {mod}`fragmentmnp.examples` module:

```python
>>> from fragmentmnp.examples import minimal_data, full_data
>>> minimal_data
{'initial_concs': [42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0], 'density': 1380, 'k_frag': 0.01}
>>> full_data
{'initial_concs': [42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0], 'density': 1380, 'k_frag': 0.01,
'theta_1': 0.0, 'k_diss': 0.0, 'k_diss_gamma': 1.0}
```

The `minimal_data` contains only required variables, whilst `full_data` includes variables that have defaults. Here we take a look at the schema for this data dict.

`initial_concs`
: *Required, list of positive floats with length equal to `n_size_classes`, units: kg/m3.*
: The intial mass concentration in each model size class. 

`density`
: *Required, float, units: kg/m3.*
: The density of the polymer being modelled.

`k_frag`
: *Required, float or iterable with length equal to `n_size_classes`, units: s<sup>-1</sup>.*
: Either a scalar representing the average fragmentation rate $k_\text{frag}$ across the size classes, or the full $k_\text{frag}$ distribution for all size classes. $k_\text{frag}$ is defined as the fraction of the *mass* of a particular size class that fragments each second. If a scalar (average) is provided, `theta_1` is used to calculate the distribution, otherwise any value given for `theta_1` is ignored.

`theta_1`
: *Optional, float, default: 0.*
: Surface energy empirical parameter $\theta_1$, which dictates how $k_\text{frag}$ varies with particle size. $k_\text{frag}$ varies as $d_{2\theta_1}$, meaning if $\theta_1 = 0$, the same $k_\text{frag}$ is used across size classes, and if $\theta_1 > 0$, $k_\text{frag}$ is larger for larger size particles. `theta_1` is ignored if a distribution is given for `k_frag`.

`k_diss`
: *Optional, float, default: 0.*
: Either the average (median) dissolution rate $k_\text{diss}$ across the size classes, or the full $k_\text{diss}$ distribution for all size classes. $k_\text{diss}$ describes the rate of loss to "dissolved" organics from each size class. If an average (scalar) is given, how this average is distributed across size classes is determined by the `k_diss_scaling_method` in the [model configuration](config). If a distribution is given, the `k_diss_scaling_method` is ignored.

`k_diss_gamma`
: *Optional, float, default: 1.*
: If the config option `k_diss_scaling_method` is set to `surface_area` and an average (scalar) `k_diss` is given in input data, then `k_diss_gamma` $\gamma$ is an empirical parameter that scales how dependent the distribution of `k_diss` across size classes is on the particle surface area to volume ratio. In short, $k_\text{diss} \propto s^\gamma$, where $s$ is the surface area to volume ratio. In other words, if $\gamma = 1$, `k_diss` scales directly with the surface area to volume ratio, and if $\gamma = 0$, `k_diss` is constant across size classes (i.e. it has the same effect as setting `k_diss_scaling_method` to `constant`). If a distribution is given for `k_diss`, then `k_diss_gamma` is ignored.


```{admonition} Units
:class: tip
Strictly speaking, the FRAGMENT-MNP model is unit-agnostic. As long as the units of data specified in input data are consistent, then the units of variables output by the model will match these. For example, if you specify `initial_concs` and `density` with units of mg/m<sup>3</sup>, then output concentrations will be in mg/m<sup>3</sup>.

**However**, we strongly recommend you stick to SI units to avoid the chance of mistakes being made. Hence this page gives unit recommendations in SI units.
```