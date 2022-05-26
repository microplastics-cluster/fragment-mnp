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
: *Required, list of floats with length equal to `n_size_classes`, units: particles/m3.*
: The intial particle number concentration in each model size class. 

`density`
: *Required, float, units: kg/m3.*
: The density of the polymer being modelled.

`k_frag`
: *Required, float, units: s<sup>-1</sup>.*
: The average fragmentation rate $k_\text{frag}$ across the size classes.

`theta_1`
: *Optional, float, default: 0.*
: Surface energy empirical parameter $\theta_1$, which dictates how $k_\text{frag}$ varies with particle size. $k_\text{frag}$ varies as $d_{2\theta_1}$, meaning if $\theta_1 = 0$, the same $k_\text{frag}$ is used across size classes, and if $\theta_1 > 0$, $k_\text{frag}$ is larger for larger size particles.

`k_diss`
: *Optional, float, default: 0.*
: The average (median) dissolution rate $k_\text{diss}$, which describes the rate of loss to dissolved organics from each size class. How this average is distributed across size classes is determined by the `k_diss_scaling_method` in the [model configuration](config).

`k_diss_gamma`
: *Optional, float, default: 1.*
: If the config option `k_diss_scaling_method` is set to `surface_area`, then `k_diss_gamma` $\gamma$ is an empirical parameter that scales how dependent the distribution of `k_diss` across size classes is on the particle surface area to volume ratio. In short, $k_\text{diss} \propto s^\gamma$, where $s$ is the surface area to volume ratio. In other words, if $\gamma = 1$, `k_diss` scales directly with the surface area to volume ratio, and if $\gamma = 0$, `k_diss` is constant across size classes (i.e. it has the same effect as setting `k_diss_scaling_method` to `constant`).
