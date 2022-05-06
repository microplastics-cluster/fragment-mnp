# Input data

Input data must be passed to the {class}`fragmentmnp.FragmentMNP` model class. Example are given in the {mod}`fragmentmnp.examples` module:

```python
>>> from fragmentmnp.examples import minimal_data, full_data
>>> minimal_data
{'initial_concs': [42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0], 'k_frag': 0.01}
>>> full_data
{'initial_concs': [42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0], 'k_frag': 0.01, 'theta_1': 0.0}
```

The `minimal_data` contains only required variables, whilst `full_data` includes variables that have defaults. Here we take a look at the schema for this data dict.

`initial_concs`
: *Required, list of floats with length equal to `n_size_classes`, units: particles/m3.*
: The intial particle number concentration in each model size class. 

`k_frag`
: *Required, float, units: s<sup>-1</sup>.*
: The average fragmentation rate $k_\text{frag}$ across the size classes.

`theta_1`
: *Optional, float, default: 0.*
: Surface energy empirical parameter $\theta_1$, which dictates how $k_\text{frag}$ varies with particle size. $k_\text{frag}$ varies as $d_{2\theta_1}$, meaning if $\theta_1 = 0$, the same $k_\text{frag}$ is used across size classes, and if $\theta_1 > 0$, $k_\text{frag}$ is larger for larger size particles.