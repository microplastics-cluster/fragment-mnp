# Configuration

When initialising the {class}`fragmentmnp.FragmentMNP` model class, a config dict must be provided. Examples are given in the {mod}`fragmentmnp.examples` module:

```python
>>> from fragmentmnp.examples import minimal_config, full_config
>>> minimal_config
{'n_size_classes': 7, 'particle_size_range': [-9, -3], 'n_timesteps': 100}
>>> full_config
{'n_size_classes': 7, 'particle_size_range': [-9, -3], 'n_timesteps': 100,
'dt': 1, 'allow_loss': False}
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

`allow_loss`
: *Optional, bool, default: False*.
: Whether loss (fragmentation) from the smallest size class is allowed.