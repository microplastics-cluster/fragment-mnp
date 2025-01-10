# FRAGMENT-MNP

FRAGMENT-MNP is a mechanistic model of Micro and NanoPlastic FRAGMentation in the ENvironmenT.

[See the full model documentation here.](https://microplastics-cluster.github.io/fragment-mnp)

## Usage

Install the model with `pip`:

```bash
$ pip install fragmentmnp
```

Then run the model with example data and plot the results:

```python
from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data

# Create the model, pass it example config and data, then run it
fmnp = FragmentMNP(minimal_config, minimal_data)
output = fmnp.run()
# Plot the mass concentration time series
output.plot()
```

## Developing FRAGMENT-MNP

[See the developers documentation.](https://microplastics-cluster.github.io/fragment-mnp/developers/quickstart.html)

## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO59.
