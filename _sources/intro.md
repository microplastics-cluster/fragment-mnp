# FRAGMENT-MNP

FRAGMENT-MNP is a mechanistic model of Micro and NanoPlastic FRAGMentation in the ENvironmenT.

## Quickstart

The easiest way to install FRAGMENT-MNP is by using `pip`:

```bash
$ pip install fragmentmnp
```

Here is a bare minimum example. See also the [](example-usage.ipynb).

```python
from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data

# Create the model, pass it example config and data, then run it
fmnp = FragmentMNP(minimal_config, minimal_data)
output = fmnp.run()
# Plot the mass concentration time series
output.plot()
```

:::{note}
If you are running the model via a script, rather than an interactive environment (e.g. Jupyter), you can pass `show=True` to the `plot` function in order to automatically display the plot, i.e. `output.plot(show=True)`.
:::

## Issues

Issues, comments or questions? Head over to our [GitHub repository](https://github.com/microplastics-cluster/fragment-mnp) and [post an issue](https://github.com/microplastics-cluster/fragment-mnp/issues).

## Developing FRAGMENT-MNP

Check out the [](developers/quickstart.md) to get started developing the model.

## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO59.
