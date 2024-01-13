# FRAGMENT-MNP

:::{caution}
🏗️ *Currently in development* 🏗️
:::

FRAGMENT-MNP is a mechanistic model of Micro and NanoPlastic FRAGMentation in the ENvironmenT.

## Quickstart

FRAGMENT-MNP is in active development at the moment - results should be viewed with caution! That being said, you can install from GitHub using `pip`:

```bash
$ pip install git+https://github.com/microplastics-cluster/fragment-mnp.git
```

If you wish to contribute to development, then check out the [developer instructions](developers/quickstart) for instructions on how to install using Poetry.

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

## Developing FRAGMENT-MNP

Check out the [](developers/quickstart.md) to get started developing the model.

## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO59.
