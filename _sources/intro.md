# FRAGMENT-MNP

:::{caution}
🏗️ *Currently in development* 🏗️
:::

FRAGMENT-MNP is a mechanistic model of Micro and NanoPlastic FRAGMentation in the ENvironmenT.

## Quickstart

FRAGMENT-MNP is in active development at the moment - results should be viewed with caution! That being said, you are welcome to check out the package using the developer instructions below to install it. Also checkout the [](example-usage.ipynb).

```python
from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, data
import matplotlib.pyplot as plt

# Create the model, pass it config and data, then run it
fmnp = FragmentMNP(minimal_config, data)
output = fmnp.run()

# Plot the results
plt.plot(output.t, output.n.T)
plt.legend([f'{d} m' for d in fmnp.psd])
plt.show()
```

## Developing FRAGMENT-MNP

Check out the [](developers/quickstart.md) to get started developing the model.

## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO59.