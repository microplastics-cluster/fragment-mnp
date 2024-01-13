# FRAGMENT-MNP

🏗️ *Currently in development* 🏗️

FRAGMENT-MNP is a mechanistic model of Micro and NanoPlastic FRAGMentation in the ENvironmenT.

[See the full model documentation here.](https://microplastics-cluster.github.io/fragment-mnp)

## Usage

FRAGMENT-MNP is in active development at the moment - results should be viewed with caution! That being said, you are welcome to check out the package using the developer instructions below to install it. Also checkout the [example usage notebook](./docs/example-usage.ipynb).

```python
from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data
import matplotlib.pyplot as plt

# Create the model, pass it config and data, then run it
fmnp = FragmentMNP(minimal_config, minimal_data)
output = fmnp.run()

# Plot the results
plt.plot(output.t, output.n.T)
plt.legend([f'{d} m' for d in fmnp.psd])
plt.show()
```

## Developing FRAGMENT-MNP

We are using [Poetry](https://python-poetry.org/) for packaging and dependency management and this is the easiest way to get started. You can install Poetry following [these instructions](https://python-poetry.org/docs/#installation). This presumes you already have Python installed.

After you've installed Poetry, grab a copy of the code from this repo:

```shell
$ git clone https://github.com/microplastics-cluster/fragment-mnp.git
$ cd fragment-mnp
```

Then install the relevant dependencies using Poetry:

```shell
$ poetry install 
```

If you are not already in a virtual environment (i.e. you are not using Conda or virtualenv), this will create a virtual environment in which the package and its dependencies are installed. [These instructions](https://python-poetry.org/docs/basic-usage/#using-your-virtual-environment) detail how you can use this environment. If you are already in a virtual environment, the package and its dependencies will be installed into this.

You can test the installation has worked correctly by running `pytest` in the project root.

```shell
$ pytest
```

## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO59.
