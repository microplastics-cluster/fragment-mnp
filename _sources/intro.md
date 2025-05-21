# FRAGMENT-MNP

FRAGMENT-MNP is a mechanistic model of Micro and NanoPlastic FRAGMentation in the ENvironmenT.

## Quickstart

The easiest way to install FRAGMENT-MNP is by using `pip`:

```bash
$ pip install fragmentmnp
```

FRAGMENT-MNP requires a minimum Python version of 3.9. Below is a bare minimum example. See also the [](example-usage.ipynb).

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

## Statement of need

By modifying particle sizes and shapes, fragmentation influences the potential risk caused by plastics in the environment. Therefore, to assess the risk caused by plastic pollution, understanding fragmentation is crucial. Predictive models are essential in helping this assessment, enabling us to fill gaps in observation data, better understand the results of experiments, and predict theoretical scenarios, such as in prospective risk assessments. Despite this, existing models that predict plastic transport, fate and exposure to organisms either do not consider fragmentation, include it only as a loss process {cite:p}`Quik:2023`, treat fragmentation as independent of the properties and residence time in the environment {cite:p}`Koelmans:2017`, or consider fragmentation as cascading (fragmenting mass can only be partitioned to the next biggest size class, rather than allowing the formation of nano-scale fragments) {cite:p}`Kaandorp:2021,Domercq:2022`. In reality, fragmentation depends on the environmental stresses encountered in the environment, such as photolysis by sunlight, hydrolysis by water, enzymatic action and mechanical disruption causing the break-apart of particles (e.g. the action of waves or bioturbation by soil invertebrates). Fragmentation has also been shown to often occur via surface erosion into nanoscale fragments {cite:p}`Meides:2021,Menzel:2022`, rather than in a cascading manner. There is a clear need for flexible and accessible model that can account for these factors, and FRAGMENT-MNP fills this gap.

## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO59.
