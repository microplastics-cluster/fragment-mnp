# FRAGMENT-MNP

FRAGMENT-MNP is a mechanistic model of Micro and NanoPlastic FRAGMentation in the ENvironmenT.

[See the full model documentation here.](https://microplastics-cluster.github.io/fragment-mnp)

## Usage

Install the model with `pip`:

```bash
$ pip install fragmentmnp
```

FRAGMENT-MNP requires a minimum Python version of 3.9. Then run the model with example data and plot the results:

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

## Statement of need

Fragmentation influences the potential risk caused by plastics in the environment by modifying its size and shape. Therefore, to assess the risk caused by plastic pollution, understanding fragmentation is crucial. Predictive models are essential in helping this assessment, enabling us to fill gaps in observation data, better understand the results of experiments, and predict theoretical scenarios, such as in prospective risk assessments. Despite this, existing models that predict plastic transport, fate and exposure to organisms either do not consider fragmentation, include it only as a loss process, treat fragmentation as independent of the properties and residence time in the environment, or consider fragmentation as cascading (fragmenting mass can only be partitioned to the next biggest size class, rather than allowing the formation of nano-scale fragments). In reality, fragmentation depends on the environmental stresses encountered in the environment, such as photolysis by sunlight, hydrolysis by water, enzymatic action and mechanical disruption causing the break-apart of particles (e.g. the action of waves or bioturbation by soil invertebrates). Fragmentation has also been shown to often occur via surface erosion into nanoscale fragments, rather than in a cascading manner. There is a clear need for flexible and accessible model that can account for these factors, and FRAGMENT-MNP fills this gap.

## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO59.
