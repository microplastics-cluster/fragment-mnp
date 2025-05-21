---
title: 'FRAGMENT-MNP: A model of micro- and nanoplastic fragmentation in the environment'
tags:
  - microplastics
  - plastics
  - environment
  - Python
authors:
  - name: Sam Harrison
    orcid: 0000-0001-8491-4720
    affiliation: 1
  - name: Cansu Uluseker
    orcid: 0000-0001-9828-0458
    affiliation: 1
  - name: Patrizia Pfohl
    orcid: 0000-0002-1580-8833
    affiliation: 2
  - name: Katherine Santizo
    affiliation: "2, 3"
  - name: Joana Sipe
    orcid: 0000-0001-5602-0922
    affiliation: "4, 7"
  - name: Brandon Lopez
    affiliation: 4
  - name: Antonia Praetorius
    orcid: 0000-0003-0197-0116
    affiliation: 5
  - name: Richard K. Cross
    orcid: 0000-0001-5409-6552
    affiliation: 6
  - name: Gbotemi A. Adediran
    orcid: 0000-0001-6657-3336
    affiliation: "6, 8"
  - name: Claus Svendsen
    orcid: 0000-0001-7281-647X
    affiliation: 6
  - name: Wendel Wohlleben
    orcid: 0000-0003-2094-3260
    affiliation: 2
  - name: Mark Wiesner
    orcid: 0000-0002-0823-5045
    affiliation: 4
affiliations:
  - name: UK Centre for Ecology & Hydrology, Library Avenue, Bailrigg, Lancaster, LA1 4AP, UK
    index: 1
  - name: BASF SE, Carl-Bosch-Str. 38, 67056 Ludwigshafen, Germany
    index: 2
  - name: Cefic aisbl, Rue Belliard 40 (box 15) – 1040, Brussels, Belgium
    index: 3
  - name: Duke University, Durham, NC 27708, United States
    index: 4
  - name: University of Amsterdam, 1090 GE Amsterdam, Netherlands
    index: 5
  - name: UK Centre for Ecology & Hydrology, Benson Lane, Crowmarsh Gifford, Wallingford, OX10 8BB, UK
    index: 6
  - name: Ira A. Fulton Schools of Engineering, Arizona State University, 1151 S Forest Ave, Tempe, AZ, United States
    index: 7
  - name: School of Earth and Environment, University of Leeds, Leeds, LS2 9JT, UK
    index: 8
date: 17 February 2025
bibliography: paper.bib
---

# Summary

The degradation and fragmentation of plastics in the environment is an important but relatively poorly characterised process. Fragmentation leads to the formation of micro- and nanoplastics, and research has shown that particle size and shape, and thereby fragmentation, impacts a variety of processes, such as the ability of organisms to uptake plastics, and the movement of plastics around the environment [@ThorntonHampton:2022; @Thompson:2024]. In this paper, we present the FRAGMENT-MNP model as an open-source mechanistic model of micro- and nanoplastic degradation and fragmentation. FRAGMENT-MNP predicts how pieces of plastic can break part (fragment), providing invaluable insights into fragmentation dynamics to help develop our understanding of plastic degradation and fragmentation, and offering predictive capabilities to enable better prediction of the fate and exposure of plastics in the environment.

# Statement of need

By modifying particle sizes and shapes, fragmentation influences the potential risk caused by plastics in the environment. Therefore, to assess the risk caused by plastic pollution, understanding fragmentation is crucial. Predictive models are essential in helping this assessment, enabling us to fill gaps in observation data, better understand the results of experiments, and predict theoretical scenarios, such as in prospective risk assessments. Despite this, existing models that predict plastic transport, fate and exposure to organisms either do not consider fragmentation, include it only as a loss process [@Quik:2023], treat fragmentation as independent of the properties and residence time in the environment [@Koelmans:2017], or consider fragmentation as cascading (fragmenting mass can only be partitioned to the next biggest size class, rather than allowing the formation of nano-scale fragments) [@Kaandorp:2021; @Domercq:2022]. In reality, fragmentation depends on the environmental stresses encountered in the environment, such as photolysis by sunlight, hydrolysis by water, enzymatic action and mechanical disruption causing the break-apart of particles (e.g. the action of waves or bioturbation by soil invertebrates). Fragmentation has also been shown to often occur via surface erosion into nanoscale fragments [@Meides:2021; @Menzel:2022], rather than in a cascading manner. There is a clear need for flexible and accessible model that can account for these factors, and FRAGMENT-MNP fills this gap.

# Overview of the model

The model is fully documented at [https://fragmentmnp.ceh.ac.uk](https://fragmentmnp.ceh.ac.uk). Here, we provide a brief overview of its main conceptualisation and functionality.

Particle concentrations $c_k$ are represented in binned size classes $k$, with the model allowing for the fragmentation of particles from larger to smaller size classes, and the dissolution of particles into a dissolved size class with concentration $c_\text{diss}$. Conceptually, the dissolved fraction consists of oligomers, monomers and volatile organic compounds. The solutions are obtained by numerically solving the following set of differential equations:

$$
\frac{dc_k}{dt} = -k_{\text{frag},k} c_k + \sum_i f_{i,k} k_{\text{frag},i} c_i - k_{\text{diss},k} c_k
$$

$$
\frac{dc_\text{diss}}{dt} = \sum_k k_{\text{diss},k} c_k
$$

Here, $k_{\text{frag},k}$ is the fragmentation rate constant of size class $k$, $f_{i,k}$ is the fraction of daughter fragments produced from a fragmenting particle of size $i$ that are of size $k$, and $k_{\text{diss},k}$ is the dissolution rate from size class $k$. The rate constants $k_\text{frag}$ and $k_\text{diss}$ can be a function of time and particle surface area, and thus are represented internally as 2D arrays. The shape of these dependencies can be controlled by model input parameters, which allow for constant, linear, polynomial, power law, exponential, logarithmic or logistic dependencies. This flexible parameterisation allows for rate constants to model a variety of physical phenomena, such as the modulation of fragmentation as polymer particles undergo weathering in the environment.

FRAGMENT-MNP is released as a `pip` package (`pip install fragmentmnp`) and comes with example configuration and data,to enable users to get started using the model in as little time as possible. A bare minimum model run is shown below, which uses these examples to run an arbitrary fragmentation scenario (with no dissolution):

```python
from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data

# Create the model and pass it the example config and data
fmnp = FragmentMNP(minimal_config, minimal_data)
# Run the model
output = fmnp.run()
# Plot the results
output.plot()
```

![Example model output showing the time evolution of particle size distributions undergoing fragmentation](./paper/fragmentmnp_example.png){height="900px"}

# Related work

The development of FRAGMENT-MNP is part of broader efforts to further our understanding of microplastics in the environment. This includes the development of other models covering emissions, additive release, long-range transport and exposure. Ongoing work is pursuing closer integration of these models and associated experimental data in order to provide holistic predictions of microplastic fate and exposure in the environment.


# Acknowledgements

Thank you to the European Chemical Industry Council Long-Range Research Initiative (Cefic-LRI) for providing funding for this work, under project number ECO59.

# References
