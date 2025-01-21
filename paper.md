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
    affilitation: 1
  - name: Cansu Uluseker
    orcid: 0000-0001-9828-0458
    affilitation: 1
  - name: Patrizia Pfohl
    orcid: 0000-0002-1580-8833
    affiliation: 2
  - name: Katherine Santizo
    affilitation:  "2, 3"
  - name: Joana Sipe
    orcid: 0000-0001-5602-0922
    affiliation: "4, 7"
  - name: Brandon Lopez
    affiliation: 4
  - name: Antonia Praetorius
    orcid: 0000-0003-0197-0116
    affiliation: 5
  - name: Richard Cross
    orcid: 0000-0001-5409-6552
    affiliation: 6
  - name: Gbotemi Adediran
    orcid: 0000-0001-6657-3336
    affiliation: 6
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
  - name: University of Amsterdam, 1012 WP Amsterdam, Netherlands
    index: 5
  - name: UK Centre for Ecology & Hydrology, Benson Lane, Crowmarsh Gifford, Wallingford, OX10 8BB, UK
    index: 6
  - name: Arizona State University, 1151 S Forest Ave, Tempe, AZ, United States
    index: 7
date: 21 January 2025
bibliography: paper.bib
---

# Summary

The degradation and fragmentation of plastics in the environment is an important but relatively poorly characterised process. Fragmentation leads to the formation of micro- and nanoplastics, and research has shown that particle size and shape, and thereby fragmentation, impacts a variety of processes, such as the ability of organisms to uptake plastics, and the movement of plastics around the environment [@Bucci:2022; @Thompson:2024]. In this paper, we present the FRAGMENT-MNP model as an open-source mechanistic model of micro- and nanoplastic degradation and fragmentation. FRAGMENT-MNP predicts the time evolution of particle size distributions, providing invaluable insights into fragmentation dynamics to help develop our understanding of plastic degradation and fragmentation, and offering predictive capabilities to enable better prediction of the fate and exposure of plastics in the environment.

# Statement of need

By modifying particle sizes and shapes, fragmentation influences the potential risk caused by plastics in the environment. Therefore, to assess the risk caused by plastic pollution, understanding fragmentation is crucial. Predictive models are essential in helping this assessment, enabling us to fill gaps in observation data, better understand the results of experiments, and predict theoretical scenarios, such as in prospective risk assessments. Despite this, existing models that predict plastic transport, fate and exposure to organisms either do not consider fragmentation, include it only as a loss process [@Quik:2023], treat fragmentation as independent of the properties and residence time in the environment [@Koelmans:2017], or consider fragmentation as cascading (fragmenting mass can only be partitioned to the next biggest size class, rather than allowing the formation of nano-scale fragments) [@Kaandorp:2021]. In reality, fragmentation depends on the environmental stresses encountered in the environment, such as photolysis by sunlight, hydrolysis by water, enzymatic action and mechanical disruption causing the break-apart of particles (e.g. the action of waves or bioturbation by soil invertebrates). Fragmentation has also been shown to often occur via surface erosion into nanoscale fragments [@Meides:2021; @Menzel:2022], rather than in a cascading manner. There is a clear need for flexible and accessible model that can account for these factors, and FRAGMENT-MNP fills this gap.

# Acknowledgements

Thank you to the European Chemical Industry Council Long-Range Research Initiative (Cefic-LRI) for providing funding for this work, under project number ECO59.

# References