# Changelog
Notable changes to FRAGMENT-MNP will be documented here. We are using [semantic versioning](https://semver.org/), though note pre-release versions (0.0.x) may include breaking changes.


## [Unreleased]

* The creation of rate constant distributions (across time and size class) has been altered significantly, resulting in breaking changes to model config and data. In an attempt to unify the treatment of all rate constants, `k_frag` and `k_diss` are now both processed by the same function (FragmentMNP.set_k_distribution), which enables the creation of time and size class distributions that are built from a combination of power law, exponential, logarithmic and logistic regressions. The parameters that control these distributions are named differently to the power law relationships in previous versions (e.g. `theta_1`, `k_frag_tau`) and these parameters can no longer be present in input data. The new distributions give much more flexibility in controlling the time and size dependence of rate constants. The default is for rate constants to be constant in time and particle size. See the API docs for full details.
* The default solver has changed to LSODA. Fragmentation problems quite often end up being stiff, so having the automatic stiffness detection provided by LSODA is a sensible default.
* Previously, the model used timestep indices as the time grid provided to the numerical solver (i.e. `np.arange(0, n_timesteps)`), thus ignoring the config value for timestep length `dt`. As of this version, `dt` is taken into account and the real time grid is provided to the numerical solver. For very large `dt` values, you may wish to scale yourself.

## [0.1.0] - 2024-01-13

* Initial model version.


[Unreleased]: https://github.com/microplastics-cluster/fragment-mnp/tree/develop 
[0.1.0]: https://github.com/microplastics-cluster/fragment-mnp/releases/tag/0.1.0