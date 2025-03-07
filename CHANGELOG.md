# Changelog
Notable changes to FRAGMENT-MNP will be documented here. We are using [semantic versioning](https://semver.org/), though note pre-release versions (0.0.x) may include breaking changes.


## [Unreleased]

## [1.1.4] - 2025-03-07

* No changes to code, only moving Poetry dependencies to test environment.

## [1.1.3] - 2025-02-04

* No changes to code, only adding link to GitHub repo issues to docs.

## [1.1.2] - 2025-02-04

* No changes to code, only updating docs to include UKCEH logo and privacy notice.

## [1.1.1] - 2025-01-24

* No changes to code, only correcting model config and input data docs to include `k_min`.


## [1.1.0] - 2025-01-23

* Added mineralisation as a process that transfers mass from the dissolved pool to a mineralised polymer pool. This is controlled by a separate `k_min` rate constant.
* `c_diss_from_sc` and `n_diss_from_sc` are no longer stored in output data due to a change in the ODEs solved to incorporate mineralisation.
* Updated dissolution yaxis label from "Dissolution" to "Dissolved"


## [1.0.1] - 2025-01-13

* Improvements and bug fixes in documentation.


## [1.0.0] - 2025-01-10

* The creation of rate constant distributions (across time and size class) has been altered significantly, resulting in breaking changes to model config and data. In an attempt to unify the treatment of all rate constants, `k_frag` and `k_diss` are now both processed by the same function (FragmentMNP.set_k_distribution), which enables the creation of time and size class distributions that are built from a combination of power law, exponential, logarithmic and logistic regressions. The parameters that control these distributions are named differently to the power law relationships in previous versions (e.g. `theta_1`, `k_frag_tau`) and these parameters can no longer be present in input data. The new distributions give much more flexibility in controlling the time and size dependence of rate constants. The default is for rate constants to be constant in time and particle size. See the API docs for full details.
* The default solver has changed to LSODA. Fragmentation problems quite often end up being stiff, so having the automatic stiffness detection provided by LSODA is a sensible default.
* Previously, the model used timestep indices as the time grid provided to the numerical solver (i.e. `np.arange(0, n_timesteps)`), thus ignoring the config value for timestep length `dt`. As of this version, `dt` is taken into account and the real time grid is provided to the numerical solver. For very large `dt` values, you may wish to scale yourself. Additionally, the model grid uses the midpoint of the timesteps in order to prevent the first time point always being equal to 0 (which has impacts on rate constant regressions).
* The 2D `k_frag` array has been transposed such that the dimensions are now `(n_size_classes, n_timesteps)`, rather than `(n_timesteps, n_size_classes)`. This is to align the array with the model output variables.
* Initial dissolved concentrations can now be specified using the `initial_concs_diss` data parameter. The `FMNPOutput` data attributes have been changed accordingly: `c_diss` now represents the *total* dissolved fraction timeseries, summed as size classes and including the initial dissolved concentrations. `c_diss_from_sc` is concentration lost from each size class due to dissolution, but does not included the initial dissolved concentration. `n_diss_from_sc` is the same as `c_diss_from_sc` but converted to particle number concentrations.

## [0.1.0] - 2024-01-13

* Initial model version.


[Unreleased]: https://github.com/microplastics-cluster/fragment-mnp/compare/1.1.4...HEAD
[1.1.4]: https://github.com/microplastics-cluster/fragment-mnp/releases/tag/1.1.4
[1.1.3]: https://github.com/microplastics-cluster/fragment-mnp/releases/tag/1.1.3
[1.1.2]: https://github.com/microplastics-cluster/fragment-mnp/releases/tag/1.1.2
[1.1.1]: https://github.com/microplastics-cluster/fragment-mnp/releases/tag/1.1.1
[1.1.0]: https://github.com/microplastics-cluster/fragment-mnp/releases/tag/1.1.0
[1.0.1]: https://github.com/microplastics-cluster/fragment-mnp/releases/tag/1.0.1
[1.0.0]: https://github.com/microplastics-cluster/fragment-mnp/releases/tag/1.0.0
[0.1.0]: https://github.com/microplastics-cluster/fragment-mnp/releases/tag/0.1.0