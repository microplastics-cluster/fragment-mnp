"""
Main FRAGMENT-MNP model class (:mod:`fragmentmnp.fragmentmnp`)
=============================================================

This file contains the FragmentMNP model class.

COLLABORATION UPDATE
---------------------------------------------

We add optional additive tracking and analytical release modelling,
WITHOUT changing the core polymer fragmentation ODE system.

Key design choice:
- The polymer model is solved using solve_ivp exactly as before.
- After solving, if additive inputs are provided, we perform a simple
  explicit bookkeeping loop over the solver output time grid:
    (1) transport additive with fragmentation (mass-conserving)
    (2) release additive to water using an analytical model function

This "operator splitting" keeps the fragmentation core stable and makes
the collaboration easy: the analytical and/or numerical formula is isolated in one method:
    _analytical_additive_release_fraction()
    _numerical_additive_release_fraction()

"""
from typing import Tuple, Sequence
from functools import lru_cache
import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp
from scipy import interpolate
from scipy.optimize import brentq
from scipy.special import j0, j1, jn_zeros
from schema import SchemaError
from . import validation
from .output import FMNPOutput
from ._errors import FMNPNumericalError, FMNPDistributionValueError
from .geometry import make_geometry


class FragmentMNP():
    """
    The class that controls usage of the FRAGMENT-MNP model

    Parameters
    ----------
    config : dict
        Model configuration options
    data : dict
        Model input data
    validate: bool, default=True
        Should config and data be validated? It is strongly recommended to
        use validation, but this option is provided if you are *certain* your
        config and data are correct and wish to speed up model initialisation
    """

    def __init__(self,
                 config: dict,
                 data: dict,
                 validate: bool = True) -> None:
        """
        Initialise the model.

        Component/layer update
        ----------------------
        The original FRAGMENT-MNP state vector was size-resolved only.  This
        version keeps that public behaviour, but internally represents the
        particulate polymer state as a component-by-size matrix:

            c_component[j, k]

        where ``j`` is the formulation/layer/component index and ``k`` is the
        particle-size class.  When ``data['components']`` is not provided, the
        model creates one implicit component from the legacy top-level inputs,
        so existing scripts remain valid.
        """
        # Validate the config and data (if we are meant to)
        if validate:
            config, data = self._validate_inputs(config, data)

        self.config = config
        self.data = data

        # Set the number of particle size classes and timesteps, and the times
        # at which to store the computed solution.
        self.n_size_classes = self.config['n_size_classes']
        self.n_timesteps = self.config['n_timesteps']
        self.dt = self.config['dt']
        self.t_grid = np.arange(0.5 * self.dt,
                                self.n_timesteps * self.dt + 0.5 * self.dt,
                                self.dt)
        self.t_eval = self.t_grid \
            if self.config['solver_t_eval'] == 'timesteps' \
            else self.config['solver_t_eval']

        # Particle geometry shared by all components.  The legacy/default
        # geometry is a sphere, where psd values are particle diameters.  When
        # shape='fibre', the same size coordinate is interpreted as fibre
        # length and a separate fibre diameter defines the cross-section.
        self.psd = self._set_psd()
        self.geometry = make_geometry(
            self.config.get('particle_geometry', {'shape': 'sphere'}),
            n_size_classes=self.n_size_classes,
        )
        self.particle_geometry = self.geometry.metadata()
        self.surface_areas = self.geometry.surface_area(self.psd)
        self.particle_volumes = self.geometry.volume(self.psd)
        self.release_radii = self.geometry.release_radius(self.psd)
        self.fsd = self.set_fsd(self.n_size_classes,
                                self.psd,
                                self.data['fsd_beta'])

        # ------------------------------------------------------------------
        # NEW: component / formulation / layer dimension for the polymer core
        # ------------------------------------------------------------------
        component_defs = self.data.get('components', None)
        if component_defs is None:
            component_defs = [{
                'name': 'polymer_0',
                'initial_concs': list(self.data['initial_concs']),
                'initial_concs_diss': float(self.data.get('initial_concs_diss', 0.0)),
                'density': float(self.data['density']),
                'k_frag': self.data['k_frag'],
                'k_diss': self.data.get('k_diss', 0.0),
                'k_min': self.data.get('k_min', 0.0),
                'fsd_beta': self.data.get('fsd_beta', 0.0),
                'layer_thickness': None,
            }]

        self.components = component_defs
        self.n_components = len(component_defs)
        self.component_names = []
        self.component_layer_initial_thickness = np.full(self.n_components, np.nan, dtype=float)

        component_initial_concs = []
        component_initial_concs_diss = []
        component_density = []
        component_fsd = []
        component_k_frag = []
        component_k_diss = []
        component_k_min = []

        for j, component in enumerate(component_defs):
            self.component_names.append(str(component.get('name', f'component_{j}')))
            component_initial_concs.append(
                np.asarray(component['initial_concs'], dtype=float)
            )
            component_initial_concs_diss.append(
                float(component.get('initial_concs_diss', 0.0))
            )
            component_density.append(float(component.get('density', self.data['density'])))
            component_fsd.append(
                self.set_fsd(self.n_size_classes,
                             self.psd,
                             float(component.get('fsd_beta', self.data.get('fsd_beta', 0.0))))
            )
            component_k_frag.append(
                self.set_rate_constant(component.get('k_frag', self.data['k_frag']), 'k_frag')
            )
            component_k_diss.append(
                self.set_rate_constant(component.get('k_diss', self.data.get('k_diss', 0.0)), 'k_diss')
            )
            component_k_min.append(
                self.set_rate_constant(component.get('k_min', self.data.get('k_min', 0.0)), 'k_min')
            )

            layer_h = component.get('layer_thickness', None)
            if layer_h is not None:
                self.component_layer_initial_thickness[j] = float(layer_h)

        self.component_initial_concs = np.asarray(component_initial_concs, dtype=float)
        self.component_initial_concs_diss = np.asarray(component_initial_concs_diss, dtype=float)
        self.component_density = np.asarray(component_density, dtype=float)
        self.component_fsd = np.asarray(component_fsd, dtype=float)          # (C, N, N)
        self.component_k_frag = np.asarray(component_k_frag, dtype=float)    # (C, N, T)
        self.component_k_diss = np.asarray(component_k_diss, dtype=float)    # (C, N, T)
        self.component_k_min = np.asarray(component_k_min, dtype=float)      # (C, T)

        # Backward-compatible aggregate aliases used by the old API/tests.
        self.initial_concs = np.sum(self.component_initial_concs, axis=0)
        self.initial_concs_diss = float(np.sum(self.component_initial_concs_diss))
        total_initial_mass = float(np.sum(self.component_initial_concs))
        if total_initial_mass > 0.0:
            weights = np.sum(self.component_initial_concs, axis=1) / total_initial_mass
            self.density = float(np.sum(weights * self.component_density))
        else:
            self.density = float(np.mean(self.component_density))

        # Backward-compatible rate aliases.  For multi-component systems these
        # are mass-weighted effective rates used only by legacy helper logic and
        # old additive post-processing.  The core ODE uses component-specific
        # rates stored above.
        self.k_frag = self._component_weighted_rate(self.component_k_frag)
        self.k_diss = self._component_weighted_rate(self.component_k_diss)
        self.k_min = np.average(self.component_k_min, axis=0, weights=self._component_mass_weights())

        # ------------------------------------------------------------
        # OPTIONAL: multi-additive / multi-pool coupling inputs
        # ------------------------------------------------------------
        self.additives = data.get('additives', None)

        # Backward-compatible aliases for legacy single-additive inputs.
        init_A = data.get('initial_chemical_concs', data.get('initial_additive_concs', None))
        self.initial_chemical_concs = None if init_A is None else np.array(init_A, dtype=float)
        self.chemical_release_config = config.get(
            'chemical_release',
            config.get('additive_release', None)
        )
        self.chemical_release_data = data.get(
            'chemical_release',
            data.get('additive_release', None)
        )

        # Flatten particulate pools and medium pools into internal species.
        # Particulate species are size-resolved and move with fragmentation.
        # Medium species are bulk pools updated in post-processing only.
        self.chemical_species = []
        self.additive_names = []
        self.species_names = []
        self.species_additive_index = []
        self.medium_species = []
        self.medium_species_names = []
        self.medium_species_additive_index = []

        if self.additives is not None:
            for a_idx, additive in enumerate(self.additives):
                additive_name = additive['name']
                self.additive_names.append(additive_name)

                particulate_pools = list(additive['pools'])
                medium_pools = list(additive.get('medium_pools', []))
                if len(medium_pools) == 0:
                    medium_pools = [{
                        'name': 'medium',
                        'initial_mass': 0.0,
                        'fate': {'k_deg': 0.0, 'k_loss': 0.0, 'transfers': []},
                    }]

                for pool in particulate_pools:
                    self.chemical_species.append({
                        'additive_name': additive_name,
                        'pool_name': pool['name'],
                        'initial_concs': np.array(pool['initial_concs'], dtype=float),
                        'model': str(pool['release']['model']).lower(),
                        'solver': dict(pool['release'].get('solver', {})),
                        'params': dict(pool['release'].get('params', {})),
                        'release_target': str(pool['release'].get('target', f"{additive_name}:medium")),
                        'inheritance': dict(pool.get('inheritance', {'mode': 'proportional'})),
                        'fate': dict(pool.get('fate', {})),
                    })
                    self.species_names.append(f"{additive_name}:{pool['name']}")
                    self.species_additive_index.append(a_idx)

                for mpool in medium_pools:
                    self.medium_species.append({
                        'additive_name': additive_name,
                        'pool_name': mpool['name'],
                        'initial_mass': float(mpool.get('initial_mass', 0.0)),
                        'fate': dict(mpool.get('fate', {})),
                    })
                    self.medium_species_names.append(f"{additive_name}:{mpool['name']}")
                    self.medium_species_additive_index.append(a_idx)
        elif (
            self.initial_chemical_concs is not None
            and self.chemical_release_config is not None
            and self.chemical_release_data is not None
        ):
            # Legacy single-additive path normalized into the same internal
            # species representation.
            self.additive_names = ['additive_0']
            self.species_names = ['additive_0:pool_0']
            self.chemical_species = [{
                'additive_name': 'additive_0',
                'pool_name': 'pool_0',
                'initial_concs': np.array(self.initial_chemical_concs, dtype=float),
                'model': str(
                    self.chemical_release_config.get('model', 'analytical')
                ).lower(),
                'solver': dict(self.chemical_release_config.get('solver', {})),
                'params': dict(self.chemical_release_data or {}),
                'release_target': 'additive_0:medium',
                'inheritance': {'mode': 'proportional'},
                'fate': {},
            }]
            self.species_additive_index = [0]
            self.medium_species = [{
                'additive_name': 'additive_0',
                'pool_name': 'medium',
                'initial_mass': 0.0,
                'fate': {'k_deg': 0.0, 'k_loss': 0.0, 'transfers': []},
            }]
            self.medium_species_names = ['additive_0:medium']
            self.medium_species_additive_index = [0]

        self.species_additive_index = np.array(self.species_additive_index, dtype=int)
        self.medium_species_additive_index = np.array(self.medium_species_additive_index, dtype=int)
        self.n_additives = len(self.additive_names)
        self.n_chemical_species = len(self.chemical_species)
        self.n_medium_species = len(self.medium_species)

    def run(self) -> FMNPOutput:
        r"""
        Run the model with the config and data provided at initialisation.

        Returns
        -------
        :class:`fragmentmnp.output.FMNPOutput` object containing model output

        Notes
        -----
        The legacy polymer equations are now evaluated for each component
        ``j`` and size class ``k``:

        .. math::
            \frac{dc_{j,k}}{dt} =
            -k_{\mathrm{frag},j,k} c_{j,k}
            + \sum_i f_{j,i,k} k_{\mathrm{frag},j,i} c_{j,i}
            - k_{\mathrm{diss},j,k} c_{j,k}

        .. math::
            \frac{dc_{\mathrm{diss},j}}{dt} =
            \sum_k k_{\mathrm{diss},j,k} c_{j,k}
            - k_{\mathrm{min},j} c_{\mathrm{diss},j}

        .. math::
            \frac{dc_{\mathrm{min},j}}{dt} =
            k_{\mathrm{min},j} c_{\mathrm{diss},j}

        Aggregated legacy outputs are returned as sums over ``j``.  New
        component-resolved outputs are also stored in ``FMNPOutput``.
        """
        C = self.n_components
        N = self.n_size_classes

        # Build interpolators once; solve_ivp calls f many times.
        f_frag = interpolate.interp1d(self.t_grid, self.component_k_frag, axis=2,
                                      fill_value='extrapolate')
        f_diss = interpolate.interp1d(self.t_grid, self.component_k_diss, axis=2,
                                      fill_value='extrapolate')
        f_min = interpolate.interp1d(self.t_grid, self.component_k_min, axis=1,
                                     fill_value='extrapolate')

        def f(t, y):
            """
            Component-resolved initial value problem.

            State layout:
              y[:C*N]              = particulate polymer by component and size
              y[C*N:C*N+C]         = dissolved polymer by component
              y[C*N+C:C*N+2*C]     = mineralised polymer by component
            """
            c_particles = y[:C * N].reshape(C, N)
            c_dissolved = y[C * N:C * N + C]

            k_frag_t = f_frag(t)  # (C, N)
            k_diss_t = f_diss(t)  # (C, N)
            k_min_t = f_min(t)    # (C,)

            d_part = np.zeros((C, N), dtype=float)
            d_diss = np.zeros(C, dtype=float)
            d_min = np.zeros(C, dtype=float)

            for j in range(C):
                fsd_j = self.component_fsd[j]
                frag_loss = k_frag_t[j] * c_particles[j]
                diss_loss = k_diss_t[j] * c_particles[j]

                # Gains to size k from all fragmenting parent classes i.
                frag_gain = np.sum(fsd_j * frag_loss[:, np.newaxis], axis=0)

                d_part[j] = -frag_loss + frag_gain - diss_loss
                d_diss[j] = float(np.sum(diss_loss) - k_min_t[j] * c_dissolved[j])
                d_min[j] = float(k_min_t[j] * c_dissolved[j])

            return np.concatenate([d_part.ravel(), d_diss, d_min])

        y0 = np.concatenate([
            self.component_initial_concs.ravel(),
            self.component_initial_concs_diss,
            np.zeros(C, dtype=float),
        ])

        soln = solve_ivp(
            fun=f,
            method=self.config['solver_method'],
            t_span=(self.t_grid.min(), self.t_grid.max()),
            y0=y0,
            t_eval=self.t_eval,
            rtol=self.config['solver_rtol'],
            atol=self.config['solver_atol'],
            max_step=self.config['solver_max_step']
        )
        if not soln.success:
            raise FMNPNumericalError('Model solution could not be found: ' +
                                     f'{soln.message}')

        T = soln.y.shape[1]
        c_component_sol = soln.y[:C * N, :].reshape(C, N, T)
        c_diss_component_sol = soln.y[C * N:C * N + C, :]
        c_min_component_sol = soln.y[C * N + C:C * N + 2 * C, :]

        # Backward-compatible aggregate outputs.
        c_part_sol = np.sum(c_component_sol, axis=0)
        c_diss_sol = np.sum(c_diss_component_sol, axis=0)
        c_min_sol = np.sum(c_min_component_sol, axis=0)

        # Convert microparticle mass to particle number.  The legacy one-
        # component path deliberately calls mass_to_particle_number() so user
        # overloads of that method still behave as before.
        if C == 1:
            n_part_sol = self.mass_to_particle_number(c_part_sol)
            n_component_sol = n_part_sol[np.newaxis, :, :]
        else:
            n_component_sol = self._mass_to_particle_number_components(c_component_sol)
            n_part_sol = np.sum(n_component_sol, axis=0)

        layer_thickness_component, layer_thickness_by_size = self._calculate_layer_thickness(
            c_component_sol
        )

        # ------------------------------------------------------------
        # OPTIONAL: compute additive time series (post-solve)
        # ------------------------------------------------------------
        c_chem_part_species = None
        c_chem_medium_species = None
        c_chem_part_total = None
        c_chem_medium_total = None
        c_medium_pool_species = None

        if self.n_chemical_species > 0:
            (
                c_chem_part_species,
                c_chem_medium_species,
                c_chem_part_total,
                c_chem_medium_total,
                c_medium_pool_species
            ) = self._simulate_additives_postsolve(soln.t, c_part_sol)

        return FMNPOutput(
            t=soln.t,
            c=c_part_sol,
            n=n_part_sol,
            c_diss=c_diss_sol,
            c_min=c_min_sol,
            soln=soln,
            psd=self.psd,
            c_chem_part_species=c_chem_part_species,
            c_chem_medium_species=c_chem_medium_species,
            c_chem_part_total=c_chem_part_total,
            c_chem_medium_total=c_chem_medium_total,
            additive_names=self.additive_names,
            species_names=self.species_names,
            c_medium_pool_species=c_medium_pool_species,
            medium_pool_names=self.medium_species_names,
            c_component=c_component_sol,
            n_component=n_component_sol,
            c_diss_component=c_diss_component_sol,
            c_min_component=c_min_component_sol,
            component_names=self.component_names,
            layer_thickness_component=layer_thickness_component,
            layer_thickness_by_size=layer_thickness_by_size,
            particle_geometry=self.particle_geometry,
        )

    # ------------------------------------------------------------------
    # Component / multilayer helper methods
    # ------------------------------------------------------------------

    def _component_mass_weights(self) -> np.ndarray:
        """Return initial-mass weights for component aggregation."""
        m = np.sum(self.component_initial_concs, axis=1)
        total = float(np.sum(m))
        if total <= 0.0:
            return np.full(self.n_components, 1.0 / max(self.n_components, 1), dtype=float)
        return m / total

    def _component_weighted_rate(self, arr: np.ndarray) -> np.ndarray:
        """
        Convert component-resolved rates to a legacy aggregate rate array.

        ``arr`` can have shape (C, N, T) or (C, T).  The result is used for
        backward-compatible helper methods only; the core ODE never uses this
        aggregate in multi-component mode.
        """
        arr = np.asarray(arr, dtype=float)
        weights = self._component_mass_weights()
        if arr.ndim == 3:
            return np.tensordot(weights, arr, axes=(0, 0))
        if arr.ndim == 2:
            return np.tensordot(weights, arr, axes=(0, 0))
        raise ValueError('Component rate array must have shape (C,N,T) or (C,T).')

    def _mass_to_particle_number_components(self,
                                            mass_component: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """
        Convert component-resolved mass concentrations to number concentrations.

        This is density-resolved, so two layers with the same mass and size but
        different densities produce different particle-number contributions.
        """
        mass_component = np.asarray(mass_component, dtype=float)
        volume = self.volume(self.psd)
        denom = self.component_density[:, np.newaxis, np.newaxis] * volume[np.newaxis, :, np.newaxis]
        return mass_component / denom

    def _calculate_layer_thickness(self,
                                   c_component_sol: npt.NDArray[np.float64]) -> tuple[np.ndarray | None, np.ndarray | None]:
        """
        Track equivalent remaining thickness for explicitly layered components.

        For a component with initial layer thickness ``h_j(0)``, the bulk
        equivalent thickness is:

            h_j(t) = h_j(0) * M_j(t) / M_j(0)

        and the size-resolved diagnostic is:

            h_{j,k}(t) = h_j(0) * c_{j,k}(t) / c_{j,k}(0)

        where undefined ratios for initially empty size classes are returned as
        zero.  These are diagnostic layer-origin metrics: fragmentation can move
        layer-origin mass between size classes, while total component mass still
        controls the physically interpretable equivalent layer thickness.
        """
        if np.all(~np.isfinite(self.component_layer_initial_thickness)):
            return None, None

        C, N, T = c_component_sol.shape
        h_component = np.full((C, T), np.nan, dtype=float)
        h_by_size = np.full((C, N, T), np.nan, dtype=float)

        for j in range(C):
            h0 = self.component_layer_initial_thickness[j]
            if not np.isfinite(h0):
                continue

            m0 = float(np.sum(self.component_initial_concs[j]))
            mt = np.sum(c_component_sol[j], axis=0)
            if m0 > 0.0:
                h_component[j] = h0 * mt / m0
            else:
                h_component[j] = 0.0

            for k in range(N):
                c0 = float(self.component_initial_concs[j, k])
                if c0 > 0.0:
                    h_by_size[j, k] = h0 * c_component_sol[j, k] / c0
                else:
                    h_by_size[j, k] = 0.0

        return h_component, h_by_size

    # ---------------------------------------------------------------------
    # Additive coupling implementation (post-solve bookkeeping)
    # ---------------------------------------------------------------------

    def _simulate_additives_postsolve(self,
                                      t: npt.NDArray[np.float64],
                                      c_part_sol: npt.NDArray[np.float64]):
        """
        Compute additive trajectories using the already-solved polymer output.

        Phase 2 additions:
        - named medium pools per additive
        - release targets from particulate pools into named medium pools
        - medium-pool fate / transfers (for transformed-product generation etc.)

        Returns
        -------
        c_chem_part_species : array (S_part, N, T)
            particulate additive mass for each particulate species
        c_chem_medium_species : array (S_part, T)
            cumulative mass released from each particulate species to medium
            (kept for backward compatibility)
        c_chem_part_total : array (A, N, T)
            additive totals aggregated over particulate pools
        c_chem_medium_total : array (A, T)
            additive totals aggregated over all medium pools
        c_medium_pool_species : array (S_med, T)
            mass in each named medium pool
        """
        N, T = c_part_sol.shape
        S_part = self.n_chemical_species
        S_med = self.n_medium_species
        A = self.n_additives

        c_chem_part_species = np.zeros((S_part, N, T), dtype=float)
        c_chem_medium_species = np.zeros((S_part, T), dtype=float)
        c_medium_pool_species = np.zeros((S_med, T), dtype=float)

        particulate_name_to_index = {
            f"{spec['additive_name']}:{spec['pool_name']}": i
            for i, spec in enumerate(self.chemical_species)
        }
        medium_name_to_index = {
            f"{spec['additive_name']}:{spec['pool_name']}": i
            for i, spec in enumerate(self.medium_species)
        }

        normalized_particulate = []
        for spec in self.chemical_species:
            fate = dict(spec.get('fate', {}))
            transfers = []
            for tr in fate.get('transfers', []):
                target_name = tr['to']
                transfers.append({
                    'to': target_name,
                    'k': tr['k'],  # keep scalar or size-vector as provided
                    'target_species_index': int(particulate_name_to_index[target_name]),
                })

            fate['k_deg'] = fate.get('k_deg', 0.0)   # keep scalar or size-vector
            fate['k_loss'] = fate.get('k_loss', 0.0) # keep scalar or size-vector
            fate['transfers'] = transfers

            normalized_particulate.append({
                **spec,
                'fate': fate,
                'release_target_index': int(medium_name_to_index[spec['release_target']]),
            })

        normalized_medium = []
        for spec in self.medium_species:
            fate = dict(spec.get('fate', {}))
            transfers = []
            for tr in fate.get('transfers', []):
                target_name = tr['to']
                transfers.append({
                    'to': target_name,
                    'k': tr['k'],  # medium pools still usually scalar, but no need to force here
                    'target_species_index': int(medium_name_to_index[target_name]),
                })

            fate['k_deg'] = fate.get('k_deg', 0.0)
            fate['k_loss'] = fate.get('k_loss', 0.0)
            fate['transfers'] = transfers

            normalized_medium.append({**spec, 'fate': fate})

        for s, spec in enumerate(normalized_particulate):
            c_chem_part_species[s, :, 0] = spec['initial_concs']
            c_chem_medium_species[s, 0] = 0.0
        for m_idx, spec in enumerate(normalized_medium):
            c_medium_pool_species[m_idx, 0] = float(spec.get('initial_mass', 0.0))

        f_frag = interpolate.interp1d(
            self.t_grid, self.k_frag, axis=1, fill_value='extrapolate'
        )
        release_radii = self.release_radii
        eps = 1e-30
        release_cache = {}

        for ti in range(T - 1):
            dt_i = float(t[ti + 1] - t[ti])
            c = c_part_sol[:, ti]
            k_frag = f_frag(t[ti])
            L_frag = k_frag * c * dt_i
            G = (self.fsd.T * L_frag).T

            part_next = np.zeros((S_part, N), dtype=float)
            release_next = np.zeros(S_part, dtype=float)
            transfer_incoming = np.zeros((S_part, N), dtype=float)
            medium_release_incoming = np.zeros(S_med, dtype=float)

            for s, spec in enumerate(normalized_particulate):
                conc_A = np.zeros(N, dtype=float)
                mask = c > eps
                conc_A[mask] = c_chem_part_species[s, mask, ti] / c[mask]

                # Additive mass lost from each parent size class due to
                # polymer fragmentation during this timestep.
                A_loss = conc_A * L_frag

                # Redistribute that additive mass to daughter classes using the
                # selected inheritance mode.
                A_gain = np.zeros(N, dtype=float)
                inheritance = dict(spec.get('inheritance', {'mode': 'proportional'}))

                for i_parent in range(N):
                    lost_i = A_loss[i_parent]
                    if lost_i <= 0.0:
                        continue

                    w_daughters = self._fragmentation_inheritance_weights(
                        parent_size_index=i_parent,
                        inheritance=inheritance
                    )
                    A_gain += lost_i * w_daughters

                A_mid = c_chem_part_species[s, :, ti] - A_loss + A_gain

                # Evaluate physical release parameters at the current model time.
                # Backward-compatible scalar inputs become constant vectors, while
                # time-series/model dictionaries can change from one timestep to the next.
                release_params = {**spec['params'], **spec['solver']}
                evaluated_release_params = self._evaluate_release_params_for_size_classes(
                    params=release_params,
                    t_current=float(t[ti]),
                    n_size_classes=N,
                )

                cache_key = (
                    s,
                    float(dt_i),
                    spec['model'],
                    self._release_params_cache_key(evaluated_release_params),
                )
                if cache_key in release_cache:
                    rel = release_cache[cache_key]
                else:
                    if spec['model'] == 'analytical':
                        rel = np.array([
                            self._analytical_additive_release_fraction(
                                radius_m=float(release_radii[i]),
                                dt=dt_i,
                                params=self._release_params_for_size(evaluated_release_params, i),
                            ) for i in range(N)
                        ], dtype=float)
                    elif spec['model'] == 'numerical':
                        rel = np.array([
                            self._numerical_additive_release_fraction(
                                radius_m=float(release_radii[i]),
                                dt=dt_i,
                                params=self._release_params_for_size(evaluated_release_params, i),
                            ) for i in range(N)
                        ], dtype=float)
                    else:
                        raise ValueError(
                            f"Unknown chemical_release model '{spec['model']}' for "
                            f"{spec['additive_name']}:{spec['pool_name']}."
                        )
                    rel = np.clip(rel, 0.0, 1.0)
                    release_cache[cache_key] = rel

                dA_rel = rel * A_mid
                A_after_release = A_mid - dA_rel
                release_next[s] = c_chem_medium_species[s, ti] + float(dA_rel.sum())
                medium_release_incoming[spec['release_target_index']] += float(dA_rel.sum())

                A_after_fate, transfer_out = self._apply_first_order_fate_explicit(
                    A_in=A_after_release, dt=dt_i, spec=spec
                )
                part_next[s] = A_after_fate
                for target_idx, arr in transfer_out.items():
                    transfer_incoming[target_idx] += arr

            part_next += transfer_incoming
            for s in range(S_part):
                c_chem_part_species[s, :, ti + 1] = np.clip(part_next[s], 0.0, None)
                c_chem_medium_species[s, ti + 1] = max(release_next[s], 0.0)

            medium_next = c_medium_pool_species[:, ti].copy() + medium_release_incoming
            medium_transfer_incoming = np.zeros(S_med, dtype=float)
            medium_after_fate = np.zeros(S_med, dtype=float)
            for m_idx, spec in enumerate(normalized_medium):
                out_mass, transfer_out = self._apply_first_order_fate_explicit(
                    A_in=np.array([medium_next[m_idx]], dtype=float), dt=dt_i, spec=spec
                )
                medium_after_fate[m_idx] = float(out_mass[0])
                for target_idx, arr in transfer_out.items():
                    medium_transfer_incoming[target_idx] += float(np.asarray(arr)[0])
            medium_after_fate += medium_transfer_incoming
            c_medium_pool_species[:, ti + 1] = np.clip(medium_after_fate, 0.0, None)

        c_chem_part_total = np.zeros((A, N, T), dtype=float)
        c_chem_medium_total = np.zeros((A, T), dtype=float)
        for s, a_idx in enumerate(self.species_additive_index):
            c_chem_part_total[a_idx] += c_chem_part_species[s]
        for m_idx, a_idx in enumerate(self.medium_species_additive_index):
            c_chem_medium_total[a_idx] += c_medium_pool_species[m_idx]

        return (
            c_chem_part_species,
            c_chem_medium_species,
            c_chem_part_total,
            c_chem_medium_total,
            c_medium_pool_species,
        )


    def _evaluate_release_params_for_size_classes(self,
                                                  params: dict,
                                                  t_current: float,
                                                  n_size_classes: int) -> dict:
        """
        Evaluate release parameters at the current model time.

        Backward-compatible scalar parameters remain constant. New supported
        time-dependent forms include:

        1) Time series interpolation:
           {'times': [0, 10, 20], 'values': [1e-16, 2e-16, 5e-16]}
           or size-resolved values with shape (len(times), n_size_classes).

        2) Simple weathering/ageing models:
           {'type': 'linear', 'initial': 1e-16, 'rate': 1e-3}
           {'type': 'exponential', 'initial': 1e-16, 'rate': 1e-5}
           {'type': 'logistic', 'initial': 1e-16, 'factor': 10, 't_mid': 50, 'steepness': 0.1}

        3) FRAGMENT-MNP distribution dictionaries with k_f/k_0 and optional
           t/s regression terms. These are evaluated over the model t_grid and
           interpolated to t_current.
        """
        out = dict(params)
        for key in ['D_p', 'D_w', 'K_pw', 'k_m']:
            if key in out:
                out[key] = self._evaluate_release_parameter(
                    name=key,
                    value=out[key],
                    t_current=t_current,
                    n_size_classes=n_size_classes,
                )
        return out

    @staticmethod
    def _release_params_for_size(params: dict, size_index: int) -> dict:
        """Extract scalar parameter values for one particle size class."""
        out = dict(params)
        for key in ['D_p', 'D_w', 'K_pw', 'k_m']:
            if key in out:
                val = out[key]
                if isinstance(val, np.ndarray):
                    out[key] = float(val[size_index])
                elif isinstance(val, (list, tuple)):
                    out[key] = float(val[size_index])
        return out

    @staticmethod
    def _release_params_cache_key(params: dict) -> tuple:
        """Build a hashable cache key from evaluated release and solver params."""
        items = []
        for key in sorted(params.keys()):
            val = params[key]
            if isinstance(val, np.ndarray):
                items.append((key, tuple(np.round(val.astype(float), 30))))
            elif isinstance(val, (list, tuple)):
                try:
                    arr = np.asarray(val, dtype=float)
                    items.append((key, tuple(np.round(arr, 30))))
                except Exception:
                    items.append((key, tuple(val)))
            elif isinstance(val, (int, float, str, bool)) or val is None:
                items.append((key, val))
            else:
                items.append((key, repr(val)))
        return tuple(items)

    def _evaluate_release_parameter(self,
                                    name: str,
                                    value,
                                    t_current: float,
                                    n_size_classes: int) -> np.ndarray:
        """Return a length-n_size_classes vector for one release parameter."""
        if isinstance(value, (int, float)):
            out = np.full(n_size_classes, float(value), dtype=float)
            self._check_evaluated_release_parameter(name, out)
            return out

        if isinstance(value, dict):
            if 'times' in value or 'values' in value:
                out = self._evaluate_release_parameter_timeseries(
                    value=value,
                    t_current=t_current,
                    n_size_classes=n_size_classes,
                )
                self._check_evaluated_release_parameter(name, out)
                return out

            if 'k_f' in value:
                out = self._evaluate_release_parameter_distribution(
                    value=value,
                    t_current=t_current,
                    n_size_classes=n_size_classes,
                )
                self._check_evaluated_release_parameter(name, out)
                return out

            out = self._evaluate_release_parameter_model(
                value=value,
                t_current=t_current,
                n_size_classes=n_size_classes,
            )
            self._check_evaluated_release_parameter(name, out)
            return out

        arr = np.asarray(value, dtype=float)
        if arr.ndim == 1:
            if arr.size != n_size_classes:
                raise ValueError(
                    f"Release parameter '{name}' vector must have length {n_size_classes}."
                )
            out = arr.astype(float, copy=False)
            self._check_evaluated_release_parameter(name, out)
            return out

        if arr.ndim == 2:
            # Compact matrix form. One dimension must be size; the other is
            # assumed to correspond to self.t_grid.
            if arr.shape == (n_size_classes, self.t_grid.size):
                mat = arr
            elif arr.shape == (self.t_grid.size, n_size_classes):
                mat = arr.T
            else:
                raise ValueError(
                    f"Release parameter '{name}' 2D matrix must have shape "
                    f"({n_size_classes}, {self.t_grid.size}) or "
                    f"({self.t_grid.size}, {n_size_classes})."
                )
            out = np.array([
                np.interp(t_current, self.t_grid, mat[i])
                for i in range(n_size_classes)
            ], dtype=float)
            self._check_evaluated_release_parameter(name, out)
            return out

        raise ValueError(
            f"Release parameter '{name}' must be scalar, vector, matrix, or dict."
        )

    def _evaluate_release_parameter_timeseries(self,
                                               value: dict,
                                               t_current: float,
                                               n_size_classes: int) -> np.ndarray:
        """Evaluate {'times': ..., 'values': ...} release-parameter input."""
        times = np.asarray(value['times'], dtype=float)
        values = np.asarray(value['values'], dtype=float)

        if values.ndim == 1:
            val = float(np.interp(t_current, times, values))
            return np.full(n_size_classes, val, dtype=float)

        if values.shape == (times.size, n_size_classes):
            # values[time, size]
            return np.array([
                np.interp(t_current, times, values[:, i])
                for i in range(n_size_classes)
            ], dtype=float)

        if values.shape == (n_size_classes, times.size):
            # values[size, time]
            return np.array([
                np.interp(t_current, times, values[i, :])
                for i in range(n_size_classes)
            ], dtype=float)

        raise ValueError(
            "Time-dependent release parameter values must be 1D, "
            "(len(times), n_size_classes), or (n_size_classes, len(times))."
        )

    def _evaluate_release_parameter_distribution(self,
                                                 value: dict,
                                                 t_current: float,
                                                 n_size_classes: int) -> np.ndarray:
        """
        Evaluate a FRAGMENT-MNP-style t/s distribution for a release parameter.
        """
        k_f = float(value.get('k_f'))
        k_0 = float(value.get('k_0', 0.0))
        is_compound = bool(value.get('is_compound', True))
        reg_params = {
            k: v for k, v in value.items()
            if k not in ['k_f', 'k_0', 'is_compound']
        }
        grid = self.set_k_distribution(
            dims={'s': self.surface_areas, 't': self.t_grid},
            k_f=k_f,
            k_0=k_0,
            params=reg_params,
            is_compound=is_compound,
        )
        if grid.shape != (n_size_classes, self.t_grid.size):
            raise ValueError("Unexpected release-parameter distribution shape.")
        return np.array([
            np.interp(t_current, self.t_grid, grid[i])
            for i in range(n_size_classes)
        ], dtype=float)

    def _evaluate_release_parameter_model(self,
                                          value: dict,
                                          t_current: float,
                                          n_size_classes: int) -> np.ndarray:
        """Evaluate simple parametric ageing/weathering release models."""
        model = str(value.get('type', value.get('model', 'constant'))).lower()
        base = value.get('initial', value.get('base', value.get('value', None)))
        if base is None:
            raise ValueError(
                "Release parameter model dict must contain 'initial', 'base', or 'value'."
            )

        base_arr = self._expand_release_value_to_size_vector(base, n_size_classes)
        t0 = float(value.get('t0', value.get('time_origin', 0.0)))
        tau = float(t_current) - t0

        if model == 'constant':
            out = base_arr.copy()

        elif model == 'linear':
            if 'slope' in value:
                slope = self._expand_release_value_to_size_vector(value['slope'], n_size_classes)
                out = base_arr + slope * tau
            else:
                rate = self._expand_release_value_to_size_vector(value.get('rate', 0.0), n_size_classes)
                out = base_arr * (1.0 + rate * tau)

        elif model == 'exponential':
            rate = self._expand_release_value_to_size_vector(value.get('rate', 0.0), n_size_classes)
            out = base_arr * np.exp(rate * tau)

        elif model == 'logistic':
            # Smooth transition from base to base*factor around t_mid.
            factor = self._expand_release_value_to_size_vector(value.get('factor', 1.0), n_size_classes)
            t_mid = self._expand_release_value_to_size_vector(value.get('t_mid', 0.0), n_size_classes)
            steepness = self._expand_release_value_to_size_vector(value.get('steepness', 1.0), n_size_classes)
            out = base_arr * (1.0 + (factor - 1.0) / (1.0 + np.exp(-steepness * (tau - t_mid))))

        else:
            raise ValueError(
                "Release parameter model/type must be one of "
                "'constant', 'linear', 'exponential', or 'logistic'."
            )

        if 'size_factor' in value:
            out = out * self._expand_release_value_to_size_vector(value['size_factor'], n_size_classes)
        if 'minimum' in value:
            out = np.maximum(out, self._expand_release_value_to_size_vector(value['minimum'], n_size_classes))
        if 'maximum' in value:
            out = np.minimum(out, self._expand_release_value_to_size_vector(value['maximum'], n_size_classes))

        return out.astype(float, copy=False)

    @staticmethod
    def _expand_release_value_to_size_vector(value, n_size_classes: int) -> np.ndarray:
        """Expand scalar or length-n_size_classes value to a vector."""
        if isinstance(value, (int, float)):
            return np.full(n_size_classes, float(value), dtype=float)
        arr = np.asarray(value, dtype=float)
        if arr.ndim != 1 or arr.size != n_size_classes:
            raise ValueError(
                f"Expected scalar or length-{n_size_classes} value."
            )
        return arr.astype(float, copy=False)

    @staticmethod
    def _check_evaluated_release_parameter(name: str, values: np.ndarray) -> None:
        """Runtime guard for evaluated release parameters."""
        values = np.asarray(values, dtype=float)
        if np.any(~np.isfinite(values)):
            raise ValueError(f"Release parameter '{name}' evaluated to non-finite values.")
        if name in ['D_p', 'D_w', 'K_pw']:
            if np.any(values <= 0.0):
                raise ValueError(f"Release parameter '{name}' must evaluate to positive values.")
        elif name == 'k_m':
            if np.any(values < 0.0):
                raise ValueError("Release parameter 'k_m' must evaluate to non-negative values.")


    def _apply_first_order_fate_explicit(self,
                                         A_in: npt.NDArray[np.float64],
                                         dt: float,
                                         spec: dict) -> Tuple[npt.NDArray[np.float64], dict]:
        """
        Apply explicit first-order fate over one timestep.

        For particulate species:
        - k_deg, k_loss, and transfer k may be scalars or size-class vectors.

        For medium pools:
        - A_in has length 1 and rates are treated as scalars.
        """
        A_in = np.asarray(A_in, dtype=float)
        A_work = np.clip(A_in, 0.0, None).copy()
        n_local = A_work.size

        fate = dict(spec.get('fate', {}))

        k_deg = self._expand_rate_to_size_vector(fate.get('k_deg', 0.0), n_local)
        k_loss = self._expand_rate_to_size_vector(fate.get('k_loss', 0.0), n_local)
        transfers = list(fate.get('transfers', []))

        k_transfer_list = []
        for tr in transfers:
            k_tr = self._expand_rate_to_size_vector(tr.get('k', 0.0), n_local)
            k_transfer_list.append(k_tr)

        total_k = k_deg + k_loss
        for k_tr in k_transfer_list:
            total_k = total_k + k_tr

        if np.any(total_k < 0.0):
            raise SchemaError(
                f"Negative total fate rate encountered for "
                f"{spec['additive_name']}:{spec['pool_name']}."
            )

        if dt <= 0.0 or np.all(total_k == 0.0):
            return A_work, {}

        frac_deg = np.minimum(k_deg * dt, 1.0)
        frac_loss = np.minimum(k_loss * dt, 1.0)
        frac_transfers = [np.minimum(k_tr * dt, 1.0) for k_tr in k_transfer_list]

        frac_total = frac_deg + frac_loss
        for ftr in frac_transfers:
            frac_total = frac_total + ftr

        # Prevent per-size overshoot
        scale = np.ones(n_local, dtype=float)
        mask = frac_total > 1.0
        scale[mask] = 1.0 / frac_total[mask]

        frac_deg = frac_deg * scale
        frac_loss = frac_loss * scale
        frac_transfers = [ftr * scale for ftr in frac_transfers]

        transfer_out = {}
        for tr, frac in zip(transfers, frac_transfers):
            moved = frac * A_work
            target_idx = int(tr['target_species_index'])
            transfer_out[target_idx] = transfer_out.get(target_idx, 0.0) + moved

        retained_fraction = 1.0 - frac_deg - frac_loss
        for ftr in frac_transfers:
            retained_fraction = retained_fraction - ftr
        retained_fraction = np.clip(retained_fraction, 0.0, 1.0)

        A_out = retained_fraction * A_work

        return np.clip(A_out, 0.0, None), transfer_out

    def _analytical_additive_release_fraction(self, radius_m: float, dt: float, params: dict) -> float:
        """
        Return the fraction of additive released from a particle of radius_m
        over timestep dt.

        This is what your post-solve bookkeeping loop calls.

        Required params in params dict after time evaluation:
            - "D_p" : polymer diffusivity at this timestep [m^2/s]
            - "D_w" : water diffusivity at this timestep [m^2/s]
            - "K_pw": polymer-water partition coefficient at this timestep [-]

        Optional:
            - "n_terms": int number of series terms (Bi>=100)

        Returns
        -------
        float
            Release fraction in [0,1]
        """
        D_p = float(params.get("D_p"))
        D_w = float(params.get("D_w"))
        K_pw = float(params.get("K_pw"))
        n_terms = int(params.get("n_terms", 50))

        if self.geometry.release_geometry == 'sphere':
            F_remain = self._analytical_additive_remaining_fraction(
                t=float(dt),
                r=float(radius_m),
                D_p=D_p,
                D_w=D_w,
                K_pw=K_pw,
                n_terms=n_terms
            )
        elif self.geometry.release_geometry == 'cylinder':
            F_remain = self._analytical_cylindrical_remaining_fraction(
                t=float(dt),
                r=float(radius_m),
                D_p=D_p,
                D_w=D_w,
                K_pw=K_pw,
                n_terms=n_terms
            )
        else:
            raise ValueError(
                f"Unsupported additive release geometry: {self.geometry.release_geometry}"
            )

        # Release fraction = 1 - remaining fraction (bounded)
        rel = 1.0 - F_remain
        return float(np.clip(rel, 0.0, 1.0))
    
    
    def _numerical_additive_release_fraction(self, radius_m: float, dt: float, params: dict) -> float:
        """
        Release fraction computed by numerically solving diffusion in a sphere.

        Required after time evaluation:
            D_p, K_pw
        Either after time evaluation:
            k_m directly
        OR:
            D_w (and we compute k_m = K_pw * D_w / r for consistency with analytical assumptions)

        Optional:
            n_r (radial resolution, default 40)
        """
        D_p = float(params.get("D_p"))
        K_pw = float(params.get("K_pw"))
        n_r = int(params.get("n_r", 60))
        n_substeps = int(params.get("n_substeps", 20))
        theta = float(params.get("theta", 1.0))

        # Preferred: k_m explicitly supplied
        k_m = params.get("k_m", None)

        if k_m is None:
            # Backward-compatible derivation from your previous notes:
            # k_m = K_pw * D_w / r
            D_w = float(params.get("D_w"))
            k_m = K_pw * D_w / max(radius_m, 1e-30)
        else:
            k_m = float(k_m)

        n_substeps = int(params.get("n_substeps", 20))

        if self.geometry.release_geometry == 'sphere':
            F_remain = self._numerical_additive_remaining_fraction(
                t=float(dt),
                r=float(radius_m),
                D_p=D_p,
                K_pw=K_pw,
                k_m=k_m,
                n_r=n_r,
                n_substeps=n_substeps,
                theta=theta
            )
        elif self.geometry.release_geometry == 'cylinder':
            F_remain = self._numerical_cylindrical_remaining_fraction(
                t=float(dt),
                r=float(radius_m),
                D_p=D_p,
                K_pw=K_pw,
                k_m=k_m,
                n_r=n_r,
                n_substeps=n_substeps,
                theta=theta
            )
        else:
            raise ValueError(
                f"Unsupported additive release geometry: {self.geometry.release_geometry}"
            )
        rel = 1.0 - F_remain
        return float(np.clip(rel, 0.0, 1.0))
     

    def mass_to_particle_number(self, mass):
        """
        Convert mass (concentration) to particle number (concentration).
        """
        return mass / (self.density * self.volume(self.psd))[:, np.newaxis]

    def _set_psd(self) -> npt.NDArray[np.float64]:
        """
        Calculate the particle size distribution based on
        the config options passed in `self.config`
        """
        if 'particle_size_classes' in self.config:
            # If a particle size distribution has been specified, use that
            psd = np.array(self.config['particle_size_classes'])
        elif 'particle_size_range' in self.config:
            # Else, construct the size distribution from the given size range
            psd = np.logspace(*self.config['particle_size_range'],
                              self.n_size_classes)
        else:
            raise KeyError(
                "Missing particle size configuration. Expected either "
                "'particle_size_classes' or 'particle_size_range'."
            )
        return psd

    def set_rate_constant(self, params: dict, name: str) \
            -> npt.NDArray[np.float64]:
        """
        Set a rate constant either for fragmentation, dissolution or
        mineralisation. The function largely parses input data and
        deals with defaults etc, before calling `set_k_distribution`
        to create the rate constant distribution.

        Parameters
        ----------
        params : dict
            The param dict from the model input data, from which the rate
            constant parameters are extracted. See the notes in
            `set_k_distribution` for more details on the expected format
        name : str
            The name of the rate constant to set. This should be one of
            `k_frag`, `k_diss` or `k_min`

        Returns
        -------
        np.ndarray
            The rate constant distribution
        """
        # If we've been given a dict in data, use the contained params to
        # create the distribution. Else, presume that we've been given a
        # scalar (validation will make sure this is so) and use that as
        # the average.
        if isinstance(params, dict):
            k_f = params['k_f']
            k_0 = params['k_0']
            is_compound = params['is_compound']
            # Get the params from the dict, excluding the average
            reg_params = {n: p for n, p in params.items()
                            if n not in ['k_f', 'k_0']}
        else:
            k_f = params
            k_0 = 0.0
            is_compound = True
            reg_params = {}
        # If k_frag or k_diss, we want to calculate a 2D (s, t) distribution,
        # and if k_min, we just want a 1D (t) distribution
        if name in ['k_frag', 'k_diss']:
            k_dist = self.set_k_distribution(dims={'s': self.surface_areas,
                                                   't': self.t_grid},
                                             k_f=k_f, k_0=k_0,
                                             params=reg_params,
                                             is_compound=is_compound)
            # If the rate constant is k_frag, then no fragmentation is
            # allowed from the smallest size class and therefore we
            # manually set this to zero
            if name == 'k_frag':
                k_dist[0, :] = 0.0
        else:
            k_dist = self.set_k_distribution(dims={'t': self.t_grid},
                                             k_f=k_f, k_0=k_0,
                                             params=reg_params,
                                             is_compound=is_compound)
        # Check no values are less than zero
        if np.any(k_dist < 0.0):
            msg = (f'Value for {name} distribution calculated from input '
                   'data resulted in negative values. Ensure '
                   f'distribution params are such that all {name} values '
                   'are positive.')
            raise FMNPDistributionValueError(msg)

        # Return this distribution
        return k_dist
    
    def _fragmentation_inheritance_weights(self,
                                           parent_size_index: int,
                                           inheritance: dict) -> np.ndarray:
        """
        Return daughter weighting factors for additive inheritance from one
        parent size class during a fragmentation event.

        The returned vector has length n_size_classes and is only non-zero for
        valid daughter classes (< parent_size_index). It is normalised to sum
        to 1 over valid daughters.

        Modes
        -----
        proportional:
            Uses the polymer fragment size distribution row directly.

        size_biased:
            Reweights daughters by d^beta on top of the polymer FSD row.

        surface_enriched:
            Reweights daughters by surface_area^gamma on top of the polymer
            FSD row.
        """
        mode = str(inheritance.get("mode", "proportional")).lower()

        w = np.array(self.fsd[parent_size_index], dtype=float)
        if parent_size_index <= 0 or np.sum(w) <= 0.0:
            return w

        daughter_mask = np.arange(self.n_size_classes) < parent_size_index

        if mode == "proportional":
            pass

        elif mode == "size_biased":
            beta = float(inheritance.get("beta", 0.0))
            bias = np.zeros_like(w)
            bias[daughter_mask] = self.psd[daughter_mask] ** beta
            w = w * bias

        elif mode == "surface_enriched":
            gamma = float(inheritance.get("gamma", 1.0))
            bias = np.zeros_like(w)
            bias[daughter_mask] = self.surface_areas[daughter_mask] ** gamma
            w = w * bias

        else:
            raise ValueError(f"Unknown inheritance mode '{mode}'.")

        s = float(np.sum(w))
        if s > 0.0:
            w = w / s
        return w

    @staticmethod
    def set_k_distribution(dims: dict, k_f: float, k_0: float = 0.0,
                           params: dict = {},
                           is_compound: bool = True) -> npt.NDArray[np.float64]:
        r"""
        Create a distribution based on the rate constant scaling factor ``k_f``
        and baseline adjustment factor ``k_0``. The distribution will be a
        compound or additive combination of power law / polynomial,
        exponential, logarithmic and logistic regressions, encapsulated in the
        function :math:`X(x)`, and have dimensions given by `dims`. For a
        distribution with `D` dimensions:

        .. math::
            k(\mathbf{x}) = k_f \prod_{d=1}^D X(x_d) + k_0

        or

        .. math::
            k(\mathbf{x}) = k_f \sum_{d=1}^D X(x_d) + k_0

        :math:`X(x)` is then given either by:

        .. math::
            X(x) = A_x \hat{x}^{\alpha_x} \cdot B_x e^{-\beta_x \hat{x}}
            \cdot C_x \ln (\gamma_x \hat{x}) \cdot
            \frac{D_x}{1 + e^{-\delta_{x,1}(\hat{x} - \delta_{x,2})}}

        or the user can specify a polynomial instead of the power law term:

        .. math::
            X(x) = \sum_{n=1}^N A_{x,n} \hat{x}^n \cdot
            B_x e^{-\beta_x \hat{x}} \cdot
            C_x \ln (\gamma_x \hat{x}) \cdot
            \frac{D_x}{1 + e^{-\delta_{x,1}(\hat{x} - \delta_{x,2})}}

        In the above, the dimension value :math:`\hat{x}` is normalised such
        that the median value is equal to 1: :math:`\hat{x} = x/\tilde{x}`.

        Parameters
        ----------
        dims : dict
            A dictionary that maps dimension names to their grids, e.g. to
            create a distribution of time `t` and particle surface area `s`,
            `dims` would equal `{'t': t, 's': s}`, where `t` and `s` are the
            timesteps and particle surface area bins over which to create this
            distribution. The dimension names must correspond to the subscripts
            used in `params`. The values are normalised such that the median
            of each dimension is 1.
        k_f : float
            Rate constant scaling factor
        k_0 : float, default=0
            Rate constant baseline adjustment factor
        params : dict, default={}
            A dictionary of values to parameterise the distribution with. See
            the notes below.
        is_compound : bool, default=True
            Whether the regression for each dimension are combined by
            multiplying (compound) or adding.

        Returns
        -------
        k = np.ndarray
            Distribution array over the dims provided

        Notes
        -----
        `k` is modelled as a function of the dims provided, and the model
        builds this distribution as a combination of power law / polynomial,
        exponential, logarithmic and logistic regressions, enabling a broad
        range of dependencies to be accounted for. This distribution is
        intended to be applied to rate constants used in the model, such as
        `k_frag` and `k_diss`. The `params` dict gives the parameters used
        to construct this distribution using the equation above. That is,
        :math:`A_{x}` (where `x` is the dimension), :math:`\alpha_{x_i}` etc
        are given in the `params` dict as e.g. `A_t`, `alpha_t`, where the
        subscript (`t` in this case) is the name of the dimension corresponding
        to the `dims` dict.

        This function does not require any parameters to be present in
        `params`. Non-present values are defaulted to values that remove the
        influence of that particular expression, and letting all parameters
        default results in a constant `k` distribution.

        More specifically, the params that can be specified are:

        A_x : array-like or float, default=1
            Power law coefficient(s) for dim `x` (e.g. ``A_t`` for dim `t`).
            If a scalar is provided, this is used as the coefficient for a
            power law expression with ``alpha_x`` as the exponent. If a list is
            provided, these are used as coefficients in a polynomial
            expression, where the order of the polynomial is given by the
            length of the list. For example, if a length-2 list ``A_t=[2, 3]``
            is given, then the resulting polynomial will be :math:`3t^2 + 2t`
            (note the list is in *ascending* order of polynomials).
        alpha_x : float, default=0
            If `A_x` is a scalar, `alpha_x` is the exponent for this power law
            expression. For example, if ``A_t=2`` and ``alpha_t=0.5``, the
            resulting power law will be :math:`2t^{0.5}`.
        B_x : float, default=1
            Exponential coefficient.
        beta_x : float, default=0
            Exponential scaling factor.
        C_x : float or None, default=None
            If a scalar is given, this is the coefficient for the logarithmic
            expression. If ``None`` is given, the logarithmic expression is
            set to 1 (i.e. it is ignored).
        gamma_x : float, default=1
            Logarithmic scaling factor.
        D_x : float or None, default=None
            If a scalar is given, this is the coefficient for the logistic
            expression. If `None` is given, the logistic expression is set
            to 1.
        delta1_x : float, default=1
            Logistic growth rate (steepness of the logistic curve).
        delta2_x : float or None, default=None
            Midpoint of the logistic curve, which denotes the `x` value where
            the logistic curve is at its midpoint. If `None` is given, the
            midpoint is assumed to be the at the midpoint of the `x` range. For
            example, for the time dimension `t`, if the model timesteps go
            from 1 to 100, then the default is ``delta2_t=50``.

        If any dimension values are equal to 0, the logarithmic term returns 0
        rather than being undefined.

        .. warning:: The parameters used for calculating distributions such as
            `k_frag` and `k_diss` have changed from previous versions, which
            only allowed for a power law relationship. This causes breaking
            changes after v0.1.0.
        """
        if dims == {}:
            raise Exception('Trying to create k distribution but `dims` dict',
                            'is empty. You must provide at least one',
                            'dimension.')
        else:
            # Create a grid out of our dimensions (and preserve the matrix
            # indexing order with 'ij')
            grid = np.meshgrid(*[x for x in dims.values()],
                               indexing='ij')
            # Assign the relevant params to each dimension and get a
            # list of regressions across the dimensions, X
            X = FragmentMNP._assign_regression_params(params,
                                                      grid,
                                                      list(dims.keys()))
            # Calculate the final distribution by multiplying or summing
            # X across the dimensions
            if is_compound:
                k = k_f * np.prod(X, axis=0) + k_0
            else:
                k = k_f * np.sum(X, axis=0) + k_0
            return k
        
    @staticmethod
    def _factor_tridiagonal(lower: np.ndarray,
                            diag: np.ndarray,
                            upper: np.ndarray):
        """
        Factor a tridiagonal matrix for repeated solves using the Thomas algorithm.

        Returns modified copies of lower/diag/upper containing the factorization.
        """
        n = len(diag)
        a = lower.astype(float, copy=True)
        b = diag.astype(float, copy=True)
        c = upper.astype(float, copy=True)

        for i in range(1, n):
            if b[i - 1] == 0.0:
                raise ZeroDivisionError("Tridiagonal factorization encountered zero pivot.")
            w = a[i] / b[i - 1]
            a[i] = w
            b[i] = b[i] - w * c[i - 1]

        return a, b, c

    @staticmethod
    def _solve_tridiagonal_factored(a_fact: np.ndarray,
                                    b_fact: np.ndarray,
                                    c_fact: np.ndarray,
                                    rhs: np.ndarray) -> np.ndarray:
        """
        Solve a tridiagonal system using a precomputed Thomas factorization.
        """
        n = len(b_fact)
        d = rhs.astype(float, copy=True)

        for i in range(1, n):
            d[i] = d[i] - a_fact[i] * d[i - 1]

        x = np.zeros(n, dtype=float)
        if b_fact[-1] == 0.0:
            raise ZeroDivisionError("Tridiagonal solver encountered zero pivot at last row.")
        x[-1] = d[-1] / b_fact[-1]

        for i in range(n - 2, -1, -1):
            if b_fact[i] == 0.0:
                raise ZeroDivisionError("Tridiagonal solver encountered zero pivot.")
            x[i] = (d[i] - c_fact[i] * x[i + 1]) / b_fact[i]

        return x

    @staticmethod
    def set_fsd(n: int,
                psd: npt.NDArray[np.float64],
                beta: float) -> npt.NDArray[np.float64]:
        r"""
        Set the fragment size distribution matrix, assuming that
        fragmentation events result in a split in mass between daughter
        fragments that scales proportionally to :math:`d^\beta`,
        where :math:`d` is the particle diameter and :math:`\beta`
        is an empirical fragment size distribution parameter. For
        example, if :math:`\beta` is negative, then a larger
        proportion of the fragmenting mass goes to smaller size
        classes than larger.

        For an equal split between daughter size classes, set
        :math:`\beta` to 0.

        Parameters
        ----------
        n : int
            Number of particle size classes
        psd : np.ndarray
            Particle size distribution
        beta : float
            Fragment size distribution empirical parameter

        Returns
        -------
        np.ndarray
            Matrix of fragment size distributions for all size classes
        """
        # Start with a zero-filled array of shape (N,N)
        fsd = np.zeros((n, n))
        # Fill with the split to daughter size classes scaled
        # proportionally to d^beta
        for i in np.arange(1, n):
            fsd[i, :-(n-i)] = (psd[:-(n-i)] ** beta
                               / np.sum(psd[:-(n-i)] ** beta))
        return fsd

    def surface_area(self, psd: npt.NDArray[np.float64]) -> \
            npt.NDArray[np.float64]:
        """Return particle surface area using the configured geometry."""
        return self.geometry.surface_area(psd)

    def volume(self, psd: npt.NDArray[np.float64]) -> \
            npt.NDArray[np.float64]:
        """Return particle volume using the configured geometry."""
        return self.geometry.volume(psd)

    @staticmethod
    def _assign_regression_params(params: dict,
                                  grid: Sequence[npt.NDArray],
                                  dim_names: list) \
            -> Sequence[npt.NDArray]:
        """
        Given a list of regression parameters, a grid and dimension names,
        return a list of the regressions for each dimensions.

        Parameters
        ----------
        params : dict
            The dict of parameters provided as input data
        grid : list(np.ndarray)
            The grid of values for each dimension, broadcast to the
            correct number of dimensions
        dim_names : list
            The names of the dimensions

        Returns
        -------
        X : list(np.ndarray)
            List of regressions for each dimension
        """
        # List of the regressions for each dimension, which will be
        # populated when we loop over the dimensions
        X = []
        # Loop over the dimensions
        for i, x in enumerate(grid):
            # Normalise the values
            x_norm = x / np.median(x)
            # Get the name of this dim
            name = dim_names[i]
            # Pull out the params for convenience
            A = params.get(f'A_{name}', 1.0)
            alpha = float(params.get(f'alpha_{name}', 0.0))
            B = float(params.get(f'B_{name}', 1.0))
            beta = float(params.get(f'beta_{name}', 0.0))
            C = params.get(f'C_{name}', None)
            gamma = float(params.get(f'gamma_{name}', 1.0))
            D = params.get(f'D_{name}', None)
            delta_1 = float(params.get(f'delta1_{name}', 1.0))
            delta_2 = params.get(f'delta2_{name}', None)
            # Let users specify a polynomial by listing A coefficients,
            # presuming the exponents (alpha) will be 1, 2, 3 etc
            # corresponding to the list elements in A. Otherwise, use
            # A and alpha as a power law A*x**alpha
            if isinstance(A, (list, tuple, np.ndarray)):
                powers = np.arange(1, len(A) + 1)
                power_x = 0.0
                for i, power in enumerate(powers):
                    power_x = power_x + A[i] * x_norm**power
            else:
                power_x = float(A) * x_norm**alpha
            # Exponential term
            exp_x = B * np.exp(-beta*x_norm)
            # Only calculate the ln term if C is not None, and if x=0,
            # then set ln(x) to 0
            if C is not None:
                arg = gamma * x_norm
                ln_term = np.log(arg,
                                 out=np.zeros_like(arg, dtype=np.float64),
                                 where=(arg != 0))
                ln_x = float(C) * ln_term
            else:
                ln_x = 1.0
            # If the logistic delta_2 term (the x value at the midpoint
            # along the k axis) isn't specified, # then calculate it as
            # halfway along the x axis
            if (delta_2 is None) and (D is not None):
                delta_2 = (x_norm.max() - x_norm.min()) / 2
            # Calculate the logistic contribution, only if D is not None
            logit_x = float(D) / \
                (1 + np.exp(-delta_1 * (x_norm - float(delta_2)))) \
                if D is not None else 1.0
            # Multiply all the expressions together for this dimension
            X.append(power_x * exp_x * ln_x * logit_x)
        # Return the list of regressions for each dimension
        return X

    def _f_surface_area(self,
                        psd: npt.NDArray[np.float64],
                        gamma: float = 1.0) -> npt.NDArray[np.float64]:
        r"""
        Geometry-aware surface-area-to-volume scaling factor.

        The factor is the particle surface-area-to-volume ratio normalized to
        the median size class and raised to ``gamma``.  For the default sphere
        geometry this is exactly the legacy FRAGMENT-MNP result.  For fibres
        it uses cylindrical surface area and volume, including endcaps when
        requested in ``particle_geometry``.
        """
        psd = np.asarray(psd, dtype=float)
        ratio = self.geometry.surface_area_to_volume(psd)
        med = float(np.median(ratio))
        if med <= 0.0:
            return np.ones_like(ratio)
        return (ratio / med) ** gamma

    @staticmethod
    def _validate_inputs(config: dict, data: dict) -> Tuple[dict, dict]:
        """
        Validate the config and data dicts passed to the model

        Parameters
        ----------
        config : dict
            Model config dict
        data : dict
            Model input data dict

        Returns
        -------
        Tuple[dict, dict]
             Config and data dicts validated and filled with defaults
        """
        # Try and validate the config
        try:
            # Returns the config dict with defaults filled
            config = validation.validate_config(config)
        except SchemaError as err:
            raise SchemaError('Model config did not pass validation!') from err
        # Try and validate data
        try:
            # Returns the data dict with defaults filled
            data = validation.validate_data(data, config)
        except SchemaError as err:
            raise SchemaError('Input data did not pass validation!') from err
        # Return the config and data with filled defaults
        return config, data
    
    #############################
    # Additive analytical
    #############################
    
    @staticmethod
    def _analytical_additive_remaining_fraction(
        t: float,
        r: float,
        D_p: float,
        D_w: float,
        K_pw: float,
        n_terms: int = 50
    ) -> float:
        """
        Return the fraction of additive remaining in a spherical particle:
            F(t) = M(t) / M0

        This implements the "selection" approach based on a Biot number (Bi) regime:

        - Bi <= 1:     external transfer dominated -> simple approximation
        - 1 < Bi <100: intermediate -> smooth bridge approximation
        - Bi >= 100:   internal diffusion dominated -> classical infinite series
                       solution for radial diffusion in a sphere

        Parameters
        ----------
        t : float
            Time since (re-)initialisation of a uniform additive profile [s].
            NOTE: In our operator-splitting implementation we apply this over a
            timestep dt, which implicitly assumes the particle's additive profile
            is "effectively reset" to uniform each step (or that the analytical
            expression is used as a Markov-like step operator). This is the
            minimal-change way to integrate analytical release.
        r : float
            Particle radius [m]
        D_p : float
            Additive diffusivity in polymer [m^2/s]
        D_w : float
            Additive diffusivity in water [m^2/s]
        K_pw : float
            Polymer-water partition coefficient [-]
        n_terms : int
            Number of terms used for truncating the infinite series (Bi>=100).
            50 is usually plenty for numerical stability and accuracy.

        Returns
        -------
        float
            Remaining mass fraction F in [0, 1]
        """
        # Defensive programming: keep this stable and bounded.
        if t <= 0.0:
            return 1.0
        if r <= 0.0:
            # Zero radius is non-physical; treat as instantaneous release
            return 0.0
        if D_p <= 0.0 or D_w <= 0.0 or K_pw <= 0.0:
            raise ValueError("D_p, D_w, and K_pw must be > 0 for analytical release model.")

        # Fourier number: internal diffusion timescale for a sphere
        Fo = D_p * t / (r * r)

        # definitions:
        # Km = (K_pw * D_w) / r  [m/s]
        # Bi = Km * r / D_p = (K_pw * D_w) / D_p  [-]
        #
        # Note: Bi becomes independent of radius once Km is expressed this way.
        Bi = (K_pw * D_w) / D_p

        # -------------------------
        # Case I: low Bi (<= 1)
        # -------------------------
        # Simple approximation:
        #   F = exp(-Bi * Fo)
        if Bi <= 1.0:
            F = float(np.exp(-Bi * Fo))
            return float(np.clip(F, 0.0, 1.0))

        # -------------------------
        # Case II: intermediate Bi
        # -------------------------
        # A smooth, semi-empirical bridge is often used between limiting cases.
        # Here we use a stable exponential bridge:
        #   F = exp(-Fo * (Bi/(1+Bi)))
        # This has correct limiting behavior:
        #   Bi->0 => exp(-Bi*Fo)
        #   Bi->inf => exp(-Fo) (a reasonable bridge trend)
        if Bi < 100.0:
            F = float(np.exp(-Fo * (Bi / (1.0 + Bi))))
            return float(np.clip(F, 0.0, 1.0))

        # -------------------------
        # Case III: high Bi (>=100)
        # -------------------------
        # Classical solution for diffusion in a sphere (Dirichlet boundary):
        #   F = sum_{n=1}^\infty (6/(n^2*pi^2)) * exp(-n^2*pi^2*Fo)
        #
        # We truncate safely. This series is strictly positive and decreasing.
        n = np.arange(1, n_terms + 1, dtype=float)
        lam2 = (n * np.pi) ** 2
        terms = (6.0 / lam2) * np.exp(-lam2 * Fo)
        F = float(np.sum(terms))
        return float(np.clip(F, 0.0, 1.0))
    
    @staticmethod
    @lru_cache(maxsize=512)
    def _cylindrical_robin_eigenvalues(Bi: float, n_terms: int) -> tuple[float, ...]:
        """
        Eigenvalues for radial diffusion in a cylinder with a Robin boundary.

        Roots satisfy

            alpha * J1(alpha) = Bi * J0(alpha)

        where J0 and J1 are Bessel functions of the first kind.
        """
        Bi = float(Bi)
        n_terms = int(n_terms)
        if Bi <= 0.0:
            return tuple()
        if n_terms < 1:
            raise ValueError("n_terms must be >= 1.")

        z0 = jn_zeros(0, n_terms)

        # At very large Bi the Robin roots are indistinguishable from J0
        # zeros at double precision. Returning the limiting roots is both
        # numerically robust and physically exact for this regime.
        if Bi >= 1e8:
            return tuple(float(x) for x in z0)

        z1 = jn_zeros(1, max(n_terms - 1, 1))
        roots = []

        def eigen_eq(x):
            return x * j1(x) - Bi * j0(x)

        # First root lies between zero and the first J0 zero. Keep the upper
        # endpoint slightly inside the zero so roundoff in J0 is not amplified
        # by Bi.
        lo = np.finfo(float).eps
        hi = float(z0[0]) * (1.0 - 1e-12)
        roots.append(float(brentq(eigen_eq, lo, hi, xtol=1e-13, rtol=1e-12)))

        # Root n>=2 lies between J1 zero n-1 and J0 zero n.
        for n in range(1, n_terms):
            lo = float(z1[n - 1]) * (1.0 + 1e-12)
            hi = float(z0[n]) * (1.0 - 1e-12)
            roots.append(float(brentq(eigen_eq, lo, hi, xtol=1e-13, rtol=1e-12)))

        return tuple(roots)

    @staticmethod
    def _analytical_cylindrical_remaining_fraction(
        t: float,
        r: float,
        D_p: float,
        D_w: float,
        K_pw: float,
        n_terms: int = 50,
    ) -> float:
        r"""
        Additive remaining fraction for radial diffusion from a long cylinder.

        This is the fibre analytical release kernel.  The fibre is represented
        as an infinitely long circular cylinder for release, so ``r`` is the
        fibre radius, not half the fibre length.  Axial/end diffusion is
        intentionally neglected in this first fibre implementation.

        Governing equation
        ------------------

        .. math::
            \frac{\partial C}{\partial t}
            = D_p\frac{1}{r}\frac{\partial}{\partial r}
              \left(r\frac{\partial C}{\partial r}\right)

        with symmetry at the centre and a Robin surface boundary.  With the
        same FRAGMENT-MNP external-transfer convention used by the spherical
        release model,

        .. math::
            Fo = \frac{D_p t}{R^2},\qquad
            Bi = \frac{K_{pw}D_w}{D_p}.

        The volume-averaged remaining fraction is

        .. math::
            F(t)=\sum_n A_n\exp(-\alpha_n^2 Fo),

        where

        .. math::
            \alpha_n J_1(\alpha_n)=Bi J_0(\alpha_n)

        and

        .. math::
            A_n=\frac{4}{\alpha_n^2\left[1+(\alpha_n/Bi)^2\right]}.

        In the perfect-sink/high-Bi limit this reduces to the standard
        cylindrical Bessel series with J0 zeros and coefficients 4/alpha^2.
        """
        if t <= 0.0:
            return 1.0
        if r <= 0.0:
            return 0.0
        if D_p <= 0.0 or D_w <= 0.0 or K_pw <= 0.0:
            raise ValueError("D_p, D_w, and K_pw must be > 0 for analytical release model.")
        if n_terms < 1:
            raise ValueError("n_terms must be >= 1.")

        Fo = D_p * float(t) / (float(r) ** 2)
        Bi = (K_pw * D_w) / D_p

        # Very small Bi is accurately represented by the cylindrical lumped
        # capacitance limit, avoiding a poorly conditioned first root near zero.
        if Bi < 1e-8:
            return float(np.clip(np.exp(-2.0 * Bi * Fo), 0.0, 1.0))

        roots = np.asarray(
            FragmentMNP._cylindrical_robin_eigenvalues(float(Bi), int(n_terms)),
            dtype=float,
        )
        coeff = 4.0 / (roots**2 * (1.0 + (roots / Bi) ** 2))
        F = float(np.sum(coeff * np.exp(-(roots**2) * Fo)))
        return float(np.clip(F, 0.0, 1.0))

    #############################
    # Additive numerical 
    #############################
    @staticmethod
    def _numerical_additive_remaining_fraction(
        t: float,
        r: float,
        D_p: float,
        K_pw: float,
        k_m: float,
        n_r: int = 60,
        n_substeps: int = 20,
        theta: float = 1.0
    ) -> float:
        """
        Numerical solution for additive remaining fraction in a sphere:
            F(t) = M(t) / M0

        We solve polymer-phase diffusion in spherical coordinates using a
        *finite-volume* discretisation (mass-conservative), and a Robin
        (mass-transfer-limited) boundary at the polymer-water interface.

        Governing equation (polymer phase):
            ∂C/∂t = D_p * (1/r^2) ∂/∂r ( r^2 ∂C/∂r )

        Center symmetry:
            ∂C/∂r = 0 at r = 0

        Surface Robin BC (sink water, C_w ≈ 0):
            -D_p ∂C/∂r |_{r=R} = (k_m / K_pw) * C_s

        where:
            k_m  : mass transfer coefficient in water [m/s]
            K_pw : polymer-water partition coefficient [-]
            C_s  : polymer concentration at the surface

        Numerical notes / design choices
        -------------------------------
        - Finite-volume cells in radius ensure mass conservation.
        - We use a theta-scheme in time:
              (I - θ Δt L) C^{n+1} = (I + (1-θ) Δt L) C^n
          Default θ=1 (Backward Euler) is unconditionally stable.
        - Because the fragmentation model can have relatively large dt,
          we add n_substeps internal steps to improve accuracy while keeping
          stability. This makes the solver “more truly numerical” and less
          sensitive to dt.

        Initial condition for each call:
            C(r,0) = 1 (uniform)
        This matches your operator-splitting assumption (each timestep is
        treated as a new “release operator” acting on a mixed particle).

        Returns
        -------
        float
            Remaining mass fraction F in [0,1]
        """
        # -------------------------
        # Defensive checks
        # -------------------------
        if t <= 0.0:
            return 1.0
        if r <= 0.0:
            return 0.0
        if D_p <= 0.0:
            raise ValueError("D_p must be > 0.")
        if K_pw <= 0.0:
            raise ValueError("K_pw must be > 0.")
        if k_m < 0.0:
            raise ValueError("k_m must be >= 0.")
        if n_r < 5:
            raise ValueError("n_r should be >= 5 for a meaningful radial grid.")
        if n_substeps < 1:
            raise ValueError("n_substeps must be >= 1.")
        if not (0.0 <= theta <= 1.0):
            raise ValueError("theta must be in [0, 1].")
        
        # If particle is extremely small, diffusion timescale is tiny.
        # Treat as instantaneous release governed by boundary control.
        # This avoids ill-conditioned grids when r is ~nanometers.
        if r < 1e-8:
            # If there is any coupling to water, release ~fully in this step.
            # If k_m==0, no release.
            return 1.0 if k_m == 0.0 else 0.0

        # -------------------------
        # Radial finite-volume grid
        # -------------------------
        # Faces at: 0, dr, 2dr, ..., R
        # Cell centers at: dr/2, 3dr/2, ..., R - dr/2
        dr = r / n_r
        r_faces = np.arange(n_r + 1, dtype=float) * dr
        r_centers = (np.arange(n_r, dtype=float) + 0.5) * dr

        # In spherical FV, geometric factors (4π cancels in mass fractions):
        # - cell "volume weights" proportional to ∫ r^2 dr over the cell
        # - face "area weights" proportional to r^2 at the face
        V = (r_faces[1:]**3 - r_faces[:-1]**3) / 3.0   # proportional volume weights
        A = r_faces**2                                 # proportional face areas

        # Start uniform profile: C = 1 everywhere.
        C = np.ones(n_r, dtype=float)

        # Total time stepping (substeps for accuracy)
        dt_total = float(t)
        dt = dt_total / float(n_substeps)

        # -------------------------
        # Build diffusion operator L as a tridiagonal (FV form)
        # -------------------------
        # FV update:
        #   dC_i/dt = (D/V_i) * [ A_R (C_{i+1}-C_i)/dr - A_L (C_i - C_{i-1})/dr ] / dr
        #          = (D/(V_i*dr^2)) * [ A_R C_{i+1} - (A_R + A_L) C_i + A_L C_{i-1} ]
        #
        # So L has:
        #   L_lower[i] =  D * A_L / (V_i*dr^2)
        #   L_upper[i] =  D * A_R / (V_i*dr^2)
        #   L_diag[i]  = -(L_lower[i] + L_upper[i])  (plus boundary sink term)
        #
        # Center symmetry is naturally handled because A_L at r=0 is zero.

        beta = D_p / dr
        L_lower = np.zeros(n_r, dtype=float)
        L_diag  = np.zeros(n_r, dtype=float)
        L_upper = np.zeros(n_r, dtype=float)

        for i in range(n_r):
            Vi = V[i]
            A_L = A[i]  # face area at r=i*dr

            # Interior right face area:
            # For the LAST cell, the right face is the boundary face.
            # We must NOT treat it as an interior diffusive connection to a non-existent cell.
            if i < n_r - 1:
                A_R = A[i + 1]
            else:
                A_R = 0.0

            cL = beta * (A_L / Vi)
            cR = beta * (A_R / Vi)

            if i > 0:
                L_lower[i] = cL
            if i < n_r - 1:
                L_upper[i] = cR

            L_diag[i] = -(cL + cR)

        # -------------------------
        # Robin BC at r=R as an *effective sink* on the last cell
        # -------------------------
        # Approximate gradient between last cell center and surface:
        #   -D_p ∂C/∂r |_{r=R} = k_m * C_s
        #
        # Robin:
        #   -D (C_s - C_last)/(dr/2) = (k_m) C_s
        #
        # Solve for C_s in terms of C_last:
        #   C_s = C_last / (1 + (k_m*(dr/(2D)))
        #
        # Flux to water:
        #   J = (k_m) * C_s
        #
        # FV sink term in last cell:
        #   dC_last/dt includes -(A_face/V_last) * J
        #
        # So define an effective boundary "loss velocity":
        #   k_eff = (k_m) / (1 + (k_m)*(dr/(2D)))
        #
        # Then sink = (A_face/V_last) * k_eff
        #
        
        if k_m > 0.0:
            denom = 1.0 + k_m * (dr / (2.0 * D_p))
            k_eff = k_m / denom  # [m/s]
        else:
            k_eff = 0.0

        sink = (A[-1] / V[-1]) * k_eff   # [1/s] in the last cell equation
        L_diag[-1] -= sink               # more negative diagonal => more loss

        # -------------------------
        # Helper: tridiagonal matvec y = L x
        # -------------------------
        def L_dot(x: np.ndarray) -> np.ndarray:
            y = L_diag * x
            y[1:] += L_lower[1:] * x[:-1]
            y[:-1] += L_upper[:-1] * x[1:]
            return y

        # -------------------------
        # Time stepping (theta scheme)
        # -------------------------
        # (I - θ dt L) C^{n+1} = (I + (1-θ) dt L) C^n
        #
        # Build constant tridiagonal system matrix for each substep:
        A_lower = -theta * dt * L_lower
        A_diag = 1.0 - theta * dt * L_diag
        A_upper = -theta * dt * L_upper

        # Factor once, solve many times
        A_lower_fact, A_diag_fact, A_upper_fact = FragmentMNP._factor_tridiagonal(
            A_lower, A_diag, A_upper
        )

        for _ in range(n_substeps):
            if theta == 1.0:
                rhs = C.copy()
            else:
                rhs = C + (1.0 - theta) * dt * L_dot(C)

            C = FragmentMNP._solve_tridiagonal_factored(
                A_lower_fact, A_diag_fact, A_upper_fact, rhs
            )

            # Safety clip (numerical roundoff can create tiny negatives)
            C = np.clip(C, 0.0, None)

        # -------------------------
        # Mass remaining fraction
        # -------------------------
        M0 = float(np.sum(V * 1.0))      # initial uniform profile
        M1 = float(np.sum(V * C))        # mass-weighted remaining
        return float(np.clip(M1 / M0, 0.0, 1.0))
    
    @staticmethod
    def _numerical_cylindrical_remaining_fraction(
        t: float,
        r: float,
        D_p: float,
        K_pw: float,
        k_m: float,
        n_r: int = 60,
        n_substeps: int = 20,
        theta: float = 1.0,
    ) -> float:
        r"""
        Numerical radial diffusion from a long cylindrical fibre.

        The finite-volume and theta-scheme structure mirrors the established
        spherical numerical release solver.  Only the radial geometry factors
        change:

        - cylindrical cell volume weight: integral r dr
        - cylindrical face area weight: r

        The external boundary convention is deliberately identical to the
        existing spherical implementation: ``k_m`` is the effective surface
        loss velocity used in the Robin sink.  When ``D_w`` is supplied by the
        calling release wrapper, ``k_m = K_pw * D_w / R``.
        """
        if t <= 0.0:
            return 1.0
        if r <= 0.0:
            return 0.0
        if D_p <= 0.0:
            raise ValueError("D_p must be > 0.")
        if K_pw <= 0.0:
            raise ValueError("K_pw must be > 0.")
        if k_m < 0.0:
            raise ValueError("k_m must be >= 0.")
        if n_r < 5:
            raise ValueError("n_r should be >= 5 for a meaningful radial grid.")
        if n_substeps < 1:
            raise ValueError("n_substeps must be >= 1.")
        if not (0.0 <= theta <= 1.0):
            raise ValueError("theta must be in [0, 1].")

        if r < 1e-8:
            return 1.0 if k_m == 0.0 else 0.0

        dr = r / n_r
        r_faces = np.arange(n_r + 1, dtype=float) * dr

        # Circular-cylinder geometry; common 2*pi*L factors cancel.
        V = (r_faces[1:]**2 - r_faces[:-1]**2) / 2.0
        A = r_faces.copy()

        C = np.ones(n_r, dtype=float)
        dt = float(t) / float(n_substeps)

        beta = D_p / dr
        L_lower = np.zeros(n_r, dtype=float)
        L_diag = np.zeros(n_r, dtype=float)
        L_upper = np.zeros(n_r, dtype=float)

        for i in range(n_r):
            Vi = V[i]
            A_L = A[i]
            A_R = A[i + 1] if i < n_r - 1 else 0.0

            cL = beta * (A_L / Vi)
            cR = beta * (A_R / Vi)

            if i > 0:
                L_lower[i] = cL
            if i < n_r - 1:
                L_upper[i] = cR
            L_diag[i] = -(cL + cR)

        if k_m > 0.0:
            denom = 1.0 + k_m * (dr / (2.0 * D_p))
            k_eff = k_m / denom
        else:
            k_eff = 0.0

        sink = (A[-1] / V[-1]) * k_eff
        L_diag[-1] -= sink

        def L_dot(x: np.ndarray) -> np.ndarray:
            y = L_diag * x
            y[1:] += L_lower[1:] * x[:-1]
            y[:-1] += L_upper[:-1] * x[1:]
            return y

        A_lower = -theta * dt * L_lower
        A_diag = 1.0 - theta * dt * L_diag
        A_upper = -theta * dt * L_upper

        A_lower_fact, A_diag_fact, A_upper_fact = FragmentMNP._factor_tridiagonal(
            A_lower, A_diag, A_upper
        )

        for _ in range(n_substeps):
            rhs = C.copy() if theta == 1.0 else C + (1.0 - theta) * dt * L_dot(C)
            C = FragmentMNP._solve_tridiagonal_factored(
                A_lower_fact, A_diag_fact, A_upper_fact, rhs
            )
            C = np.clip(C, 0.0, None)

        M0 = float(np.sum(V))
        M1 = float(np.sum(V * C))
        return float(np.clip(M1 / M0, 0.0, 1.0))

    @staticmethod
    def _solve_tridiagonal(lower: np.ndarray,
                           diag: np.ndarray,
                           upper: np.ndarray,
                           rhs: np.ndarray) -> np.ndarray:
        """
        Solve a tridiagonal linear system Ax = rhs with Thomas algorithm.

        lower[i] = A[i,i-1] for i>=1 (lower[0] unused or 0)
        diag[i]  = A[i,i]
        upper[i] = A[i,i+1] for i<=n-2 (upper[n-1] unused or 0)
        """
        n = len(diag)
        a = lower.astype(float, copy=True)
        b = diag.astype(float, copy=True)
        c = upper.astype(float, copy=True)
        d = rhs.astype(float, copy=True)

        # Forward elimination
        for i in range(1, n):
            if b[i - 1] == 0.0:
                raise ZeroDivisionError("Tridiagonal solver encountered zero pivot.")
            w = a[i] / b[i - 1]
            b[i] = b[i] - w * c[i - 1]
            d[i] = d[i] - w * d[i - 1]

        # Back substitution
        x = np.zeros(n, dtype=float)
        if b[-1] == 0.0:
            raise ZeroDivisionError("Tridiagonal solver encountered zero pivot at last row.")
        x[-1] = d[-1] / b[-1]
        for i in range(n - 2, -1, -1):
            if b[i] == 0.0:
                raise ZeroDivisionError("Tridiagonal solver encountered zero pivot.")
            x[i] = (d[i] - c[i] * x[i + 1]) / b[i]

        return x
    
    @staticmethod
    def _expand_rate_to_size_vector(rate, n_size_classes: int) -> np.ndarray:
        """
        Convert a scalar or length-n_size_classes iterable into a size-vector.
        """
        if isinstance(rate, (int, float)):
            return np.full(n_size_classes, float(rate), dtype=float)

        arr = np.asarray(rate, dtype=float)
        if arr.ndim != 1 or len(arr) != n_size_classes:
            raise ValueError(
                f"Expected scalar or length-{n_size_classes} rate vector."
            )
        return arr.astype(float, copy=False)