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
import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp
from scipy import interpolate
from schema import SchemaError
from . import validation
from .output import FMNPOutput
from ._errors import FMNPNumericalError, FMNPDistributionValueError


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
        Initialise the model
        """
        # Validate the config and data (if we are meant to)
        if validate:
            config, data = self._validate_inputs(config, data)
        # If we passed validation (or aren't validating), save attributes
        self.config = config
        self.data = data
        # Set the number of particle size classes and timesteps, and
        # the times at which to store the computed solution
        self.n_size_classes = self.config['n_size_classes']
        self.n_timesteps = self.config['n_timesteps']
        self.dt = self.config['dt']
        self.t_grid = np.arange(0.5*self.dt,
                                self.n_timesteps*self.dt + 0.5*self.dt,
                                self.dt)
        self.t_eval = self.t_grid \
            if self.config['solver_t_eval'] == 'timesteps' \
            else self.config['solver_t_eval']
        # Initial concentrations
        self.initial_concs = np.array(data['initial_concs'], dtype=float)
        self.initial_concs_diss = data['initial_concs_diss']
        # Set the particle phys-chem properties
        self.psd = self._set_psd()
        self.surface_areas = self.surface_area(self.psd)
        self.fsd = self.set_fsd(self.n_size_classes,
                                self.psd,
                                self.data['fsd_beta'])
        self.density = float(data['density'])
        # Stop Pylance complaining about k_frag, k_diss and k_min not being
        # present (or being the wrong type) by declaring them as NumPy arrays
        self.k_frag = np.empty((self.n_size_classes, self.n_timesteps))
        self.k_diss = np.empty((self.n_size_classes, self.n_timesteps))
        self.k_min = np.empty((self.n_timesteps,))
        # Calculate the rate constant distributions
        for k in ['k_frag', 'k_diss', 'k_min']:
            k_dist = self.set_rate_constant(data[k], k)
            setattr(self, k, k_dist)

        # ------------------------------------------------------------
        # OPTIONAL: additive coupling inputs (collaboration feature)
        # ------------------------------------------------------------
        init_A = data.get('initial_additive_concs', None)
        self.initial_additive_concs = None if init_A is None else np.array(init_A, dtype=float)

        # Example shape:
        #   {'model': 'analytical', 'params': {'Dp':..., 'Dw':..., 'Kpw':..., 'Km':...}}
        self.additive_release = data.get('additive_release', None)

    def run(self) -> FMNPOutput:
        r"""
        Run the model with the config and data provided at initialisation.

        Returns
        -------
        :class:`fragmentmnp.output.FMNPOutput` object containing model output

        Notes
        -----
        The model numerically solves the following set of differential
        equations to give a time series of mass concentrations of particles
        `c`, dissolved polymer `c_diss` and mineralised polymer `c_min`.
        `k` is the current size class, `i` are the daughter size classes.

        .. math::
            \frac{dc_k}{dt} = -k_{\text{frag},k} c_k +
            \sum_i f_{i,k} k_{\text{frag},i} c_i - k_{\text{diss},k} c_k

        .. math::
            \frac{dc_\text{diss}}{dt} = \sum_k k_{\text{diss},k} c_k -
            k_\text{min} c_\text{diss}

        .. math::
            \frac{dc_\text{min}}{dt} = k_\text{min} c_\text{diss}

        Here, :math:`k_{\text{frag},k}` is the fragmentation rate of size
        class `k`, :math:`f_{i,k}` is the fraction of daughter
        fragments produced from a fragmenting particle of size `i` that are of
        size `k`, :math:`k_{\text{diss},k}` is the dissolution rate from
        size class `k` and :math:`k_\text{min}` is the mineralisation rate
        from the dissolved pool.

        Mass concentrations are converted to particle number concentrations by
        assuming spherical particles with the density given in the input data.
        """
        def f(t, c):
            """
            The initial value problem for SciPy to solve. This must satisfy
            c'(t) = f(t, c) with initial values given in data, with N+2
            solutions:
              c[:N]   = mass concentration in each of the N size classes
              c[N]    = total dissolved mass concentration
              c[N+1]  = mineralised mass concentration
            """
            # Get the number of size classes
            N = self.n_size_classes
            # Unpack the solutions
            c_particles = c[:N]
            c_dissolved = c[N]
            c_min = c[N + 1]   # Not used as a source in any ODE (only grows)
            # Interpolate the time-dependent parameters to the specific
            # timestep given (which will be a float, rather than integer index)
            f_frag = interpolate.interp1d(self.t_grid, self.k_frag, axis=1,
                                          fill_value='extrapolate')
            f_diss = interpolate.interp1d(self.t_grid, self.k_diss, axis=1,
                                          fill_value='extrapolate')
            f_min = interpolate.interp1d(self.t_grid, self.k_min, axis=0,
                                         fill_value='extrapolate')
            k_frag = f_frag(t)
            k_diss = f_diss(t)
            k_min = f_min(t)
            # Build array of d/dt
            dcdt = np.empty(N+2)
            # Particle ODEs, for each size class:
            #   Gains: fragmentation from bigger size classes
            #   Loss:  dissolution and fragmentation to smaller size classes
            for k in range(N):
                dcdt[k] = (
                    - k_frag[k] * c_particles[k]
                    + np.sum(self.fsd[:, k] * k_frag * c_particles[:N])
                    - k_diss[k] * c_particles[k]
                )
            # Dissolved mass ODE:
            #   Gains: sum of kdiss[k] * c_particles[k]
            #   Loss:  k_min * c_dissolved
            dcdt_dissolved = np.sum(k_diss * c_particles) - k_min * c_dissolved
            # Assign to last entry
            dcdt[N] = dcdt_dissolved
            # Mineralized ODE:
            #   Gains: k_min_avg * c_dissolved
            #   No loss, so it's purely accumulative
            dcdt_min = k_min * c_dissolved
            dcdt[N + 1] = dcdt_min
            # Final differential to return
            return dcdt

        # Build the new N+2 initial conditions: N particulate states,
        # +1 dissolved, +1 mineralized
        y0 = np.concatenate([
           self.initial_concs,          # microplastic mass per size class
           [self.initial_concs_diss],   # dissolved mass
           [0.0]                        # mineralized (initially zero)
        ])
        # Solve the ODE
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
        # Extract solution
        c_part_sol = soln.y[:self.n_size_classes, :]
        c_diss_sol = soln.y[self.n_size_classes, :]
        c_min_sol = soln.y[self.n_size_classes + 1, :]
        # Convert microparticle mass to particle number
        n_part_sol = self.mass_to_particle_number(c_part_sol)

        # ------------------------------------------------------------
        # OPTIONAL: compute additive time series (post-solve)
        # ------------------------------------------------------------
        A_part = None
        A_aq = None

        if (self.initial_additive_concs is not None) and (self.additive_release is not None):
            A_part, A_aq = self._simulate_additives_postsolve(soln.t, c_part_sol)

        # Build the FMNPOutput object from the solution
        return FMNPOutput(
           t=soln.t,
           c=c_part_sol,
           n=n_part_sol,
           c_diss=c_diss_sol,
           c_min=c_min_sol,
           soln=soln,
           psd=self.psd,
           A_part=A_part,
           A_aq=A_aq
        )
    
    # ---------------------------------------------------------------------
    # Additive coupling implementation (post-solve bookkeeping)
    # ---------------------------------------------------------------------

    def _simulate_additives_postsolve(self,
                                      t: npt.NDArray[np.float64],
                                      c_part_sol: npt.NDArray[np.float64]):
        """
        Compute additive trajectories using the already-solved polymer output.

        This method does NOT change polymer dynamics.

        We do:
          1) transport additive with polymer fragmentation
          2) release additive to water using an analytical fraction

        Inputs
        ------
        t : array (T,)
            solver output times
        c_part_sol : array (N, T)
            polymer mass per size class over time

        Returns
        -------
        A_part : array (N, T)
            additive mass in particulate bins
        A_aq : array (T,)
            additive mass in aqueous pool
        """
        N, T = c_part_sol.shape
        A_part = np.zeros((N, T), dtype=float)
        A_aq = np.zeros((T,), dtype=float)

        A_part[:, 0] = self.initial_additive_concs
        A_aq[0] = 0.0

        # Interpolate k_frag/k_diss from internal t_grid onto solver times.
        f_frag = interpolate.interp1d(self.t_grid, self.k_frag, axis=1,
                                      fill_value='extrapolate')
        # k_diss not used yet, but kept for future coupling
        # f_diss = interpolate.interp1d(self.t_grid, self.k_diss, axis=1,
        #                               fill_value='extrapolate')

        radii = self.psd / 2.0

        release_model = str(self.additive_release.get('model', 'analytical')).lower()
        release_params = self.additive_release.get('params', {})

        # Optional: fail fast if the release model is unknown
        allowed = {"analytical", "numerical"}
        if release_model not in allowed:
            raise ValueError(
                f"Unknown additive_release model '{release_model}'. "
                f"Allowed values are {sorted(allowed)}."
            )

        # Optional: if numerical is selected, ensure required params exist
        if release_model == "numerical":
            required = ["D_p", "K_pw"]
            missing = [k for k in required if k not in release_params]
            if missing:
                raise ValueError(
                    f"additive_release.params missing required keys for numerical model: {missing}. "
                    f"Provided keys: {sorted(release_params.keys())}"
                )

        eps = 1e-30  # avoid division by zero when polymer mass is ~0

        for ti in range(T - 1):
            dt_i = float(t[ti + 1] - t[ti])

            # polymer mass in each bin at current time
            c = c_part_sol[:, ti]

            # fragmentation rate at current time
            k_frag = f_frag(t[ti])

            # -------------------------
            # (1) Additive transport with fragmentation
            # -------------------------
            # Polymer mass lost by fragmentation in dt_i:
            #   L_frag[i] = k_frag[i] * c[i] * dt_i
            L_frag = k_frag * c * dt_i

            # Polymer mass transferred from parent i to daughter k:
            #   G[i,k] = fsd[i,k] * L_frag[i]
            G = (self.fsd.T * L_frag).T  # (N,N)

            # Assume additive is uniformly mixed within each size class:
            # additive per polymer mass:
            conc_A = A_part[:, ti] / np.maximum(c, eps)

            # Additive transferred i -> k:
            A_to_k = (G.T * conc_A).T  # (N,N)
            A_gain = A_to_k.sum(axis=0)
            A_loss = A_to_k.sum(axis=1)

            # particulate additive after fragmentation redistribution
            A_mid = A_part[:, ti] - A_loss + A_gain

            # -------------------------
            # (2) Additive release to water (analytical or numerical)
            # -------------------------
            if release_model == 'analytical':
                rel = np.array([
                    self._analytical_additive_release_fraction(
                        radius_m=float(radii[i]),
                        dt=dt_i,
                        params=release_params
                    )
                    for i in range(N)
                ], dtype=float)

            elif release_model == 'numerical':
                rel = np.array([
                    self._numerical_additive_release_fraction(
                        radius_m=float(radii[i]),
                        dt=dt_i,
                        params=release_params
                    )
                    for i in range(N)
                ], dtype=float)

            else:
                # Unknown release model -> no release (safe default)
                rel = np.zeros((N,), dtype=float)

            rel = np.clip(rel, 0.0, 1.0)

            dA_rel = rel * A_mid
            A_part[:, ti + 1] = A_mid - dA_rel
            A_aq[ti + 1] = A_aq[ti] + float(dA_rel.sum())

        return A_part, A_aq

    def _analytical_additive_release_fraction(self, radius_m: float, dt: float, params: dict) -> float:
        """
        Return the fraction of additive released from a particle of radius_m
        over timestep dt.

        This is what your post-solve bookkeeping loop calls.

        Required params in params dict:
            - "D_p" : polymer diffusivity [m^2/s]
            - "D_w" : water diffusivity [m^2/s]
            - "K_pw": polymer-water partition coefficient [-]

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

        F_remain = self._analytical_additive_remaining_fraction(
            t=float(dt),
            r=float(radius_m),
            D_p=D_p,
            D_w=D_w,
            K_pw=K_pw,
            n_terms=n_terms
        )

        # Release fraction = 1 - remaining fraction (bounded)
        rel = 1.0 - F_remain
        return float(np.clip(rel, 0.0, 1.0))
    
    
    def _numerical_additive_release_fraction(self, radius_m: float, dt: float, params: dict) -> float:
        """
        Release fraction computed by numerically solving diffusion in a sphere.

        Required:
            D_p, K_pw
        Either:
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
            # TODO move this check into validation
            raise ValueError('particle_size_classes or particle_size_range ' +
                             'must be present in the model config, but ' +
                             'neither were found.')
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

    @staticmethod
    def surface_area(psd: npt.NDArray[np.float64]) -> \
            npt.NDArray[np.float64]:
        """
        Return the surface area of the particles, presuming they are
        spheres. This function can be overloaded to account for different
        shaped particles.
        """
        return 4.0 * np.pi * (psd / 2.0) ** 2

    @staticmethod
    def volume(psd: npt.NDArray[np.float64]) -> \
            npt.NDArray[np.float64]:
        """
        Return the volume of the particles, presuming they are spheres.
        This function can be overloaded to account for different shaped
        particles.
        """
        return (4.0/3.0) * np.pi * (psd / 2.0) ** 3

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

    @staticmethod
    def _f_surface_area(psd: npt.NDArray[np.float64],
                        gamma: float = 1.0) -> npt.NDArray[np.float64]:
        r"""
        Calculate the scaling factor for surface area, which is defined
        as the ratio of the surface area to volume ratio of the polymer
        for each size class to the median size class, such that ``f`` is
        1 for the median size class, larger for the smaller size classes
        (because there are more particles per unit volume), and smaller
        for larger size classes. An empirical parameter ``gamma`` linearly
        scales the factor by :math:`f^\gamma`.

        Parameters
        ----------
        psd : np.ndarray
            The particle size distribution
        gamma: float
            Empirical scaling factor that scales ``f`` as :math:`s^\gamma`,
            where ``s`` is the surface area to volume ratio of each size
            class. Therefore, if ``gamma`` is 1, then ``k_diss`` scales
            directly with ``s``.

        Returns
        -------
        np.ndarray
            Surface area scaling factor

        Notes
        -----
        By assuming spherical particles, calculating their volumes and
        surface areas and simplifying the algebra, ``f`` can be defined as

        .. math::
            f_\text{s} = \left(\frac{s}{\hat{s}}\right)^\gamma

        where :math:`s` is the surface area to volume ratio, and
        :math:`\hat{s}` is the median of :math:`s`:

        .. math::
            s = \frac{4 \pi r_\text{max}}{\textbf{r}}

        Here, :math:`r_\text{max}` is the radius of the largest particle size
        class, and :math:`\textbf{r}` is an array of the particle size class
        radii.
        """
        # Calculate ratio of surface area to the largest volume and scale
        # to the median (so f = 1 for the median size class)
        surface_area_volume_ratio = (4 * np.pi * (psd.max() / 2) ** 3) / psd
        f = surface_area_volume_ratio / np.median(surface_area_volume_ratio)
        return f ** gamma

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
        #   ∂C/∂r|R ≈ (C_s - C_last) / (dr/2)
        #
        # Robin:
        #   -D (C_s - C_last)/(dr/2) = (k_m/K_pw) C_s
        #
        # Solve for C_s in terms of C_last:
        #   C_s = C_last / (1 + (k_m/K_pw)*(dr/(2D)))
        #
        # Flux to water:
        #   J = (k_m/K_pw) * C_s
        #
        # FV sink term in last cell:
        #   dC_last/dt includes -(A_face/V_last) * J
        #
        # So define an effective boundary "loss velocity":
        #   k_eff = (k_m/K_pw) / (1 + (k_m/K_pw)*(dr/(2D)))
        #
        # Then sink = (A_face/V_last) * k_eff
        #
        
        if k_m > 0.0:
            km_over_K = k_m / K_pw
            denom = 1.0 + km_over_K * (dr / (2.0 * D_p))
            k_eff = km_over_K / denom  # [m/s]
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
        A_diag  = 1.0 - theta * dt * L_diag
        A_upper = -theta * dt * L_upper

        for _ in range(n_substeps):
            if theta == 1.0:
                rhs = C.copy()
            else:
                rhs = C + (1.0 - theta) * dt * L_dot(C)

            C = FragmentMNP._solve_tridiagonal(A_lower, A_diag, A_upper, rhs)

            # Safety clip (numerical roundoff can create tiny negatives)
            C = np.clip(C, 0.0, None)

        # -------------------------
        # Mass remaining fraction
        # -------------------------
        M0 = float(np.sum(V * 1.0))      # initial uniform profile
        M1 = float(np.sum(V * C))        # mass-weighted remaining
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

