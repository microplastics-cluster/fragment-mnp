from typing import Tuple
import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp
from scipy import interpolate
from schema import SchemaError
from . import validation
from .output import FMNPOutput
from ._errors import FMNPNumericalError


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
        self.t_eval = np.arange(0, self.n_timesteps) \
            if self.config['solver_t_eval'] == 'integer' \
            else self.config['solver_t_eval']
        # Initial concentrations
        self.initial_concs = np.array(data['initial_concs'])
        # Set the particle phys-chem properties
        self.psd = self._set_psd()
        self.fsd = self.set_fsd(self.n_size_classes,
                                self.psd,
                                self.data['fsd_beta'])
        self.theta_1 = data['theta_1']
        self.density = data['density']
        self.k_diss = self._set_k_diss(data['k_diss'],
                                       config['k_diss_scaling_method'],
                                       self.psd,
                                       self.n_size_classes,
                                       data['k_diss_gamma'])
        self.k_frag = self._set_k_frag(data['k_frag'], self.theta_1,
                                       self.data['k_frag_tau'],
                                       self.psd, self.n_timesteps)

    def run(self) -> FMNPOutput:
        r"""
        Run the model with the config and data provided at initialisation.

        Returns
        -------
        :class:`fragmentmnp.output.FMNPOutput` object containing model output
        data.

        Notes
        -----
        Internally the model numerically solves the following differential
        equation for each size class, to give a time series of mass
        concentrations `c`. `k` is the current size class, `i` are the
        daughter size classes.

        .. math::
            \frac{dc_k}{dt} = -k_{\text{frag},k} c_k +
            \Sigma_i f_{i,k} k_{\text{frag},i} c_i - k_{\text{diss},k} c_k

        Here, :math:`k_{\text{frag},k}` is the fragmentation rate of size
        class `k`, :math:`f_{i,k}` is the fraction of daughter
        fragments produced from a fragmenting particle of size `i` that are of
        size `k`, and :math:`k_{\text{diss},k}` is the dissolution rate from
        size class `k`.

        Mass concentrations are converted to particle number concentrations by
        assuming spherical particles with the density given in the input data.
        """

        # Define the initial value problem to pass to SciPy to solve.
        # This must satisfy c'(t) = f(t, c) with initial values given in data.
        def f(t, c):
            # Get the number of size classes and create results to be filled
            N = self.n_size_classes
            dcdt = np.empty(N)
            # Interpolate the time-dependent parameters to the specific
            # timestep given (which will be a float, rather than integer index)
            t_model = np.arange(self.n_timesteps)
            f = interpolate.interp1d(t_model, self.k_frag, axis=0,
                                     fill_value='extrapolate')
            k_frag = f(t)
            # Loop over the size classes and perform the calculation
            for k in np.arange(N):
                # The differential equation that is being solved
                dcdt[k] = - k_frag[k] * c[k] \
                    + np.sum(self.fsd[:, k] * k_frag * c[:N]) \
                    - self.k_diss[k] * c[k]
            # Return the solution for all of the size classes
            return dcdt

        # Numerically solve this given the initial values for c
        soln = solve_ivp(fun=f,
                         method=self.config['solver_method'],
                         t_span=(0, self.n_timesteps),
                         y0=self.initial_concs,
                         t_eval=self.t_eval,
                         rtol=self.config['solver_rtol'],
                         atol=self.config['solver_atol'],
                         max_step=self.config['solver_max_step'])
        # If we didn't find a solution, raise an error
        if not soln.success:
            raise FMNPNumericalError('Model solution could not be ' +
                                     f'found: {soln.message}')
        # Calculate the timeseries of mass concentration lost to dissolution
        # from the solution
        j_diss = self.k_diss[:, None] * soln.y
        # Use this to calculate the cumulative mass concentration
        # lost to dissolution
        c_diss = np.cumsum(j_diss, axis=1)
        # Now we have the solutions as mass concentrations, we can
        # convert to particle number concentrations
        n = self.mass_to_particle_number(soln.y)
        n_diss = self.mass_to_particle_number(c_diss)

        # Return the solution in an FMNPOutput object
        return FMNPOutput(soln.t, soln.y, n, c_diss, n_diss, soln, self.psd)

    def mass_to_particle_number(self, mass):
        """
        Convert mass (concentration) to particle number (concentration).
        Here, particles are assumed to be spherical, but this function
        can be overloaded to account for different shape particles.
        """
        n = mass / (self.density * (4.0/3.0) * np.pi
                    * (self.psd[:, None] / 2) ** 3)
        return n

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

    @staticmethod
    def _set_k_frag(k_frag: float, theta_1: float, tau: float,
                    psd: npt.NDArray[np.float64],
                    n_timesteps: int) -> npt.NDArray[np.float64]:
        r"""
        Set the fragmentation rate `k_frag` based on either the average
        `k_frag` for the median particle size bin and `theta_1` (surface
        energy empirical parameter), or directly if a distribution is
        provided. `tau` (time dependence empirical parameter) then scales
        `k_frag` to be dependent on time.

        Parameters
        ----------
        k_frag : float or iterable
            Either the average :math:`k_frag` for the median particle
            size bin, or a distribution of :math:`k_frag` values
        theta_1 : float
            The surface energy empirical parameter :math:`\theta_1`
        tau : float
            The time-dependence parameter :math:`\tau`
        psd : np.ndarray
            The particle size distribution
        n_timesteps : int
            The number of model timesteps

        Returns
        -------
        k_frag = np.ndarray (n_timesteps, n_size_classes)
            Fragmentation rate array over timesteps and size classes
        """
        # Check if k_frag is a scalar value, in which case we need
        # to calculate a distribution based on theta1
        if isinstance(k_frag, (int, float)):
            # Get the proportionality constant
            k_prop = k_frag / (np.median(psd) ** (2 * theta_1))
            # Now create the array of k_frags
            k_frag_dist = k_prop * psd ** (2 * theta_1)
            # We presume fragmentation from the smallest size class
            # can't happen, and the only loss from this size class
            # is from dissolution
            k_frag_dist[0] = 0.0
        # Else just set k_frag directly from the provided array.
        # Validation makes sure this is the correct length
        else:
            k_frag_dist = np.array(k_frag)
        # Now set the time dependence of k_frag using tau
        t = np.arange(1, n_timesteps + 1)
        k_frag_2d = k_frag_dist * t[:, np.newaxis] ** tau / np.median(t) ** tau
        return k_frag_2d

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
    def _set_k_diss(k_diss: float,
                    scaling_method: str,
                    psd: npt.NDArray[np.float64],
                    n_size_classes: int,
                    gamma: float = 1.0) -> npt.NDArray[np.float64]:
        """
        Set the dissolution rate for each of the size classes,
        based on an average dissolution rate

        Parameters
        ----------
        k_diss : float
            Either average dissolution rate across size classes, or the
            full distribution.
        scaling_method: str
            How to scale ``k_diss`` across size classes? Either `constant`
            or `surface_area`. If `constant`, the same ``k_diss`` is used
            for all size classes. If `surface_area`, ``k_diss`` is scaled
            according to particle surface area per unit volume of polymer.
        psd : np.ndarray
            The particle size distribution
        n_size_classes : int
            The number of particle size classes

        Returns
        -------
        np.ndarray
            Dissolution rate distribution

        Notes
        -----
        At the moment, we are assuming spherical particles when scaling
        by surface area. This might change.
        """
        # Check if k_diss is a scalar value, in which case we need
        # to calculate a distribution based on gamma
        if isinstance(k_diss, (int, float)):
            # What scaling method has been chosen?
            if scaling_method == 'constant':
                # Use the average k_diss across all size classes
                k_diss_dist = np.full((n_size_classes,), k_diss)
            elif scaling_method == 'surface_area':
                # Scale the dissolution according to surface area per unit
                # volume, assuming our particles are spheres
                k_diss_dist = k_diss * FragmentMNP._f_surface_area(psd, gamma)
            else:
                # We shouldn't get here, if validation has been performed!
                raise ValueError('Invalid k_diss_scaling_factor provided: ',
                                 {scaling_method})
        # Otherwise we will have been given a distribution, so use
        # that directly
        else:
            k_diss_dist = k_diss
        return k_diss_dist

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
