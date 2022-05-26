from typing import NamedTuple, Tuple
import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp
from schema import SchemaError
from . import validation
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

    __slots__ = ['config', 'data', 'n_size_classes', 'psd', 'fsd',
                 'n_timesteps', 'density', 'k_frag', 'theta_1', 'k_diss',
                 'initial_concs']

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
        # Set the number of particle size classes and number of timesteps
        self.n_size_classes = self.config['n_size_classes']
        self.n_timesteps = self.config['n_timesteps']
        # Initial concentrations
        self.initial_concs = np.array(data['initial_concs'])
        # Set the particle phys-chem properties
        self.psd = self._set_psd()
        self.fsd = self._set_fsd(self.n_size_classes)
        self.theta_1 = data['theta_1']
        self.density = data['density']
        self.k_diss = self._set_k_diss(data['k_diss'],
                                       config['k_diss_scaling_method'],
                                       self.psd,
                                       self.n_size_classes,
                                       data['k_diss_gamma'])
        self.k_frag = self._set_k_frag(data['k_frag'], self.theta_1,
                                       self.psd,
                                       self.n_size_classes)

    def run(self) -> NamedTuple:
        r"""
        Run the model with the config and data provided at initialisation.

        Returns
        -------
        Named tuple with the following fields defined:
        t : np.ndarray, shape (n_timesteps,)
            Time series over which the model was run
        n : np.ndarray, shape (n_size_classes, n_timesteps)
            Particle number concentrations for each size class over
            this time series
        n_diss : np.ndarray, shape (n_size_classes, n_timesteps)
            Particle number concentrations lost from size classes due
            to dissolution
        c_diss : np.ndarray, shape (n_size_classes, n_timesteps)
            Mass concentrations of dissolved organics

        Notes
        -----
        Internally the model numerically solves the following differential
        equation for each size class, to give a time series of particle number
        concentrations `n`. `k` is the current size class, `i` are the
        daughter size classes.

        .. math::
            \frac{dn_k}{dt} = -k_{\text{frag},k} n_k +
            \Sigma_i f_{i,k} k_{\text{frag},i} n_i - k_{\text{diss},k} n_k

        Here, :math:`k_{\text{frag},k}` is the fragmentation rate of size class
        `k`, :math:`f_{i,k}` is the fraction of daughter fragments produced
        from a fragmenting particle of size `i` that are of size `k`, and
        :math:`k_{\text{diss},k}` is the dissolution rate from size class `k`.
        """

        # Define the initial value problem to pass to SciPy to solve.
        # This must satisfy n'(t) = f(t, n) with initial values given in data.
        # In our case, n is an array of equations covering the differential
        # equation for each size class, n[:n_size_classes], and particle
        # numbers lost from each size class due to dissolution in
        # n[n_size_classes:]. The actual IVP that is solved is the first
        # n_size_classes equations in n, and the dissolution elements after
        # this are only used to return dissolution pools to be output from
        # the model
        def f(t, n):
            # Get the number of size classes and create results to be filled
            N = self.n_size_classes
            dndt = np.empty(N)
            n_diss = np.empty(N)
            # Loop over the size classes and perform the calculation
            for k in np.arange(N):
                # Particle number concentration lost from this size class due
                # to dissolution
                n_diss[k] = self.k_diss[k] * n[k]
                # The differential equation that is being solved
                dndt[k] = - self.k_frag[k] * n[k] \
                    + np.sum(self.fsd[:, k] * self.k_frag * n[:N]) \
                    - n_diss[k]
            # Return the solution for all of the size classes, including
            # particle number concentration in each size class, and particle
            # number concentrations lost due to dissolution
            return [*dndt, *n_diss]

        # Numerically solve this given the initial values for n
        soln = solve_ivp(fun=f,
                         t_span=(0, self.n_timesteps),
                         y0=[*self.initial_concs,
                             *np.zeros(self.n_size_classes)],
                         t_eval=np.arange(0, self.n_timesteps))
        # If we didn't find a solution, raise an error
        if not soln.success:
            raise FMNPNumericalError('Model solution could not be ' +
                                     f'found: {soln.message}')
        # Extract the particle number concentrations and dissolution losses
        n = soln.y[:self.n_size_classes]
        n_diss = soln.y[self.n_size_classes:]
        # Now we need to convert dissolution loss from a number concentration
        # to a mass concentration using the polymer density and assuming
        # spherical particles
        c_diss = n_diss * self.density * (4/3) * np.pi \
                 * (self.psd[:, None] / 2) ** 3

        # Create a named tuple to return the results in
        FMNPOutput = NamedTuple('FMNPOutput', [('t', npt.NDArray),
                                               ('n', npt.NDArray),
                                               ('c_diss', npt.NDArray),
                                               ('n_diss', npt.NDArray)])
        # Return the solution in this named tuple
        return FMNPOutput(soln.t, n, c_diss, n_diss)

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
    def _set_k_frag(k_frag_av: float, theta_1: float,
                    psd: npt.NDArray[np.float64],
                    n_size_classes: int) -> npt.NDArray[np.float64]:
        r"""
        Set the fragmentation rate :math:`k_frag` based on
        the average :math:`k_frag` for the median particle size bin,
        and :math:`\theta_1` (surface energy empirical parameter)

        Parameters
        ----------
        k_frag_av : float
            The average :math:`k_frag` for the median particle size bin
        theta_1 : float
            The surface energy empirical parameter :math:`\theta_1`
        psd : np.ndarray
            The particle size distribution
        n_size_classes : int
            The number of particle size classes

        Returns
        -------
        k_frag = np.ndarray
            Fragmentation rate array over particle size classes
        """
        # Get the proportionality constant
        k_prop = k_frag_av / (np.median(psd) ** (2 * theta_1))
        # Now create the array of k_frags
        k_frag = k_prop * psd ** (2 * theta_1)
        # We presume fragmentation from the smallest size class
        # can't happen, and the only loss from this size class
        # is from dissolution
        k_frag[0] = 0.0
        return k_frag

    @staticmethod
    def _set_fsd(n_size_classes: int) -> npt.NDArray[np.float64]:
        """
        Set the fragment size distribution matrix, assuming that
        fragmentation events result in even split between size classes
        of daughter particles

        Parameters
        ----------
        n_size_classes : int
            Number of particle size classes

        Returns
        -------
        np.ndarray
            Matrix of fragment size distributions for all size classes
        """
        # Start with a zero-filled array of shape (N,N)
        fsd = np.zeros((n_size_classes, n_size_classes))
        # Fill with an equal split between daughter size classes
        for i in np.arange(n_size_classes):
            fsd[i, :] = 1 / i if i != 0 else 0
        # Get the lower triangle of this matrix, which effectively sets fsd
        # to zero for size classes larger or equal to the current one
        return np.tril(fsd, k=-1)

    @staticmethod
    def _set_k_diss(k_diss_av: float,
                    scaling_method: str,
                    psd: npt.NDArray[np.float64],
                    n_size_classes: int,
                    gamma: float = 1.0) -> npt.NDArray[np.float64]:
        """
        Set the dissolution rate for each of the size classes,
        based on an average dissolution rate

        Parameters
        ----------
        k_diss_av : float
            Average dissolution rate across size classes.
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
            Dissolution rates for all size classes

        Notes
        -----
        At the moment, we are assuming spherical particles when scaling
        by surface area. This might change.
        """
        # What scaling method has been chosen?
        if scaling_method == 'constant':
            # Use the average k_diss across all size classes
            k_diss = np.full((n_size_classes,), k_diss_av)
        elif scaling_method == 'surface_area':
            # Scale the dissolution according to surface area per unit
            # volume, assuming our particles are spheres
            k_diss = k_diss_av * FragmentMNP._f_surface_area(psd, gamma)
        else:
            # We shouldn't get here, if validation has been performed!
            raise ValueError('Invalid k_diss_scaling_factor provided: ',
                             {scaling_method})
        return k_diss

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
            data = validation.validate_data(data)
        except SchemaError as err:
            raise SchemaError('Input data did not pass validation!') from err
        # Return the config and data with filled defaults
        return config, data
