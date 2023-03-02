import uuid
import numpy.typing as npt
from scipy.integrate import OdeSolution
import matplotlib.pyplot as plt


class FMNPOutput():
    """
    Class that holds output data from the model, and
    includes a number of useful plotting functions

    Parameters
    ----------
    t : np.ndarray, shape (n_timesteps,)
        Time series over which the model was run
    c: np.ndarray, shape (n_size_classes, n_timesteps)
        Mass concentrations for each size class over the time
        series
    n : np.ndarray, shape (n_size_classes, n_timesteps)
        Particle number concentrations for each size class over
        the time series
    c_diss : np.ndarray, shape (n_size_classes, n_timesteps)
        Mass concentrations of dissolved organics
    n_diss : np.ndarray, shape (n_size_classes, n_timesteps)
        Particle number concentrations lost from size classes due
        to dissolution
    soln : Bunch object from scipy.integrate.solve_ivp return
        The solution to the model ODE, passed directly from the
        scipy.integrate.solve_ivp method
    psd : np.ndarray, shape (n_size_classes, )
        Particle size distribution - the average diameters of
        each of the particle size classes
    """

    __slots__ = ['t', 'c', 'n', 'c_diss', 'n_diss', 'n_timesteps',
                 'n_size_classes', 'soln', 'psd', 'id']

    def __init__(self,
                 t: npt.NDArray,
                 c: npt.NDArray,
                 n: npt.NDArray,
                 c_diss: npt.NDArray,
                 n_diss: npt.NDArray,
                 soln, psd,
                 id=None) -> None:
        """
        Initialise the output data object
        """
        # Set the main output data
        self.t = t
        self.c = c
        self.n = n
        self.c_diss = c_diss
        self.n_diss = n_diss
        self.soln = soln
        self.psd = psd
        # Save the number of timesteps and size classes
        self.n_timesteps = self.t.shape[0]
        self.n_size_classes = self.c.shape[0]
        # Set the ID based on what we've been given, or
        # give a unique ID
        if id is None:
            self.id = uuid.uuid4()
        else:
            self.id = id

    def plot(self,
             type : str = 'particle_number_conc',
             options : dict = {}):
        """
        Plot the output data by choosing from a number of
        pre-defined plot types 

        Parameters
        ----------
        type : str, default='particle_number_conc'
            Tells the function what type of plot to produce.
            Either `particle_number_conc`, `mass_conc`,
            or `dissolution_mass_conc`.
        options : dict, default={}
            Options that control the different plots:
            * `'plot_dissolution': True` plots the mass lost
              to dissolution as a separate y-axis (not on
              `dissolution_mass_conc` plot type, default=False)

        Returns
        -------
        matplotlib.figure.Figure
            A Matplotlib figure object containing the
            produced plot
        """
        if type == 'particle_number_conc':
            fig = self._plot_particle_number_conc(options)
        elif type == 'mass_conc':
            pass
        elif type == 'dissolution_mass_conc':
            pass
        else:
            raise Exception(f'Invalid type "{type}" specified for ' +
                            'FMNPOutput plot')


    def _plot_particle_number_conc(self, options):
        pass