import uuid
import numpy as np
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
             type: str = 'particle_number_conc',
             plot_dissolution: bool = False,
             log_yaxis=False):
        """
        Plot the output data by choosing from a number of
        pre-defined plot types.

        Parameters
        ----------
        type : str, default='particle_number_conc'
            Tells the function what type of plot to produce.
            Either `particle_number_conc`, `mass_conc`,
            or `dissolution_mass_conc`.
        plot_dissolution : bool, default=False
            Should dissolution be plotted on a separate y-axis
        log_yaxis: bool or str, default=False
            True and "log" plots the y axis on a log scale,
            "symlog" plots the y axis on a symlog scale (useful
            if the data contain zeros), and False plots the y
            axis on a linear scale

        Returns
        -------
        matplotlib.figure.Figure
            A Matplotlib figure object containing the
            produced plot
        """
        # Are we plotting number or mass concentrations?
        if type == 'particle_number_conc':
            ylabel = 'Particle number concentration'
            yvals = self.n.T
        elif type == 'mass_conc':
            ylabel = 'Mass concentration'
            yvals = self.c.T
        else:
            # Raise an error if an invalid plot type is given
            raise ValueError(f'Invalid option for plot `type`: {type}.'
                             'Should be `particle_number_conc` or'
                             '`mass_conc`.')
        # Change the colourmap to something useful
        plt.rcParams['axes.prop_cycle'] = \
            plt.cycler('color',
                       plt.cm.viridis(np.linspace(0, 1, self.n_size_classes)))

        # Create the figure and axes - we need to this because we need to
        # add a second axis for the dissolution data
        fig, ax1 = plt.subplots()

        # Add labels to axes
        ax1.set_xlabel('Time')
        ax1.set_ylabel(ylabel)
        # Should the y axis be logged?
        if log_yaxis in [True, 'log']:
            ax1.set_yscale('log')
        elif log_yaxis == 'symlog':
            ax1.set_yscale('symlog')

        # Plot the concentrations
        ax1.plot(self.t, yvals)
        ax1.legend([f'{d:<1g} µm' for d in self.psd])

        # Create and format the dissolution y axis
        if plot_dissolution:
            ax2 = ax1.twinx()
            ax2.set_ylabel('Dissolution mass concentration')
            ax2.set_yscale('log')
            ax2.plot(self.t, np.sum(self.c_diss, axis=0), '--')
            ax2.legend(['c_diss'])
            # Should the y axis be logged?
            if log_yaxis in [True, 'log']:
                ax2.set_yscale('log')
            elif log_yaxis == 'symlog':
                ax2.set_yscale('symlog')
            return fig, (ax1, ax2)
        else:
            return fig, ax1
