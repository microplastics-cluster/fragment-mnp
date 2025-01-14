"""
FRAGMENT-MNP Output
===================

Provides functionality for processing and visulalising model output data.
"""
import uuid
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt


class FMNPOutput():
    """
    Class that holds output data from the model and provides plotting
    functionalities.

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
    c_diss_from_sc : np.ndarray, shape (n_size_classes, n_timesteps)
        Mass concentrations of dissolved organics lost from each
        size class, not including initial dissolved concentrations
        specified by `initial_concs_diss` parameter
    c_diss : np.ndarray, shape (n_timesteps)
        Total mass concentrations of dissolved organics, including
        initial concentrations
    n_diss_from_sc : np.ndarray, shape (n_size_classes, n_timesteps)
        Particle number concentrations lost from size classes due
        to dissolution, not including initial dissolved
        concentrations specified by `initial_concs_diss` parameter
    soln : Bunch object from scipy.integrate.solve_ivp return
        The solution to the model ODE, passed directly from the
        scipy.integrate.solve_ivp method
    psd : np.ndarray, shape (n_size_classes, )
        Particle size distribution - the average diameters of
        each of the particle size classes
    """

    __slots__ = ['t', 'c', 'n', 'c_diss_from_sc', 'c_diss', 'n_diss_from_sc',
                 'n_timesteps', 'n_size_classes', 'soln', 'psd', 'id']

    def __init__(self,
                 t: npt.NDArray,
                 c: npt.NDArray,
                 n: npt.NDArray,
                 c_diss_from_sc: npt.NDArray,
                 c_diss: npt.NDArray,
                 n_diss_from_sc: npt.NDArray,
                 soln, psd,
                 id=None) -> None:
        """
        Initialise the output data object
        """
        # Set the main output data
        self.t = t
        self.c = c
        self.n = n
        self.c_diss_from_sc = c_diss_from_sc
        self.c_diss = c_diss
        self.n_diss_from_sc = n_diss_from_sc
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
             type: str = 'mass_conc',
             plot_dissolution: bool = False,
             log_yaxis=False,
             units=None,
             cmap='viridis',
             show_legend=True,
             size_classes_to_plot=None):
        """
        Plot the output data by choosing from a number of
        pre-defined plot types.

        Parameters
        ----------
        type : str, default='mass_conc'
            Tells the function what type of plot to produce.
            Either `particle_number_conc` or `mass_conc`.
        plot_dissolution : bool, default=False
            Should dissolution be plotted on a separate y-axis
        log_yaxis: bool or str, default=False
            True and "log" plots the y axis on a log scale,
            "symlog" plots the y axis on a symlog scale (useful
            if the data contain zeros), and False plots the y
            axis on a linear scale
        units: dict or str, default=None
            Units to be used in axis labels. Must either be "SI" to
            use SI units, "dim" to use dimensional labels ("mass",
            "volume", etc), or a dictionary containing elements for
            `mass`, `volume` and `time`. E.g. `units={'mass': 'kg',
            'volume': 'm3', 'time': 's', 'length': 'nm'}`. Note that
            these units are used purely for axis labelling and are
            not used to modify the model output data (which is unit
            agnostic). If `None` (the default), no units are added to
            labels.
        cmap: str, default='viridis'
            The colormap to use for the plot. Must be one of the
            colormaps `available in matplotlib
            <https://matplotlib.org/stable/gallery/color/colormap_reference.html>`.
            Note that these are case-sensitive.
        show_legend: bool, default=True
            Should size classes be shown on a legend?
        size_classes_to_plot: list of ints, default=None
            Only plot specific size classes by providing a list of
            size class indices to plot, where 0 is the index of the
            smallest size class. By default, all size classes are
            plotted.

        Returns
        -------
        (matplotlib.figure.Figure, matplotlib.axes.Axes)
            Matplotlib figure and axes objects for the plot
        """
        unit_labels = self._construct_units(units)
        # Are we plotting number or mass concentrations? Set the
        # y values and labels accordingly
        if type == 'particle_number_conc':
            # Set the yaxis label
            ylabel = 'Particle number concentration'
            if unit_labels is not None:
                ylabel += f' [{unit_labels["number_conc"]}]'
            # Set the yaxis values, transposed to pass to mpl
            yvals = self.n.T
        elif type == 'mass_conc':
            ylabel = 'Mass concentration'
            if unit_labels is not None:
                ylabel += f' [{unit_labels["mass_conc"]}]'
            yvals = self.c.T
        else:
            # Raise an error if an invalid plot type is given
            raise ValueError(f'Invalid option for plot `type`: {type}. '
                             'Should be `particle_number_conc` or'
                             '`mass_conc`.')
        # Only plot the size classes we have been asked to
        if size_classes_to_plot is not None:
            yvals = yvals[:, size_classes_to_plot]
        # Set the xaxis label
        xlabel = 'Time'
        if unit_labels is not None:
            xlabel += f' [{unit_labels["time"]}]'

        # Change the colourmap to something useful (Viridis)
        cmap_ = plt.colormaps[cmap]
        plt.rcParams['axes.prop_cycle'] = \
            plt.cycler('color',
                       cmap_(np.linspace(0, 1, yvals.shape[1])))

        # Create the figure and axes - we need to this because we need to
        # add a second axis for the dissolution data
        fig, ax1 = plt.subplots()

        # Construct the xlabel with(out) units
        xlabel = 'Time'
        if unit_labels is not None:
            xlabel += f' [{unit_labels["time"]}]'
        # Add the x and y labels
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)
        # Should the y axis be logged?
        if log_yaxis in [True, 'log']:
            ax1.set_yscale('log')
        elif log_yaxis == 'symlog':
            ax1.set_yscale('symlog')

        # Plot the concentrations
        ax1.plot(self.t, yvals)

        # Add a size class legend, if requested
        if show_legend:
            # Construct the legend with(out) units
            legend = [f'{d:<1g}' for d in self.psd]
            if unit_labels is not None:
                legend = [f'{sc} {unit_labels["length"]}' for sc in legend]
            legend = np.array(legend)
            # Only show the requested size classes on the legend
            if size_classes_to_plot is not None:
                legend = legend[size_classes_to_plot]
            ax1.legend(legend)

        # Create and format the dissolution y axis
        if plot_dissolution:
            ax2 = ax1.twinx()
            # Construct the ylabel
            ylabel_diss = 'Dissolution mass concentration'
            if unit_labels is not None:
                ylabel_diss += f' [{unit_labels["mass_conc"]}]'
            ax2.set_ylabel(ylabel_diss)
            ax2.plot(self.t, self.c_diss, '--')
            # Always show the dissolution legend to make it clear
            # which line is dissolution
            ax2.legend(['Dissolution'])
            # Should the y axis be logged?
            if log_yaxis in [True, 'log']:
                ax2.set_yscale('log')
            elif log_yaxis == 'symlog':
                ax2.set_yscale('symlog')
            return fig, (ax1, ax2)
        else:
            return fig, ax1

    def _construct_units(self, units):
        """
        Construct strings for axis labels based on the units
        dictionary or string provided. Options given in `plot()`
        method.
        """
        units_out = {}
        if (type(units) == dict
                and {'mass', 'volume', 'time', 'length'} <= units.keys()):
            # If we've been provided a dict of units, use these
            units_out = units.copy()
        elif type(units) == str and units.lower() == 'si':
            # Use SI units
            units_out['mass'] = 'kg'
            units_out['volume'] = 'm3'
            units_out['time'] = 's'
            units_out['length'] = 'm'
        elif type(units) == str and units.lower() == 'dim':
            # Use dimensional labels as the units
            units_out['mass'] = 'mass'
            units_out['volume'] = 'volume'
            units_out['time'] = 'time'
            units_out['length'] = 'length'
        elif units is None:
            # Nothing in, nothing out...
            return None
        else:
            raise ValueError('`units` parameter of `plot` function must be '
                             'a dictionary with keys {`mass`, `volume`, `time` '
                             'and `length`}, `SI` to denote use of SI '
                             'units, or `dim` to denote use of dimensions.')
        # Construct the required compound units
        units_out['mass_conc'] = f'{units_out["mass"]}/{units_out["volume"]}'
        units_out['number_conc'] = f'/{units_out["volume"]}'
        units_out['rate'] = f'/{units_out["time"]}'
        # Return the unit strings
        return units_out
