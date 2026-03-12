"""
FRAGMENT-MNP output
===================

Provides functionality for processing, summarising, and visualising model
output data.

This file optionally holds additive outputs for multi-additive / multi-pool
post-processing workflows.
"""
import uuid
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt


class FMNPOutput():
    """
    Class that holds output data from the model and provides plotting and
    post-processing utilities.
    """

    __slots__ = [
        't', 'c', 'n', 'c_diss', 'c_min',
        'n_timesteps', 'n_size_classes', 'soln', 'psd', 'id',
        'c_chem_part_species', 'c_chem_medium_species',
        'c_chem_part_total', 'c_chem_medium_total',
        'species_names', 'additive_names',
        'c_chem_part', 'c_chem_medium', 'A_part', 'A_aq'
    ]

    def __init__(self,
                 t: npt.NDArray,
                 c: npt.NDArray,
                 n: npt.NDArray,
                 c_diss: npt.NDArray,
                 c_min: npt.NDArray,
                 soln, psd,
                 id=None,
                 c_chem_part_species=None,
                 c_chem_medium_species=None,
                 c_chem_part_total=None,
                 c_chem_medium_total=None,
                 additive_names=None,
                 species_names=None) -> None:
        self.t = t
        self.c = c
        self.n = n
        self.c_diss = c_diss
        self.c_min = c_min
        self.soln = soln
        self.psd = psd

        self.c_chem_part_species = c_chem_part_species
        self.c_chem_medium_species = c_chem_medium_species
        self.c_chem_part_total = c_chem_part_total
        self.c_chem_medium_total = c_chem_medium_total
        self.additive_names = additive_names
        self.species_names = species_names

        # Backward-compatible aliases
        if c_chem_part_total is None:
            self.c_chem_part = None
            self.c_chem_medium = None
            self.A_part = None
            self.A_aq = None
        elif c_chem_part_total.shape[0] == 1:
            self.c_chem_part = c_chem_part_total[0]
            self.c_chem_medium = c_chem_medium_total[0]
            self.A_part = self.c_chem_part
            self.A_aq = self.c_chem_medium
        else:
            self.c_chem_part = c_chem_part_total
            self.c_chem_medium = c_chem_medium_total
            self.A_part = c_chem_part_total
            self.A_aq = c_chem_medium_total

        self.n_timesteps = self.t.shape[0]
        self.n_size_classes = self.c.shape[0]
        self.id = uuid.uuid4() if id is None else id

    # ------------------------------------------------------------------
    # Basic helpers
    # ------------------------------------------------------------------
    def _construct_units(self, units):
        units_out = {}
        if (isinstance(units, dict)
                and {'mass', 'volume', 'time', 'length'} <= units.keys()):
            units_out = units.copy()
        elif isinstance(units, str) and units.lower() == 'si':
            units_out['mass'] = 'kg'
            units_out['volume'] = 'm3'
            units_out['time'] = 's'
            units_out['length'] = 'm'
        elif isinstance(units, str) and units.lower() == 'dim':
            units_out['mass'] = 'mass'
            units_out['volume'] = 'volume'
            units_out['time'] = 'time'
            units_out['length'] = 'length'
        elif units is None:
            return None
        else:
            raise ValueError('`units` parameter of `plot` function must be '
                             'a dictionary with keys {`mass`, `volume`, `time` '
                             'and `length`}, `SI` to denote use of SI '
                             'units, or `dim` to denote use of dimensions.')
        units_out['mass_conc'] = f'{units_out["mass"]}/{units_out["volume"]}'
        units_out['number_conc'] = f'/{units_out["volume"]}'
        units_out['rate'] = f'/{units_out["time"]}'
        return units_out

    def _require_chemical_outputs(self):
        if (self.c_chem_part_species is None or self.c_chem_medium_species is None or
                self.c_chem_part_total is None or self.c_chem_medium_total is None):
            raise ValueError(
                'This output object does not contain additive/chemical release data.'
            )

    def _resolve_species_index(self, species):
        if isinstance(species, int):
            return species
        if self.species_names is None:
            raise ValueError('No species names are stored on this output object.')
        if species not in self.species_names:
            raise KeyError(f'Unknown species name: {species}')
        return self.species_names.index(species)

    def _resolve_additive_index(self, additive):
        if isinstance(additive, int):
            return additive
        if self.additive_names is None:
            raise ValueError('No additive names are stored on this output object.')
        if additive not in self.additive_names:
            raise KeyError(f'Unknown additive name: {additive}')
        return self.additive_names.index(additive)

    def _fraction_released(self, part_by_size, medium):
        initial_total = float(np.sum(part_by_size[:, 0]) + medium[0])
        if np.isclose(initial_total, 0.0):
            return np.zeros_like(medium)
        return medium / initial_total

    def _time_to_fraction(self, fraction_series, target):
        fraction_series = np.asarray(fraction_series, dtype=float)
        if fraction_series[0] >= target:
            return float(self.t[0])
        hits = np.where(fraction_series >= target)[0]
        if hits.size == 0:
            return np.nan
        i = int(hits[0])
        if i == 0:
            return float(self.t[0])
        x0, x1 = fraction_series[i - 1], fraction_series[i]
        t0, t1 = self.t[i - 1], self.t[i]
        if np.isclose(x1, x0):
            return float(t1)
        return float(t0 + (target - x0) * (t1 - t0) / (x1 - x0))

    def _release_flux(self, medium):
        medium = np.asarray(medium, dtype=float)
        if medium.size < 2:
            return np.zeros_like(medium)
        return np.gradient(medium, self.t)

    # ------------------------------------------------------------------
    # Public accessors
    # ------------------------------------------------------------------
    def get_species_index(self, name: str) -> int:
        self._require_chemical_outputs()
        return self._resolve_species_index(name)

    def get_additive_index(self, name: str) -> int:
        self._require_chemical_outputs()
        return self._resolve_additive_index(name)

    def get_species_timeseries(self, species):
        self._require_chemical_outputs()
        i = self._resolve_species_index(species)
        return {
            'particulate_by_size': self.c_chem_part_species[i].copy(),
            'aqueous': self.c_chem_medium_species[i].copy(),
            'particulate_total': self.c_chem_part_species[i].sum(axis=0),
            'total': self.c_chem_part_species[i].sum(axis=0) + self.c_chem_medium_species[i]
        }

    def get_additive_timeseries(self, additive):
        self._require_chemical_outputs()
        i = self._resolve_additive_index(additive)
        return {
            'particulate_by_size': self.c_chem_part_total[i].copy(),
            'aqueous': self.c_chem_medium_total[i].copy(),
            'particulate_total': self.c_chem_part_total[i].sum(axis=0),
            'total': self.c_chem_part_total[i].sum(axis=0) + self.c_chem_medium_total[i]
        }

    # ------------------------------------------------------------------
    # Summary / diagnostics tables
    # ------------------------------------------------------------------
    def summary_records(self, level: str = 'species'):
        self._require_chemical_outputs()
        level = level.lower()
        records = []

        if level == 'species':
            for i, name in enumerate(self.species_names):
                part = self.c_chem_part_species[i]
                aq = self.c_chem_medium_species[i]
                frac = self._fraction_released(part, aq)
                flux = self._release_flux(aq)

                initial_part = float(np.sum(part[:, 0]))
                final_part = float(np.sum(part[:, -1]))
                final_aq = float(aq[-1])
                final_total = final_part + final_aq
                initial_total = float(initial_part + aq[0])
                residual = final_total - initial_total

                records.append({
                    'level': 'species',
                    'name': name,
                    'initial_particulate': initial_part,
                    'final_particulate': final_part,
                    'final_aqueous': final_aq,
                    'initial_total': initial_total,
                    'final_total': final_total,
                    'fraction_released': float(frac[-1]),
                    't10': self._time_to_fraction(frac, 0.10),
                    't50': self._time_to_fraction(frac, 0.50),
                    't90': self._time_to_fraction(frac, 0.90),
                    'peak_release_flux': float(np.max(flux)),
                    't_peak_release_flux': float(self.t[np.argmax(flux)]),
                    'conservation_residual': residual,
                    'conservation_residual_rel': (residual / initial_total if not np.isclose(initial_total, 0.0) else 0.0),
                })

        elif level == 'additive':
            for i, name in enumerate(self.additive_names):
                part = self.c_chem_part_total[i]
                aq = self.c_chem_medium_total[i]
                frac = self._fraction_released(part, aq)
                flux = self._release_flux(aq)

                initial_part = float(np.sum(part[:, 0]))
                final_part = float(np.sum(part[:, -1]))
                final_aq = float(aq[-1])
                final_total = final_part + final_aq
                initial_total = float(initial_part + aq[0])
                residual = final_total - initial_total

                records.append({
                    'level': 'additive',
                    'name': name,
                    'initial_particulate': initial_part,
                    'final_particulate': final_part,
                    'final_aqueous': final_aq,
                    'initial_total': initial_total,
                    'final_total': final_total,
                    'fraction_released': float(frac[-1]),
                    't10': self._time_to_fraction(frac, 0.10),
                    't50': self._time_to_fraction(frac, 0.50),
                    't90': self._time_to_fraction(frac, 0.90),
                    'peak_release_flux': float(np.max(flux)),
                    't_peak_release_flux': float(self.t[np.argmax(flux)]),
                    'conservation_residual': residual,
                    'conservation_residual_rel': (residual / initial_total if not np.isclose(initial_total, 0.0) else 0.0),
                })
        else:
            raise ValueError("`level` must be 'species' or 'additive'.")

        return records

    def summary_dataframe(self, level: str = 'species'):
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError('summary_dataframe() requires pandas.') from exc
        return pd.DataFrame(self.summary_records(level=level))

    def mass_balance_records(self, level: str = 'species'):
        return self.summary_records(level=level)

    def mass_balance_dataframe(self, level: str = 'species'):
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError('mass_balance_dataframe() requires pandas.') from exc
        cols = [
            'level', 'name', 'initial_total', 'final_total',
            'conservation_residual', 'conservation_residual_rel'
        ]
        return pd.DataFrame(self.mass_balance_records(level=level))[cols]

    def size_class_contribution_records(self, level: str = 'species'):
        self._require_chemical_outputs()
        level = level.lower()
        records = []

        if level == 'species':
            names = self.species_names
            particulate = self.c_chem_part_species
            aqueous = self.c_chem_medium_species
        elif level == 'additive':
            names = self.additive_names
            particulate = self.c_chem_part_total
            aqueous = self.c_chem_medium_total
        else:
            raise ValueError("`level` must be 'species' or 'additive'.")

        for i, name in enumerate(names):
            part = particulate[i]
            aq = aqueous[i]
            released_total = float(aq[-1] - aq[0])
            released_by_size = part[:, 0] - part[:, -1]
            for j in range(self.n_size_classes):
                records.append({
                    'level': level,
                    'name': name,
                    'size_class_index': j,
                    'size_class_diameter': float(self.psd[j]),
                    'initial_particulate': float(part[j, 0]),
                    'final_particulate': float(part[j, -1]),
                    'released_mass_contribution': float(released_by_size[j]),
                    'released_mass_contribution_fraction': (
                        float(released_by_size[j] / released_total)
                        if not np.isclose(released_total, 0.0) else 0.0
                    )
                })
        return records

    def size_class_contribution_dataframe(self, level: str = 'species'):
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError('size_class_contribution_dataframe() requires pandas.') from exc
        return pd.DataFrame(self.size_class_contribution_records(level=level))

    # ------------------------------------------------------------------
    # Plotting helpers for additive post-processing
    # ------------------------------------------------------------------
    def plot(self,
             type: str = 'mass_conc',
             plot_dissolution: bool = False,
             plot_mineralisation: bool = False,
             plot_additive: bool = False,
             additive_log_yaxis=False,
             log_yaxis=False,
             units=None,
             cmap='viridis',
             show_legend=True,
             size_classes_to_plot=None,
             additive_index: int | None = None,
             show: bool = False):
        unit_labels = self._construct_units(units)
        if type == 'particle_number_conc':
            ylabel = 'Particle number concentration'
            if unit_labels is not None:
                ylabel += f' [{unit_labels["number_conc"]}]'
            yvals = self.n.T
        elif type == 'mass_conc':
            ylabel = 'Mass concentration'
            if unit_labels is not None:
                ylabel += f' [{unit_labels["mass_conc"]}]'
            yvals = self.c.T
        else:
            raise ValueError(f'Invalid option for plot `type`: {type}. '
                             'Should be `particle_number_conc` or `mass_conc`.')

        if size_classes_to_plot is not None:
            yvals = yvals[:, size_classes_to_plot]

        xlabel = 'Time'
        if unit_labels is not None:
            xlabel += f' [{unit_labels["time"]}]'

        cmap_ = plt.colormaps[cmap]
        plt.rcParams['axes.prop_cycle'] = plt.cycler(
            'color', cmap_(np.linspace(0, 1, yvals.shape[1]))
        )

        fig, ax1 = plt.subplots()
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)

        if log_yaxis in [True, 'log']:
            ax1.set_yscale('log')
        elif log_yaxis == 'symlog':
            ax1.set_yscale('symlog')

        ax1.plot(self.t, yvals)

        if show_legend:
            legend = [f'{d:<1g}' for d in self.psd]
            if unit_labels is not None:
                legend = [f'{sc} {unit_labels["length"]}' for sc in legend]
            legend = np.array(legend)
            if size_classes_to_plot is not None:
                legend = legend[size_classes_to_plot]
            ax1.legend(legend)

        ax2 = None
        if plot_dissolution or plot_mineralisation:
            ax2 = ax1.twinx()
            if not plot_mineralisation:
                ylabel_diss = 'Dissolved mass concentration'
                legend = ['Dissolved']
                ax2.plot(self.t, self.c_diss, c='0.4', ls='--')
            elif not plot_dissolution:
                ylabel_diss = 'Mineralised mass concentration'
                legend = ['Mineralised']
                ax2.plot(self.t, self.c_min, c='0.6', ls=':')
            else:
                ylabel_diss = 'Dissolved and mineralised mass concentration'
                legend = ['Dissolved', 'Mineralised']
                ax2.plot(self.t, self.c_diss, c='0.4', ls='--')
                ax2.plot(self.t, self.c_min, c='0.6', ls=':')
            if unit_labels is not None:
                ylabel_diss += f' [{unit_labels["mass_conc"]}]'
            ax2.set_ylabel(ylabel_diss)
            ax2.legend(legend)
            if log_yaxis in [True, 'log']:
                ax2.set_yscale('log')
            elif log_yaxis == 'symlog':
                ax2.set_yscale('symlog')

        ax3 = None
        if plot_additive:
            if self.c_chem_part is None or self.c_chem_medium is None:
                raise ValueError(
                    'plot_additive=True was requested, but this output does not '
                    'contain chemical release results.'
                )
            if additive_index is not None and self.c_chem_part.ndim == 3:
                A_part = self.c_chem_part[additive_index]
                A_aq = self.c_chem_medium[additive_index]
            else:
                A_part = self.c_chem_part
                A_aq = self.c_chem_medium

            Avals = A_part.T
            if size_classes_to_plot is not None:
                Avals = Avals[:, size_classes_to_plot]
            ax1.plot(self.t, Avals, ls='--', alpha=0.7)

            ax3 = ax1.twinx()
            if ax2 is not None:
                ax3.spines['right'].set_position(('outward', 60))
            ax3.plot(self.t, A_aq, c='k', ls='-.')
            ax3.set_ylabel('Chemical mass concentration in medium')
            if unit_labels is not None:
                ax3.set_ylabel(
                    f'Chemical mass concentration in medium [{unit_labels["mass_conc"]}]'
                )
            if additive_log_yaxis in [True, 'log']:
                ax3.set_yscale('log')
            elif additive_log_yaxis == 'symlog':
                ax3.set_yscale('symlog')

        if show:
            plt.show()

        if plot_additive:
            if ax2 is not None and ax3 is not None:
                return fig, (ax1, ax2, ax3)
            if ax3 is not None:
                return fig, (ax1, ax3)
        if ax2 is not None:
            return fig, (ax1, ax2)
        return fig, ax1

    def plot_release_by_additive(self, additive, show=False):
        self._require_chemical_outputs()
        i = self._resolve_additive_index(additive)
        name = self.additive_names[i]
        part_total = self.c_chem_part_total[i].sum(axis=0)
        aq = self.c_chem_medium_total[i]
        total = part_total + aq

        fig, ax = plt.subplots()
        ax.plot(self.t, part_total, label='Particulate')
        ax.plot(self.t, aq, label='Released to medium')
        ax.plot(self.t, total, '--', label='Total')
        ax.set_xlabel('Time')
        ax.set_ylabel('Additive mass concentration')
        ax.set_title(f'Additive release: {name}')
        ax.legend()
        if show:
            plt.show()
        return fig, ax

    def plot_by_pool(self, additive, show=False):
        self._require_chemical_outputs()
        additive_name = self.additive_names[self._resolve_additive_index(additive)]
        species_indices = [
            i for i, name in enumerate(self.species_names)
            if str(name).startswith(f'{additive_name}:')
        ]
        if not species_indices:
            raise ValueError(f'No pools found for additive {additive_name}.')

        fig, ax = plt.subplots()
        for i in species_indices:
            pool_name = self.species_names[i].split(':', 1)[1]
            ax.plot(self.t, self.c_chem_medium_species[i], label=pool_name)
        ax.set_xlabel('Time')
        ax.set_ylabel('Released mass concentration in medium')
        ax.set_title(f'Release by pool: {additive_name}')
        ax.legend()
        if show:
            plt.show()
        return fig, ax

    def plot_stacked_cumulative_release(self, level='additive', show=False):
        self._require_chemical_outputs()
        level = level.lower()
        if level == 'additive':
            y = self.c_chem_medium_total
            labels = self.additive_names
            title = 'Stacked cumulative release by additive'
        elif level == 'species':
            y = self.c_chem_medium_species
            labels = self.species_names
            title = 'Stacked cumulative release by pool/species'
        else:
            raise ValueError("`level` must be 'species' or 'additive'.")

        fig, ax = plt.subplots()
        ax.stackplot(self.t, *y, labels=labels)
        ax.set_xlabel('Time')
        ax.set_ylabel('Cumulative released mass concentration')
        ax.set_title(title)
        ax.legend(loc='best')
        if show:
            plt.show()
        return fig, ax

    def plot_size_class_contribution(self, target, level='species', show=False):
        self._require_chemical_outputs()
        level = level.lower()
        if level == 'species':
            i = self._resolve_species_index(target)
            name = self.species_names[i]
            part = self.c_chem_part_species[i]
        elif level == 'additive':
            i = self._resolve_additive_index(target)
            name = self.additive_names[i]
            part = self.c_chem_part_total[i]
        else:
            raise ValueError("`level` must be 'species' or 'additive'.")

        released_by_size = part[:, 0] - part[:, -1]
        fig, ax = plt.subplots()
        ax.bar(np.arange(self.n_size_classes), released_by_size)
        ax.set_xlabel('Size class index')
        ax.set_ylabel('Released mass contribution')
        ax.set_title(f'Size-class contribution to release: {name}')
        if show:
            plt.show()
        return fig, ax
