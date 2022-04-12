import numpy as np
from scipy.integrate import solve_ivp
from schema import SchemaError
from . import validation


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

    __slots__ = ['config', 'data', 'n_size_classes', 'psd', 'n_timesteps']

    def __init__(self,
                 config: dict,
                 data: dict,
                 validate: bool = True) -> None:
        """
        Initialise the model
        """
        # Validate the config and data (if we are meant to)
        if validate:
            config, data = self._validate_config_data(config, data)
        # If we passed validation (or aren't validating), save attributes
        self.config = config
        self.data = data
        # Set the number of particle size classes and number of timesteps
        self.n_size_classes = self.config['n_size_classes']
        self.n_timesteps = self.config['n_timesteps']
        # Set the particle size distribution
        self.psd = self._set_psd()

    def run(self):
        """
        Run the model with the config and data provided
        at initialisation.
        """
        pass

    def _set_psd(self):
        """
        Calculate the particle size distribution based on
        the config options passed in ``self.config``
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

    def _set_k_frag(self, k_frag):
        """
        Set the fragmentation rate ``k_frag`` based on
        the value(s) passed in data
        """
        # TODO....
        pass


    @staticmethod
    def _validate_config_data(config, data):
        """
        Validate the config and data dicts passed to the model
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



# How many size classes and time steps
K = 7
T = 100

# Set up a particle size distribution - not actually needed in this example
d = np.logspace(-3, -9, K)

# Set an arbitrary number concentration of plastic in each size class
n_0 = np.full(K, 42.0)

# Fragmentation rate - set the 1e-3 for each size class. Note that this
# means fragmentation still happens from smallest size class (i.e. the mass
# doesn't balance)
k_frag = np.full(K, 0.01)

# Fragment size distribution matrix - assume fragmentation event results in even
# split between size classes of daughter particles
fsd = np.zeros((K,K))
for k in np.arange(K):
    fsd[k,:] = 1 / (K - k - 1) if (K - k) != 1 else 0
# Get the upper triangle of this matrix, which effectively sets f to zero for
# size classes larger (or equal to) than the current one
fsd = np.triu(fsd, k=1)

# Now we can set up the differential equation and solve using SciPy

# Define the function that satisfies n'(t) = f(t, n)
# i.e. the RHS of our differential eq
def f(t, n):
    # Get number of size classes and create empty result to be filled
    N = n.shape[0]
    dndt = np.empty(N)
    # Loop over the size classes and perform the calculation
    for k in np.arange(N):
        dndt[k] = - k_frag[k] * n[k] + np.sum(fsd[:,k] * k_frag * n)
    # Return the solution for all of the size classes
    return dndt

# Numerically solve this given the initial values for n, over T time steps
soln = solve_ivp(fun=f,
                 t_span=(0, T),
                 y0=n_0,
                 t_eval=np.arange(0, T))

# If k_frag != 0 for the smallest size class, then there will be a loss to the
# system, so keep track of that here
n_loss = np.sum(n_0) - np.sum(soln.y, axis=0)

# Finally we can plot the results

# Create the plot with the first size class and add the size classes
# fig = go.Figure()
# for i in range(0,K):
#   fig.add_trace(go.Scatter(x=soln.t, y=soln.y[i], name=f'{d[i]} m'))
# # Plot the loss with a different style
# fig.add_trace(go.Scatter(x=soln.t, y=n_loss, name='Loss', line={'width': 3, 'dash': 'dash'}))
# fig.update_layout(xaxis_title='Time', yaxis_title='Particle number concentration')
# fig.show()