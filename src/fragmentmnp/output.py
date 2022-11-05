import uuid
import numpy.typing as npt
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
    """

    __slots__ = ['t', 'c', 'n', 'c_diss', 'n_diss', 'n_timesteps',
                 'n_size_classes', 'id']

    def __init__(self,
                 t: npt.NDArray,
                 c: npt.NDArray,
                 n: npt.NDArray,
                 c_diss: npt.NDArray,
                 n_diss: npt.NDArray,
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
        # Save the number of timesteps and size classes
        self.n_timesteps = self.t.shape[0]
        self.n_size_classes = self.c.shape[0]
        # Set the ID based on what we've been given, or
        # give a unique ID
        if id is None:
            self.id = uuid.uuid4()
        else:
            self.id = id
