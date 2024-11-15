"""
Defines exceptions that might be raised when running the model
"""

class FMNPIncorrectDistributionLength(Exception):
    """
    Raise when a size distribution is encountered that
    isn't the same length as the number of size classes
    """

class FMNPDistributionValueError(Exception):
    """
    Raise an error when a calculated distribution has
    an incorrect value, such as k_frag or k_diss having
    values less than zero 
    """

class FMNPNumericalError(Exception):
    """
    Raise when a numerical solution is not found when
    running the model
    """
