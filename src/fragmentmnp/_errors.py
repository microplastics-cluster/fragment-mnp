"""
Defines exceptions that might be raised when running the model
"""

class FMNPIncorrectDistributionLength(Exception):
    """
    Raise when a size distribution is encountered that
    isn't the same length as the number of size classes
    """

class FMNPNumericalError(Exception):
    """
    Raise when a numerical solution is not found when
    running the model
    """
