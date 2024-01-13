"""
FRAGMENT-MNP model
==================
Mechanistic model of Micro and NanoPlastic FRAGMentation in the ENvironmenT.
"""
from importlib.metadata import version

# Populate the package namespace
from .fragmentmnp import FragmentMNP

# Get the version from the installed package metadata
__version__ = version(__name__)

# Let type checkers know what is part of the package
__all__ = ['FragmentMNP']
