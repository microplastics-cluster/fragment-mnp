"""
FRAGMENT-MNP model
==================
Mechanistic model of Micro and NanoPlastic FRAGMentation in the ENvironmenT.
"""
from importlib.metadata import PackageNotFoundError, version

from .fragmentmnp import FragmentMNP
from .geometry import ParticleGeometry, SphereGeometry, FibreGeometry, make_geometry

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # Allows the source tree and test suite to run before package installation.
    __version__ = "0+local"

__all__ = [
    'FragmentMNP',
    'ParticleGeometry',
    'SphereGeometry',
    'FibreGeometry',
    'make_geometry',
]
