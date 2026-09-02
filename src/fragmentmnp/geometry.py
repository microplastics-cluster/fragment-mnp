"""
Particle geometry support for FRAGMENT-MNP.

The polymer fragmentation core remains size-resolved and geometry-agnostic.
This module maps the existing size coordinate to physical geometry:

- sphere: size = particle diameter
- fibre:  size = fibre length; fibre diameter is supplied separately

Fibre particles are represented as straight circular cylinders.  The same
geometry object is shared by all polymer components in a simulation so that
component/formulation bookkeeping remains independent of particle shape.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt


Array = npt.NDArray[np.float64]


class ParticleGeometry:
    """Small geometry interface used by the FRAGMENT-MNP core."""

    shape: str = "unknown"
    size_coordinate_name: str = "size"
    release_geometry: str = "unknown"

    def volume(self, sizes: npt.ArrayLike) -> Array:
        raise NotImplementedError

    def surface_area(self, sizes: npt.ArrayLike) -> Array:
        raise NotImplementedError

    def release_radius(self, sizes: npt.ArrayLike) -> Array:
        """
        Characteristic radial diffusion distance used by additive release.

        Sphere: particle radius = diameter / 2.
        Fibre:  fibre radius = fibre diameter / 2 (not fibre length / 2).
        """
        raise NotImplementedError

    def surface_area_to_volume(self, sizes: npt.ArrayLike) -> Array:
        v = self.volume(sizes)
        a = self.surface_area(sizes)
        return np.divide(a, v, out=np.zeros_like(a, dtype=float), where=v > 0.0)

    def metadata(self) -> dict[str, Any]:
        return {
            "shape": self.shape,
            "size_coordinate_name": self.size_coordinate_name,
            "release_geometry": self.release_geometry,
        }


@dataclass(frozen=True)
class SphereGeometry(ParticleGeometry):
    """Legacy FRAGMENT-MNP spherical particle geometry."""

    shape: str = "sphere"
    size_coordinate_name: str = "diameter"
    release_geometry: str = "sphere"

    def volume(self, sizes: npt.ArrayLike) -> Array:
        d = np.asarray(sizes, dtype=float)
        return (4.0 / 3.0) * np.pi * (d / 2.0) ** 3

    def surface_area(self, sizes: npt.ArrayLike) -> Array:
        d = np.asarray(sizes, dtype=float)
        return 4.0 * np.pi * (d / 2.0) ** 2

    def release_radius(self, sizes: npt.ArrayLike) -> Array:
        return np.asarray(sizes, dtype=float) / 2.0


@dataclass(frozen=True)
class FibreGeometry(ParticleGeometry):
    """
    Straight circular-cylinder fibre geometry.

    Parameters
    ----------
    diameter : float or 1-D array
        Fibre diameter [same length units as the configured size classes].
        A scalar applies to every length class. A vector allows diameter to
        vary by fibre-length class.
    include_endcaps : bool, default=True
        Include the two circular fibre ends in the total surface area.

    Notes
    -----
    The FRAGMENT-MNP size coordinate is interpreted as fibre length L.
    Volume and surface area are therefore

        V = pi R^2 L
        A = 2 pi R L + 2 pi R^2       (when endcaps are included)

    Additive diffusion is approximated as radial diffusion in an infinitely
    long cylinder, so the release length scale is R = diameter / 2.
    """

    diameter: float | tuple[float, ...] = 20e-6
    include_endcaps: bool = True
    shape: str = "fibre"
    size_coordinate_name: str = "length"
    release_geometry: str = "cylinder"

    def _diameter_for_sizes(self, sizes: npt.ArrayLike) -> Array:
        sizes_arr = np.asarray(sizes, dtype=float)
        d = np.asarray(self.diameter, dtype=float)
        if d.ndim == 0:
            return np.full_like(sizes_arr, float(d), dtype=float)
        if d.ndim != 1 or d.size != sizes_arr.size:
            raise ValueError(
                "Fibre diameter must be a positive scalar or a vector with "
                "the same length as the particle/fibre size classes."
            )
        return d.astype(float, copy=False)

    def volume(self, sizes: npt.ArrayLike) -> Array:
        L = np.asarray(sizes, dtype=float)
        d = self._diameter_for_sizes(L)
        R = d / 2.0
        return np.pi * R**2 * L

    def surface_area(self, sizes: npt.ArrayLike) -> Array:
        L = np.asarray(sizes, dtype=float)
        d = self._diameter_for_sizes(L)
        R = d / 2.0
        lateral = 2.0 * np.pi * R * L
        if not self.include_endcaps:
            return lateral
        return lateral + 2.0 * np.pi * R**2

    def release_radius(self, sizes: npt.ArrayLike) -> Array:
        L = np.asarray(sizes, dtype=float)
        return self._diameter_for_sizes(L) / 2.0

    def aspect_ratio(self, sizes: npt.ArrayLike) -> Array:
        L = np.asarray(sizes, dtype=float)
        d = self._diameter_for_sizes(L)
        return L / d

    def metadata(self) -> dict[str, Any]:
        md = super().metadata()
        d = np.asarray(self.diameter, dtype=float)
        md.update({
            "diameter": float(d) if d.ndim == 0 else d.tolist(),
            "include_endcaps": bool(self.include_endcaps),
            "assumption": "straight circular cylinder; radial additive diffusion",
        })
        return md


def make_geometry(config: dict | None, n_size_classes: int | None = None) -> ParticleGeometry:
    """Create the requested particle geometry from validated config."""
    cfg = {} if config is None else dict(config)
    shape = str(cfg.get("shape", "sphere")).lower()

    if shape == "sphere":
        return SphereGeometry()

    if shape in {"fibre", "fiber"}:
        if "diameter" not in cfg:
            raise ValueError("Fibre geometry requires particle_geometry.diameter.")
        diameter = cfg["diameter"]
        d_arr = np.asarray(diameter, dtype=float)
        if n_size_classes is not None and d_arr.ndim == 1 and d_arr.size != n_size_classes:
            raise ValueError(
                f"Fibre diameter vector must have length {n_size_classes}; got {d_arr.size}."
            )
        diameter_value: float | tuple[float, ...]
        if d_arr.ndim == 0:
            diameter_value = float(d_arr)
        else:
            diameter_value = tuple(float(x) for x in d_arr)
        return FibreGeometry(
            diameter=diameter_value,
            include_endcaps=bool(cfg.get("include_endcaps", True)),
        )

    raise ValueError("particle_geometry.shape must be 'sphere' or 'fibre' (alias 'fiber').")
