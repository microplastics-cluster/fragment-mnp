"""Optional particle geometries used by FRAGMENT-MNP.

The size-class coordinate means sphere diameter, fibre length, or film length.
All lengths are metres when using SI parameters. Shapes share the original
polymer-fragmentation solver; they change volumes, areas and release kernels.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np
import numpy.typing as npt

Array = npt.NDArray[np.float64]


class ParticleGeometry:
    shape: str = 'unknown'
    size_coordinate_name: str = 'size'
    release_geometry: str = 'unknown'

    def volume(self, sizes: npt.ArrayLike) -> Array:
        raise NotImplementedError

    def surface_area(self, sizes: npt.ArrayLike) -> Array:
        raise NotImplementedError

    def release_radius(self, sizes: npt.ArrayLike) -> Array:
        """Primary half-distance for release (legacy name kept for compatibility)."""
        raise NotImplementedError

    def surface_area_to_volume(self, sizes: npt.ArrayLike) -> Array:
        a = self.surface_area(sizes)
        v = self.volume(sizes)
        return np.divide(a, v, out=np.zeros_like(a), where=v > 0.0)

    def metadata(self) -> dict[str, Any]:
        return {'shape': self.shape, 'size_coordinate_name': self.size_coordinate_name,
                'release_geometry': self.release_geometry}


@dataclass(frozen=True)
class SphereGeometry(ParticleGeometry):
    shape: str = 'sphere'
    size_coordinate_name: str = 'diameter'
    release_geometry: str = 'sphere'

    def volume(self, sizes: npt.ArrayLike) -> Array:
        d = np.asarray(sizes, dtype=float)
        return (4.0 / 3.0) * np.pi * (d / 2.0)**3

    def surface_area(self, sizes: npt.ArrayLike) -> Array:
        d = np.asarray(sizes, dtype=float)
        return 4.0 * np.pi * (d / 2.0)**2

    def release_radius(self, sizes: npt.ArrayLike) -> Array:
        return np.asarray(sizes, dtype=float) / 2.0


def _size_vector(value: float | tuple[float, ...], sizes: npt.ArrayLike, name: str) -> Array:
    s = np.asarray(sizes, dtype=float)
    x = np.asarray(value, dtype=float)
    if x.ndim == 0:
        return np.full_like(s, float(x), dtype=float)
    if x.ndim != 1 or x.size != s.size:
        raise ValueError(f'{name} must be a scalar or an array aligned with size classes.')
    return x.astype(float, copy=False)


@dataclass(frozen=True)
class FibreGeometry(ParticleGeometry):
    """Straight circular cylinder; optional axial release from BOTH ends.

    ``end_diffusion=False`` preserves the exact previous radial-only behaviour.
    With end_diffusion=True the separable finite-cylinder solution is used.
    ``include_endcaps`` controls reported geometric area, and must be True if
    end diffusion is enabled. The width of a fibre is its diameter.
    """
    diameter: float | tuple[float, ...] = 20e-6
    include_endcaps: bool = True
    end_diffusion: bool = False
    shape: str = 'fibre'
    size_coordinate_name: str = 'length'
    release_geometry: str = 'cylinder'

    def _diameter_for_sizes(self, sizes: npt.ArrayLike) -> Array:
        return _size_vector(self.diameter, sizes, 'Fibre diameter')

    def volume(self, sizes: npt.ArrayLike) -> Array:
        L = np.asarray(sizes, dtype=float)
        R = self._diameter_for_sizes(L) / 2.0
        return np.pi * R**2 * L

    def surface_area(self, sizes: npt.ArrayLike) -> Array:
        L = np.asarray(sizes, dtype=float)
        R = self._diameter_for_sizes(L) / 2.0
        return 2.0 * np.pi * R * L + (2.0 * np.pi * R**2 if self.include_endcaps else 0.0)

    def release_radius(self, sizes: npt.ArrayLike) -> Array:
        return self._diameter_for_sizes(sizes) / 2.0

    def aspect_ratio(self, sizes: npt.ArrayLike) -> Array:
        L = np.asarray(sizes, dtype=float)
        return L / self._diameter_for_sizes(L)

    def metadata(self) -> dict[str, Any]:
        md = super().metadata()
        d = np.asarray(self.diameter, dtype=float)
        md.update({'diameter': float(d) if d.ndim == 0 else d.tolist(),
                   'include_endcaps': bool(self.include_endcaps),
                   'end_diffusion': bool(self.end_diffusion),
                   'assumption': 'circular cylinder; separable radial and optional axial diffusion'})
        return md


@dataclass(frozen=True)
class FilmGeometry(ParticleGeometry):
    """Rectangular film fragment with independent width and thickness.

    ``sizes`` is film-fragment LENGTH L (not equivalent sphere diameter).
    Width W and thickness h can be scalar or size-class vectors.
    ``exposed_faces``: 2 (both broad faces) or 1 (one broad face; opposite
    face impermeable). ``include_edges`` adds all four narrow edge surfaces to
    the geometric area and corresponding orthogonal diffusion directions.
    Film dimensions are taken as constant during additive release; the model
    still changes the LENGTH class on fragmentation, not the thickness.
    """
    width: float | tuple[float, ...] = 1e-3
    thickness: float | tuple[float, ...] = 20e-6
    exposed_faces: int = 2
    include_edges: bool = True
    shape: str = 'film'
    size_coordinate_name: str = 'length'
    release_geometry: str = 'slab'

    def width_for_sizes(self, sizes: npt.ArrayLike) -> Array:
        return _size_vector(self.width, sizes, 'Film width')

    def thickness_for_sizes(self, sizes: npt.ArrayLike) -> Array:
        return _size_vector(self.thickness, sizes, 'Film thickness')

    def volume(self, sizes: npt.ArrayLike) -> Array:
        L = np.asarray(sizes, dtype=float)
        return L * self.width_for_sizes(L) * self.thickness_for_sizes(L)

    def surface_area(self, sizes: npt.ArrayLike) -> Array:
        L = np.asarray(sizes, dtype=float)
        W = self.width_for_sizes(L)
        h = self.thickness_for_sizes(L)
        broad = self.exposed_faces * L * W
        return broad + (2.0 * h * (L + W) if self.include_edges else 0.0)

    def release_radius(self, sizes: npt.ArrayLike) -> Array:
        """Half-slab distance: h/2 for two faces, h for one face."""
        h = self.thickness_for_sizes(sizes)
        return h / 2.0 if self.exposed_faces == 2 else h

    def metadata(self) -> dict[str, Any]:
        md = super().metadata()
        for key in ('width', 'thickness'):
            x = np.asarray(getattr(self, key), dtype=float)
            md[key] = float(x) if x.ndim == 0 else x.tolist()
        md.update({'exposed_faces': int(self.exposed_faces), 'include_edges': bool(self.include_edges),
                   'assumption': 'rectangular planar film; separable diffusion into exposed faces/edges'})
        return md


def make_geometry(config: dict | None, n_size_classes: int | None = None) -> ParticleGeometry:
    cfg = {} if config is None else dict(config)
    shape = str(cfg.get('shape', 'sphere')).lower()
    if shape == 'sphere':
        return SphereGeometry()

    def dim(name: str):
        if name not in cfg:
            raise ValueError(f'{shape} geometry requires particle_geometry.{name}.')
        a = np.asarray(cfg[name], dtype=float)
        if (a.ndim > 1 or (a.ndim == 1 and n_size_classes is not None
                           and a.size != n_size_classes) or np.any(~np.isfinite(a))
                or np.any(a <= 0)):
            raise ValueError(f'{shape} {name} must be positive scalar or size vector.')
        return float(a) if a.ndim == 0 else tuple(float(v) for v in a)

    if shape in ('fibre', 'fiber'):
        g = FibreGeometry(diameter=dim('diameter'),
                          include_endcaps=cfg.get('include_endcaps', True),
                          end_diffusion=cfg.get('end_diffusion', False))
        if g.end_diffusion and not g.include_endcaps:
            raise ValueError('fibre end_diffusion requires include_endcaps=True.')
        return g
    if shape in ('film', 'mulch', 'mulch_film'):
        exposed = cfg.get('exposed_faces', 2)
        if type(exposed) is not int or exposed not in (1, 2):
            raise ValueError('film exposed_faces must be 1 or 2.')
        return FilmGeometry(width=dim('width'), thickness=dim('thickness'),
                            exposed_faces=exposed, include_edges=cfg.get('include_edges', True))
    raise ValueError("particle_geometry.shape must be 'sphere', 'fibre', or 'film'.")
