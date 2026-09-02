import copy
import numpy as np
import pytest
from schema import SchemaError

from fragmentmnp import FragmentMNP
from fragmentmnp.geometry import FibreGeometry, SphereGeometry
from fragmentmnp.examples import (
    minimal_config,
    minimal_config_with_fibre,
    minimal_data,
    minimal_data_with_fibre,
)
from fragmentmnp.validation import validate_config


def test_default_geometry_is_sphere_and_legacy_formulas_are_preserved():
    fmnp = FragmentMNP(minimal_config, minimal_data)
    assert fmnp.particle_geometry['shape'] == 'sphere'
    np.testing.assert_allclose(
        fmnp.volume(fmnp.psd),
        (4.0 / 3.0) * np.pi * (fmnp.psd / 2.0) ** 3,
    )
    np.testing.assert_allclose(
        fmnp.surface_area(fmnp.psd),
        4.0 * np.pi * (fmnp.psd / 2.0) ** 2,
    )


def test_fibre_geometry_volume_surface_and_release_radius():
    L = np.array([100e-6, 1e-3, 5e-3])
    d = 20e-6
    g = FibreGeometry(diameter=d, include_endcaps=True)
    R = d / 2.0

    np.testing.assert_allclose(g.volume(L), np.pi * R**2 * L)
    np.testing.assert_allclose(
        g.surface_area(L),
        2.0 * np.pi * R * L + 2.0 * np.pi * R**2,
    )
    np.testing.assert_allclose(g.release_radius(L), np.full_like(L, R))


def test_fibre_geometry_without_endcaps_uses_lateral_area_only():
    L = np.array([1e-3, 2e-3])
    d = 10e-6
    R = d / 2.0
    g = FibreGeometry(diameter=d, include_endcaps=False)
    np.testing.assert_allclose(g.surface_area(L), 2.0 * np.pi * R * L)


def test_fibre_diameter_can_be_size_resolved():
    L = np.array([100e-6, 200e-6, 300e-6])
    d = np.array([10e-6, 20e-6, 30e-6])
    g = FibreGeometry(diameter=tuple(d))
    np.testing.assert_allclose(g.release_radius(L), d / 2.0)
    np.testing.assert_allclose(g.volume(L), np.pi * (d / 2.0)**2 * L)


def test_fibre_config_requires_diameter():
    cfg = copy.deepcopy(minimal_config)
    cfg['particle_geometry'] = {'shape': 'fibre'}
    with pytest.raises(SchemaError):
        validate_config(cfg)


def test_fiber_alias_is_canonicalised_to_fibre():
    cfg = copy.deepcopy(minimal_config_with_fibre)
    cfg['particle_geometry']['shape'] = 'fiber'
    validated = validate_config(cfg)
    assert validated['particle_geometry']['shape'] == 'fibre'


def test_fibre_mass_to_particle_number_uses_cylindrical_volume():
    fmnp = FragmentMNP(minimal_config_with_fibre, minimal_data_with_fibre)
    out = fmnp.run()

    L = fmnp.psd
    d = minimal_config_with_fibre['particle_geometry']['diameter']
    V = np.pi * (d / 2.0)**2 * L
    expected = out.c / (minimal_data_with_fibre['density'] * V[:, None])
    np.testing.assert_allclose(out.n, expected)


def test_fibre_fragmentation_polymer_mass_is_conserved():
    out = FragmentMNP(minimal_config_with_fibre, minimal_data_with_fibre).run()
    total = out.c.sum(axis=0) + out.c_diss + out.c_min
    assert np.allclose(total, total[0], rtol=5e-2, atol=1e-10)


def test_fibre_geometry_metadata_is_stored_on_output():
    out = FragmentMNP(minimal_config_with_fibre, minimal_data_with_fibre).run()
    assert out.particle_geometry['shape'] == 'fibre'
    assert out.particle_geometry['size_coordinate_name'] == 'length'
    assert out.particle_geometry['release_geometry'] == 'cylinder'
