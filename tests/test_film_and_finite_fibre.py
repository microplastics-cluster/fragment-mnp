"""Regression, physics and mass-balance tests for finite fibres and mulch film."""
import copy
import numpy as np
import pytest
from scipy.special import jn_zeros
from schema import SchemaError
from fragmentmnp import FragmentMNP
from fragmentmnp.geometry import FilmGeometry, FibreGeometry, make_geometry
from fragmentmnp.validation import validate_config
from fragmentmnp.examples import (
    minimal_config_with_fibre, minimal_data_with_fibre_additive_analytical,
    minimal_data_with_fibre_additive_numerical,
    minimal_config_with_film, minimal_data_with_film_additive_analytical,
    minimal_data_with_film_additive_numerical,
)


def test_film_volume_area_and_release_span():
    sizes = np.array([1e-3, 5e-3, 1e-2])
    g = FilmGeometry(width=2e-3, thickness=20e-6)
    np.testing.assert_allclose(g.volume(sizes), sizes * 2e-3 * 20e-6)
    np.testing.assert_allclose(g.surface_area(sizes), 2*sizes*2e-3 + 2*20e-6*(sizes+2e-3))
    np.testing.assert_allclose(g.release_radius(sizes), np.full(3, 10e-6))
    g1 = FilmGeometry(width=2e-3, thickness=20e-6, exposed_faces=1, include_edges=False)
    np.testing.assert_allclose(g1.release_radius(sizes), np.full(3, 20e-6))
    np.testing.assert_allclose(g1.surface_area(sizes), sizes*2e-3)


def test_film_width_and_thickness_can_vary_by_length_class():
    sizes = [1e-3, 2e-3]
    g = FilmGeometry(width=(1e-3, 4e-3), thickness=(10e-6, 30e-6))
    np.testing.assert_allclose(g.volume(sizes), [1e-3*1e-3*10e-6, 2e-3*4e-3*30e-6])
    np.testing.assert_allclose(g.release_radius(sizes), [5e-6, 15e-6])


@pytest.mark.parametrize('bad', [
    {'shape':'film', 'width':1e-3},
    {'shape':'film', 'thickness':20e-6},
    {'shape':'film', 'width':0.0, 'thickness':20e-6},
    {'shape':'film', 'width':1e-3, 'thickness':-20e-6},
    {'shape':'film', 'width':[1e-3], 'thickness':20e-6},
    {'shape':'film', 'width':1e-3, 'thickness':20e-6, 'exposed_faces':0},
    {'shape':'film', 'width':1e-3, 'thickness':20e-6, 'include_edges':'yes'},
    {'shape':'film', 'width':1e-3, 'thikness':20e-6},
    {'shape':'fibre', 'diameter':20e-6, 'include_endcaps':False, 'end_diffusion':True},
])
def test_invalid_geometry_rejected(bad):
    config = copy.deepcopy(minimal_config_with_film)
    config['particle_geometry'] = bad
    with pytest.raises(SchemaError):
        validate_config(config)


def test_fibre_aliases_and_legacy_default():
    cfg = copy.deepcopy(minimal_config_with_fibre)
    assert validate_config(cfg)['particle_geometry']['end_diffusion'] is False
    cfg['particle_geometry'].update(shape='fiber', end_diffusion=True)
    assert validate_config(cfg)['particle_geometry']['shape'] == 'fibre'
    assert make_geometry({'shape':'mulch', 'width':1e-3, 'thickness':20e-6}).shape == 'film'


def test_analytical_slab_dirichlet_exact_series():
    t=1e5; a=20e-6; D=1e-16
    Fo=D*t/(a*a)
    n=np.arange(1,3001,dtype=float)
    ref=np.sum(8/(np.pi*np.pi*(2*n-1)**2)*np.exp(-(2*n-1)**2*np.pi*np.pi*Fo/4))
    actual=FragmentMNP._analytical_slab_remaining_fraction(t,a,D,k_m=1e8,n_terms=250)
    assert np.isclose(actual,ref,rtol=3e-3,atol=2e-4)


@pytest.mark.parametrize('km',[0.0, 1e-12, 1e-9, 1e-5])
def test_analytical_numerical_slab_agree_and_no_release(km):
    t=1e5; a=20e-6; D=1e-16
    ref=FragmentMNP._analytical_slab_remaining_fraction(t,a,D,km,n_terms=200)
    got=FragmentMNP._numerical_slab_remaining_fraction(t,a,D,K_pw=1e3,k_m=km,
                       n_r=180,n_substeps=90,theta=1.0)
    assert np.isclose(ref,got,rtol=0.025,atol=0.006), (km,ref,got)


def test_finite_fibre_has_more_release_than_radial_only_when_short():
    cfg0=copy.deepcopy(minimal_config_with_fibre)
    cfg0['particle_size_classes']=[40e-6, 80e-6, 1e-3, 3e-3, 4e-3, 5e-3, 10e-3]
    cfg1=copy.deepcopy(cfg0)
    cfg0['particle_geometry']['end_diffusion']=False
    cfg1['particle_geometry']['end_diffusion']=True
    # This test is about the geometry-aware RELEASE operator, not fragmentation.
    fm0=FragmentMNP(cfg0, copy.deepcopy(minimal_data_with_fibre_additive_analytical))
    fm1=FragmentMNP(cfg1, copy.deepcopy(minimal_data_with_fibre_additive_analytical))
    p={'D_p':1e-16, 'D_w':1e-9, 'K_pw':1.0, 'n_terms':160}
    radial=fm0._analytical_additive_release_fraction(10e-6, 200.0, p, size_index=0)
    finite_short=fm1._analytical_additive_release_fraction(10e-6, 200.0, p, size_index=0)
    finite_long=fm1._analytical_additive_release_fraction(10e-6, 200.0, p, size_index=6)
    assert finite_short > radial
    assert finite_short > finite_long


def test_finite_fibre_numerical_analytical_comparable():
    cfg = copy.deepcopy(minimal_config_with_fibre)
    cfg['particle_geometry']['end_diffusion'] = True
    ma = FragmentMNP(cfg, copy.deepcopy(minimal_data_with_fibre_additive_analytical))
    mn = FragmentMNP(cfg, copy.deepcopy(minimal_data_with_fibre_additive_numerical))
    p = {'D_p':1e-16,'D_w':1e-12,'K_pw':10.,'n_terms':150,
         'n_r':120,'n_substeps':80,'theta':1.0}
    a = ma._analytical_additive_release_fraction(10e-6,3000,p,size_index=1)
    n = mn._numerical_additive_release_fraction(10e-6,3000,p,size_index=1)
    assert np.isclose(a,n,rtol=0.10,atol=0.02), (a,n)


def test_film_width_independent_of_fraction_without_edges():
    p = {'D_p':1e-16,'D_w':1e-9,'K_pw':1e3,'n_terms':120}
    cfg1=copy.deepcopy(minimal_config_with_film)
    cfg2=copy.deepcopy(cfg1)
    cfg1['particle_geometry'].update(width=0.5e-3,include_edges=False)
    cfg2['particle_geometry'].update(width=5e-3,include_edges=False)
    a=FragmentMNP(cfg1,copy.deepcopy(minimal_data_with_film_additive_analytical))
    b=FragmentMNP(cfg2,copy.deepcopy(minimal_data_with_film_additive_analytical))
    r=a.release_radii[3]
    assert a._analytical_additive_release_fraction(r,1000.,p,3) == b._analytical_additive_release_fraction(r,1000.,p,3)


def test_film_smaller_width_increases_release_when_edges_exposed():
    p={'D_p':1e-16,'D_w':1e-12,'K_pw':10.0,'n_terms':120}
    cfg=copy.deepcopy(minimal_config_with_film)
    cfg['particle_geometry']['include_edges']=True
    cfg['particle_geometry']['width']=100e-6
    a=FragmentMNP(cfg,copy.deepcopy(minimal_data_with_film_additive_analytical))
    cfg['particle_geometry']['width']=2e-3
    b=FragmentMNP(cfg,copy.deepcopy(minimal_data_with_film_additive_analytical))
    ia=1; time=2500.
    narrow=a._analytical_additive_release_fraction(a.release_radii[ia],time,p,ia)
    wide=b._analytical_additive_release_fraction(b.release_radii[ia],time,p,ia)
    assert narrow > wide


@pytest.mark.parametrize('name,config,data',[
    ('film analytical',minimal_config_with_film,minimal_data_with_film_additive_analytical),
    ('film numerical',minimal_config_with_film,minimal_data_with_film_additive_numerical),
])
def test_film_end_to_end_additive_mass_balances(name,config,data):
    cfg=copy.deepcopy(config)
    cfg['n_timesteps']=15
    out=FragmentMNP(cfg,copy.deepcopy(data)).run()
    initial=out.c_chem_part_total[0,:,0].sum()+out.c_chem_medium_total[0,0]
    final=out.c_chem_part_total[0,:,-1].sum()+out.c_chem_medium_total[0,-1]
    assert np.isclose(initial,final,rtol=1e-10,atol=1e-12),name
    assert out.particle_geometry['shape']=='film'
    assert out.c.shape[1]==cfg['n_timesteps']


def test_finite_fibre_end_to_end_additive_mass_balances():
    cfg=copy.deepcopy(minimal_config_with_fibre)
    cfg['n_timesteps']=12
    cfg['particle_geometry']['end_diffusion']=True
    data=copy.deepcopy(minimal_data_with_fibre_additive_numerical)
    out=FragmentMNP(cfg,data).run()
    initial=out.c_chem_part_total[0,:,0].sum()+out.c_chem_medium_total[0,0]
    final=out.c_chem_part_total[0,:,-1].sum()+out.c_chem_medium_total[0,-1]
    assert np.isclose(initial,final,rtol=1e-10,atol=1e-12)


def test_film_one_face_releases_more_slowly_than_two_faces():
    cfg_two = copy.deepcopy(minimal_config_with_film)
    cfg_one = copy.deepcopy(cfg_two)
    cfg_one['particle_geometry'].update(exposed_faces=1, include_edges=False)
    cfg_two['particle_geometry']['include_edges'] = False
    data = copy.deepcopy(minimal_data_with_film_additive_analytical)
    p = {'D_p':1e-17,'D_w':1e-10,'K_pw':10.0,'n_terms':150}
    two = FragmentMNP(cfg_two, data)
    one = FragmentMNP(cfg_one, data)
    i=3;time=86400.0
    f_two=two._analytical_additive_release_fraction(two.release_radii[i],time,p,i)
    f_one=one._analytical_additive_release_fraction(one.release_radii[i],time,p,i)
    assert 0 <= f_one < f_two <= 1


def test_film_combines_with_components():
    from fragmentmnp.examples import minimal_data_with_components
    cfg=copy.deepcopy(minimal_config_with_film)
    cfg['n_timesteps']=6
    out=FragmentMNP(cfg,copy.deepcopy(minimal_data_with_components)).run()
    assert out.c_component.shape[0]==3
    np.testing.assert_allclose(out.c_component.sum(axis=0), out.c)
    assert out.particle_geometry['shape']=='film'


def test_film_release_with_time_dependent_diffusivity():
    cfg = copy.deepcopy(minimal_config_with_film)
    cfg['n_timesteps']=12
    d=copy.deepcopy(minimal_data_with_film_additive_analytical)
    d['additives'][0]['pools'][0]['release']['params']['D_p'] = {
        'times': [0.0, 6.0, 12.0],
        'values':[1e-20, 1e-18, 1e-16],
    }
    out=FragmentMNP(cfg,d).run()
    ts=out.get_additive_timeseries(0)
    assert np.all(np.diff(ts['aqueous'])>=-1e-14)
    assert ts['aqueous'][-1]>0
