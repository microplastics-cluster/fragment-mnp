"""
Integration tests for the full model
"""
import copy
import numpy as np
from types import MethodType
from fragmentmnp import FragmentMNP
from fragmentmnp.examples import minimal_config, minimal_data, full_data


def test_model_init():
    """
    Test the model is correctly initialised by checking a few
    attributes after creation
    """
    fmnp = FragmentMNP(minimal_config, minimal_data)
    assert (
        fmnp.n_size_classes == minimal_config['n_size_classes'] and
        fmnp.n_timesteps == minimal_config['n_timesteps']
    )


def test_model_run():
    """
    Test the model run by check the outputs
    """
    fmnp = FragmentMNP(minimal_config, minimal_data)
    output = fmnp.run()
    assert (
        np.array_equal(output.t,
                       np.arange(0.5, minimal_config['n_timesteps'] + 0.5))
        and np.allclose(output.c.sum(), 29400.0)
    )


def test_model_init_timestep_t_eval():
    """
    Test that specifying an integer t_eval results in the correct
    monotonically spaced t_eval timesteps
    """
    config_t_eval = copy.deepcopy(minimal_config)
    config_t_eval['solver_t_eval'] = 'timesteps'
    fmnp = FragmentMNP(config_t_eval, minimal_data)
    np.testing.assert_array_equal(fmnp.t_eval,
                                  np.arange(0.5, config_t_eval['n_timesteps'] + 0.5))


def test_fsd_equal_split():
    """
    Test the fragment size distribution is calculated correctly
    when an even split amongst daughter size classes is assumed
    """
    N = minimal_config['n_size_classes']
    fsd = np.zeros((N, N))
    for i in np.arange(N):
        fsd[i, :] = 1 / i if i != 0 else 0
    fsd = np.tril(fsd, k=-1)
    fmnp = FragmentMNP(minimal_config, minimal_data)
    # Do the assertion
    np.testing.assert_array_equal(fmnp.fsd, fsd)


def test_fsd_rows_sum_to_unity():
    """
    Test the fragment size distribution is calculated correctly
    and that rows sum to unity when a beta parameter other than
    zero is used
    """
    # Set a beta parameter that isn't zero
    psd = np.logspace(-9, -3, 7)
    fsd = FragmentMNP.set_fsd(n=7, psd=psd, beta=-0.1)
    # Check the calculated fsd by summing each row (except the
    # first, i.e. the smallest size class) and checking they 
    # all sum to unity
    assert np.all(np.sum(fsd, axis=1)[1:] == 1.0)


def test_f_surface_area():
    """
    Test the surface area scaling factor is calculated correctly
    """
    # With default parameters, f_surface_area should equal this
    f_surface_area = np.logspace(3, -3, 7)
    fmnp = FragmentMNP(minimal_config, minimal_data)
    model_f_surface_area = fmnp._f_surface_area(fmnp.psd)
    # Do the assertion
    np.testing.assert_allclose(model_f_surface_area, f_surface_area)


def test_k_dists_default_to_constant():
    """
    Test that not providing (t,s) distribution params results
    in constant k_frag and k_diss values over time
    """
    fmnp = FragmentMNP(minimal_config, full_data)
    # Only check all but the smallest size class, as the smallest size
    # class should be zero
    np.testing.assert_array_equal(fmnp.k_frag[1:, :],
                                  full_data['k_frag']['k_f'])
    np.testing.assert_array_equal(fmnp.k_diss, full_data['k_diss']['k_f'])


def test_smallest_size_class_doesnt_fragment():
    """
    Test that the smallest size class has a k_frag equal to
    0, such that there is no fragmentation from it
    """
    fmnp = FragmentMNP(minimal_config, minimal_data)
    np.testing.assert_equal(fmnp.k_frag[0, :], 0.0)


def test_k_frag_linear():
    """
    Test that providing only linear params for the k_frag distribution
    results in linearly increasing k_frag values
    """
    data = copy.deepcopy(full_data)
    data['k_frag']['A_t'] = [1.0]
    fmnp = FragmentMNP(minimal_config, data)
    # Create the same distribution that k_frag should now have,
    # which based on just A_t being set, should be k_frag = 1 * t^1 = t,
    # making sure to set the smallest size class to k_frag=0
    _, t = np.meshgrid(fmnp.surface_areas, fmnp.t_grid, indexing='ij')
    k_frag = data['k_frag']['k_f'] * t
    k_frag[0, :] = 0.0
    # Check its the same
    np.testing.assert_equal(fmnp.k_frag, k_frag)


def test_k_frag_logistic():
    """
    Test that providing logistic params for the k_frag distribution
    results in a logistic regression
    """
    data = copy.deepcopy(full_data)
    data['k_frag']['D_t'] = 1.0
    fmnp = FragmentMNP(minimal_config, data)
    # Create the same distribution that k_frag should now have,
    _, t = np.meshgrid(fmnp.surface_areas, fmnp.t_grid, indexing='ij')
    delta2 = (t.max() - t.min()) / 2
    k_frag = data['k_frag']['k_f'] * (1.0 / (1 + np.exp(-t + delta2)))
    k_frag[0, :] = 0.0
    # Check its the same
    np.testing.assert_equal(fmnp.k_frag, k_frag)


def test_mass_to_particle_number_overload():
    """
    Testing that overloaded the fragmentmnp.mass_to_particle_number
    function works as expected.
    """
    fmnp = FragmentMNP(minimal_config, minimal_data)

    def cube_mass_to_particle_number(self, mass):
        # Presuming particles are cubes instead and the particle
        # size distribution is the length of each face
        n = mass / (self.density * self.psd[:, None] ** 3)
        return n

    # Assign this new function
    setattr(fmnp, "mass_to_particle_number",
            MethodType(cube_mass_to_particle_number, fmnp))
    # Run the model and check the mass to number conversion is
    # done correctly
    output = fmnp.run()
    np.testing.assert_array_equal(output.n,
                                  cube_mass_to_particle_number(fmnp, output.c))
