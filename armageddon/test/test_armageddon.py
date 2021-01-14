from collections import OrderedDict
import pandas as pd
import numpy as np
import operator

from pytest import fixture, mark

# Use pytest fixtures to generate objects we know we'll reuse.
# This makes sure tests run quickly


@fixture(scope='module')
def armageddon():
    import armageddon
    return armageddon


@fixture(scope='module')
def planet(armageddon):
    return armageddon.Planet(atmos_func='constant')


@fixture(scope='module')
def loc(armageddon):
    return armageddon.PostcodeLocator(
        postcode_file='armageddon/resources/full_postcodes.csv',
        census_file='armageddon/resources/population_by_postcode_sector.csv'
    )


@fixture(scope='module')
def result(planet):
    input = {'radius': 1.,
             'velocity': 2.0e4,
             'density': 3000.,
             'strength': 1e5,
             'angle': 30.0,
             'init_altitude': 0.0,
             }

    result = planet.solve_atmospheric_entry(**input)

    return result


@fixture(scope='module')
def outcome(planet, result):
    outcome = planet.analyse_outcome(result=result)
    return outcome


def test_import(armageddon):
    assert armageddon


def test_planet_signature(armageddon):
    inputs = OrderedDict(atmos_func='constant',
                         atmos_filename=None,
                         Cd=1., Ch=0.1, Q=1e7, Cl=1e-3,
                         alpha=0.3, Rp=6371e3,
                         g=9.81, H=8000., rho0=1.2)

    # call by keyword
    armageddon.Planet(**inputs)

    # call by position
    armageddon.Planet(*inputs.values())


def test_attributes(planet):
    for key in ('Cd', 'Ch', 'Q', 'Cl',
                'alpha', 'Rp', 'g', 'H', 'rho0'):
        assert hasattr(planet, key)


def test_solve_atmospheric_entry(result):

    assert type(result) is pd.DataFrame

    for key in ('velocity', 'mass', 'angle', 'altitude',
                'distance', 'radius', 'time'):
        assert key in result.columns


def test_analytical_solution(armageddon):
    # condition for analytical solution
    inputs = OrderedDict(atmos_func='exponential',
                         atmos_filename=None,
                         Cd=1., Ch=0, Q=1e7, Cl=0,
                         alpha=0, Rp=np.inf,
                         g=0, H=8000., rho0=1.2)

    planet = armageddon.Planet(**inputs)

    input = {'radius': 10.,
             'velocity': 19e3,
             'density': 3000.,
             'strength': 1e6,
             'angle': 20.0,
             'init_altitude': 100e3,
             'dt': 0.05
             }

    result = planet.solve_atmospheric_entry(**input)
    m0 = 4/3*np.pi*10**3*3000
    coeff = -1.2 * np.pi * 10**2 * 8000 / 2 / m0 / np.sin(20*np.pi/180)
    C = 19000 / np.exp(coeff*np.exp(-100000/8000))
    z = result.loc[:, 'altitude']

    # analytical solution equation (see in documentation)
    v_analyical = C*np.exp(coeff*np.exp(-z/8000))
    v = result.loc[:, 'velocity']
    assert np.allclose(v_analyical, v)

    m = result.loc[:, 'mass']
    assert np.allclose(m0, m)

    angle = result.loc[:, 'angle']
    assert np.allclose(20, angle)

    radius = result.loc[:, 'radius']
    assert np.allclose(10, radius)


def test_calculate_energy(planet, result):

    energy = planet.calculate_energy(result=result)

    assert type(energy) is pd.DataFrame

    for key in ('velocity', 'mass', 'angle', 'altitude',
                'distance', 'radius', 'time', 'dedz'):
        assert key in energy.columns


def test_analyse_outcome(planet, outcome):

    assert type(outcome) is dict

    for key in ('outcome', 'burst_peak_dedz', 'burst_altitude',
                'burst_distance', 'burst_energy'):
        assert key in outcome.keys()


def test_damage_zones(armageddon):

    outcome = {'burst_peak_dedz': 1000.,
               'burst_altitude': 9000.,
               'burst_distance': 90000.,
               'burst_energy': 6000.,
               'outcome': 'Airburst'}

    blat, blon, damrad = armageddon.damage_zones(
        outcome, 55.0, 0., 135., [27e3, 43e3])

    assert type(blat) is np.float64
    assert type(blon) is np.float64
    assert type(damrad) is list
    assert len(damrad) == 2


@mark.xfail
def test_great_circle_distance(armageddon):

    pnts1 = np.array([[54.0, 0.0], [55.0, 1.0], [54.2, -3.0]])
    pnts2 = np.array([[55.0, 1.0], [56.0, -2.1], [54.001, -0.003]])

    data = np.array([[1.28580537e+05, 2.59579735e+05, 2.25409117e+02],
                     [0.00000000e+00, 2.24656571e+05, 1.28581437e+05],
                     [2.72529953e+05, 2.08175028e+05, 1.96640630e+05]])

    dist = armageddon.great_circle_distance(pnts1, pnts2)

    assert np.allclose(data, dist, rtol=1.0e-4)


def test_locator_postcodes(loc):

    latlon = (51.4981, -0.1773)
    result = loc.get_postcodes_by_radius(latlon, [0.13e3])
    expect = [['SW7 2AZ', 'SW7 2BT', 'SW7 2BU', 'SW7 2DD', 'SW7 5HF', 'SW7 5HG', 'SW7 5HQ']]

    assert type(result) is list
    if len(result) > 0:
        for element in result:
            assert type(element) is list

    assert operator.eq(expect, result)


def test_locator_sectors(loc):

    latlon = (51.4981, -0.1773)
    result = loc.get_postcodes_by_radius(latlon, [0.4e3, 0.2e3], True)
    expect = [['SW7 1', 'SW7 2', 'SW7 3', 'SW7 4', 'SW7 5', 'SW7 9'], ['SW7 1', 'SW7 2', 'SW7 3', 'SW7 4', 'SW7 5', 'SW7 9']]

    assert type(result) is list
    if len(result) > 0:
        for element in result:
            assert type(element) is list

    assert operator.eq(expect, result)


def test_population_units(loc):

    result = loc.get_population_of_postcode([['SW7 2AZ', 'SW7 2BT', 'SW7 2BU', 'SW7 2DD'], ['SA8 3AB', 'SA8 3AD', 'SA8 3AE']])
    expect = [[18, 18, 18, 18], [44, 44, 44]]

    assert type(result) is list
    if len(result) > 0:
        for element in result:
            assert type(element) is list

    assert operator.eq(expect, result)


def test_population_sectors(loc):

    result = loc.get_population_of_postcode([['SW7  2']], True)
    expect = [[2283]]

    assert type(result) is list
    if len(result) > 0:
        for element in result:
            assert type(element) is list

    assert operator.eq(expect, result)
