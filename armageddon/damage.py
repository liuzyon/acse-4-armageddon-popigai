import numpy as np
import pandas as pd

from armageddon import locator


def p(r, E, z, p_t=0):
    """
    The exact function to evaluate pressure at given radius, burst energy
    and burst altitude

    Parameters
    ----------
    r: array-like, float
        surface distance from surface zero location, in meter
    E: float
        explosive/burst energy, in kiloton of TNT
    z: float
        burst altitude, in meter
    p_t: float
        the target pressure that wish to find r for,
        =0 means exact function for p(r, E, z)

    Returns
    -------
    pressure: float
        the pressure under given conditions, in Pa
    """
    im = (r**2 + z**2) / E**(2/3)  # intermediate value
    return 3.14e11 * im**-1.3 + 1.8e7 * im**-0.565 - p_t


def bisection(a, b, E, z, p_t, atol=1e-6, max_iter=100):
    """
    Implement the bisection root-finding method on p()

    Parameters
    ----------
    a, b are the initial interval
    E, z, p_t are the same as p()
    atol is a stopping tolerance
    max_iter is the maximum number of iterations allowed

    Returns
    -------
    the final estimate of the root
    """
    # ensure a and b has different sign
    d_ab = np.abs(b - a)
    while np.sign(p(a, E, z, p_t)) == np.sign(p(b, E, z, p_t)):
        a += d_ab
        b += d_ab

    n = 0
    while n <= max_iter:
        c = (a+b) / 2.
        p_c = p(c, E, z, p_t)

        if p_c == 0. or (b-a)/2. < atol:
            return c
        n += 1

        if np.sign(p_c) == np.sign(p(a, E, z, p_t)):
            a = c
        else:
            b = c

    raise RuntimeError('Hit maximum number of iterations with no root found')


def find_r_bisect(E, z, p_l):
    """
    Return radius for each input damage level using bisection method

    Parameters
    ----------
    E, z: float
        same as function p()
    p_l: array-like, float
        the threshold value for each inspecting damage level

    Returns
    -------
    radii: list
        radii/radius corresponds to the provided levels
    """
    # parse and check input parameters
    try:
        E = float(E)
    except ValueError:
        print('E must be a number')
    assert E >= 0, 'burst energy cannot be negative'

    try:
        z = float(z)
    except ValueError:
        print('z must be a number')
    assert z >= 0, 'burst altitude cannot be negative'

    # define constant
    min_r = 1e-6
    dr = 2**15

    # if energy is zero, all impact radii is zero
    if E == 0:
        return [0 for _ in p_l]

    max_p = p(min_r, E, z)
    return [bisection(min_r, min_r+dr, E, z, p_t)
            if max_p >= p_t else 0. for p_t in p_l]


def surf_zero_loc(r, elat, elon, bearing, Rp=6.371e6):
    """
    Calculate the latitude and longitude of the surface zero location

    Parameters
    ----------
    r: float
        the burst distance (m)
    elat, elon, bearing are the same as damage_zones()
    Rp: float
        the radius of the spherical Earth

    Return
    ------
    blat, blon is latitude and longitude of the surface zero location (radian)
    """
    # pre-computation
    im = r / Rp  # intermidiate value
    elat, elon, bearing = np.radians([elat, elon, bearing])  # to radian

    blat = np.arcsin(np.sin(elat) * np.cos(im)
                     + np.cos(elat) * np.sin(im) * np.cos(bearing))
    blon = np.arctan(np.sin(bearing) * np.sin(im) * np.cos(elat)
                     / (np.cos(im) - np.sin(elat) * np.sin(blat))) + elon

    return np.degrees([blat, blon])  # to degree


def damage_zones(outcome, lat, lon, bearing, pressures):
    """
    Calculate the latitude and longitude of the surface zero location and the
    list of airblast damage radii (m) for a given impact scenario.

    Parameters
    ----------

    outcome: Dict
        the outcome dictionary from an impact scenario
    lat: float
        latitude of the meteoroid entry point (degrees)
    lon: float
        longitude of the meteoroid entry point (degrees)
    bearing: float
        Bearing (azimuth) relative to north of meteoroid trajectory (degrees)
    pressures: float, arraylike
        List of threshold pressures to define airblast damage levels

    Returns
    -------

    blat: float
        latitude of the surface zero point (degrees)
    blon: float
        longitude of the surface zero point (degrees)
    damrad: arraylike, float
        List of or one blast radius for the input damage levels

    Examples
    --------

    >>> import armageddon
    >>> outcome = {'burst_altitude': 8e3, 'burst_energy': 7e3,
                   'burst_distance': 90e3, 'burst_peak_dedz': 1e3,
                   'outcome': 'Airburst'}
    >>> armageddon.damage_zones(outcome, 52.79, -2.95, 135,
                                pressures=[1e3, 3.5e3, 27e3, 43e3])
    """

    blat, blon = surf_zero_loc(outcome['burst_distance'], lat, lon, bearing)
    damrad = find_r_bisect(outcome['burst_energy'],
                           outcome['burst_altitude'],
                           pressures)

    return blat, blon, damrad


fiducial_means = {'radius': 10, 'angle': 20, 'strength': 1e6,
                  'density': 3000, 'velocity': 19e3,
                  'lat': 51.5, 'lon': 1.5, 'bearing': -45.}
fiducial_stdevs = {'radius': 1, 'angle': 1, 'strength': 5e5,
                   'density': 500, 'velocity': 1e3,
                   'lat': 0.025, 'lon': 0.025, 'bearing': 0.5}


def impact_risk(planet, means=fiducial_means, stdevs=fiducial_stdevs,
                pressure=27.e3, nsamples=100, sector=True):
    """
    Perform an uncertainty analysis to calculate the risk for each affected
    UK postcode or postcode sector

    Parameters
    ----------
    planet: armageddon.Planet instance
        The Planet instance from which to solve the atmospheric entry

    means: dict
        A dictionary of mean input values for the uncertainty analysis. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    stdevs: dict
        A dictionary of standard deviations for each input value. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    pressure: float
        The pressure at which to calculate the damage zone for each impact

    nsamples: int
        The number of iterations to perform in the uncertainty analysis

    sector: logical, optional
        If True (default) calculate the risk for postcode sectors, otherwise
        calculate the risk for postcodes

    Returns
    -------
    risk: DataFrame
        A pandas DataFrame with columns for postcode (or postcode sector) and
        the associated risk. These should be called ``postcode`` or ``sector``,
        and ``risk``.
    """
    gaussian_param = {}
    # generate parameters through Gaussian Distribution ramdomly
    for key in fiducial_means:
        gaussian_param[key] = np.random.normal(fiducial_means.get(key), fiducial_stdevs.get(key), nsamples)

    postcode_all1 = []
    my_postcodelocator = None
    for i in range(nsamples):
        # Solve the system of differential equations for a given impact scenario
        result = planet.solve_atmospheric_entry(radius=gaussian_param['radius'][i], angle=gaussian_param['angle'][i],
                                                strength=gaussian_param['strength'][i],
                                                density=gaussian_param['density'][i],
                                                velocity=gaussian_param['velocity'][i])

        result = planet.calculate_energy(result)  # calculate the kinetic energy lost per unit altitude
        outcome = planet.analyse_outcome(result)  # calculate the impact and airburst stats
        # Calculate surface zero location and the list of airblast damage radii
        blast_lat, blast_lon, damage_rad = damage_zones(outcome, lat=gaussian_param['lat'][i],
                                                        lon=gaussian_param['lon'][i],
                                                        bearing=gaussian_param['bearing'][i],
                                                        pressures=[pressure])
        print(damage_rad)
        my_postcodelocator = locator.PostcodeLocator()  # instance
        # (unit or sector) postcodes within specific distances of input location
        postcodes = my_postcodelocator.get_postcodes_by_radius([blast_lat, blast_lon], radii=damage_rad, sector=sector)
        postcode_all1 += postcodes  # store total postcodes of samples
    postcode_all2 = []
    # change list of list to list
    for item in postcode_all1:
        postcode_all2 += item
    index = set(postcode_all2)
    for i in set(postcode_all2):
        postcode_all2.count(i)  # calculate the postcode appearing times
    probability = {}  # key is postcode/sector, value is probability
    for i in index:
        probability[i] = postcode_all2.count(i) / nsamples

    # postcode/ postcode sectors that had been hit at least once inside the n times simulation
    probability = {key: value for key, value in probability.items() if value != 0}
    risk = []
    sector = []
    unit = []
    for key in probability:
        if sector:
            # sector is true
            sector += [key]
        else:
            unit += [key]
        # calculate risk
        risk += [probability[key] * my_postcodelocator.get_population_of_postcode([[key]], sector=sector)[-1][-1]]
    if sector:
        # sector is true
        return pd.DataFrame({'sector': sector, 'risk': risk})
    else:
        return pd.DataFrame({'postcode': unit, 'risk': risk})
