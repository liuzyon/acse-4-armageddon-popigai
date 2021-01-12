import numpy as np
import pandas as pd


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
        List of distances specifying the blast radii for the input damage levels

    Examples
    --------
    >>> import armageddon
    >>> outcome = {'burst_altitude': 8e3, 'burst_energy': 7e3,
                   'burst_distance': 90e3, 'burst_peak_dedz': 1e3,
                   'outcome': 'Airburst'}
    >>> armageddon.damage_zones(outcome, 52.79, -2.95, 135, pressures=[1e3, 3.5e3, 27e3, 43e3])
    """

    # Replace this code with your own. For demonstration we return lat, lon and 1000 m
    blat = lat
    blon = lon
    damrad = find_r(p, outcome['burst_energy'], outcome['burst_altitude'], pressures)

    return blat, blon, damrad

    sin(phi_2) = sin(phi_1)*cos(r/Rp) + cos(phi_1)*sin(r/Rp)*cos(beta)
    tan(lambda_2-lambda_1) = sin(beta)*sin(r/Rp)*cos(phi_1)/(cos(r/Rp)-sin(phi_1)*sin(phi_2))
    # tbc
    phi_2 = arcsin(sin(phi_2))
    lambda_2 = arctan(tan(lambda_2-lambda_1)) + lambda_1


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
    
    return pd.DataFrame({'sector': '', 'risk': 0}, index=range(1))
