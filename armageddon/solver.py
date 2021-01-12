import numpy as np
import pandas as pd


class Planet():
    """
    The class called Planet is initialised with constants appropriate
    for the given target planet, including the atmospheric density profile
    and other constants
    """

    def __init__(self, atmos_func='constant', atmos_filename=None,
                 Cd=1., Ch=0.1, Q=1e7, Cl=1e-3, alpha=0.3, Rp=6371e3,
                 g=9.81, H=8000., rho0=1.2):
        """
        Set up the initial parameters and constants for the target planet

        Parameters
        ----------
        atmos_func : string, optional
            Function which computes atmospheric density, rho, at altitude, z.
            Default is the exponential function rho = rho0 exp(-z/H).
            Options are 'exponential', 'tabular' and 'constant'

        atmos_filename : string, optional
            Name of the filename to use with the tabular atmos_func option

        Cd : float, optional
            The drag coefficient

        Ch : float, optional
            The heat transfer coefficient

        Q : float, optional
            The heat of ablation (J/kg)

        Cl : float, optional
            Lift coefficient

        alpha : float, optional
            Dispersion coefficient

        Rp : float, optional
            Planet radius (m)

        rho0 : float, optional
            Air density at zero altitude (kg/m^3)

        g : float, optional
            Surface gravity (m/s^2)

        H : float, optional
            Atmospheric scale height (m)

        """

        # Input constants
        self.Cd = Cd
        self.Ch = Ch
        self.Q = Q
        self.Cl = Cl
        self.alpha = alpha
        self.Rp = Rp
        self.g = g
        self.H = H
        self.rho0 = rho0

        # set function to define atmoshperic density
        if atmos_func == 'exponential':
            raise NotImplementedError
        elif atmos_func == 'tabular':
            raise NotImplementedError
        elif atmos_func == 'constant':
            self.rhoa = lambda x: rho0
        else:
            raise NotImplementedError(
                "atmos_func must be 'exponential', 'tabular' or 'constant'")

    def solve_atmospheric_entry(
            self, radius, velocity, density, strength, angle,
            init_altitude=100e3, dt=0.05, radians=False):
        """
        Solve the system of differential equations for a given impact scenario

        Parameters
        ----------
        radius : float
            The radius of the asteroid in meters

        velocity : float
            The entery speed of the asteroid in meters/second

        density : float
            The density of the asteroid in kg/m^3

        strength : float
            The strength of the asteroid (i.e. the maximum pressure it can
            take before fragmenting) in N/m^2

        angle : float
            The initial trajectory angle of the asteroid to the horizontal
            By default, input is in degrees. If 'radians' is set to True, the
            input should be in radians

        init_altitude : float, optional
            Initial altitude in m

        dt : float, optional
            The output timestep, in s

        radians : logical, optional
            Whether angles should be given in degrees or radians. Default=False
            Angles returned in the dataframe will have the same units as the
            input

        Returns
        -------
        Result : DataFrame
            A pandas dataframe containing the solution to the system.
            Includes the following columns:
            'velocity', 'mass', 'angle', 'altitude',
            'distance', 'radius', 'time'
        """

        # Enter your code here to solve the differential equations

        return pd.DataFrame({'velocity': velocity,
                             'mass': np.nan,
                             'angle': angle,
                             'altitude': init_altitude,
                             'distance': 0.0,
                             'radius': radius,
                             'time': 0.0}, index=range(1))

    def calculate_energy(self, result):
        """
        Function to calculate the kinetic energy lost per unit altitude in
        kilotons TNT per km, for a given solution.

        Parameters
        ----------
        result : DataFrame
            A pandas dataframe with columns for the velocity, mass, angle,
            altitude, horizontal distance and radius as a function of time

        Returns : DataFrame
            Returns the dataframe with additional column ``dedz`` which is the
            kinetic energy lost per unit altitude

        """

        # Replace these lines with your code to add the dedz column to
        # the result DataFrame
        result = result.copy()
        result.insert(len(result.columns),
                      'dedz', np.array(np.nan))
        result['dedz'] = abs(((1 / 2) * result['mass'] * result['velocity'] ** 2).diff() / (result['altitude'] / 1000).diff()) / (4.184e12)
        return result

    def analyse_outcome(self, result):
        """
        Inspect a pre-found solution to calculate the impact and airburst stats

        Parameters
        ----------
        result : DataFrame
            pandas dataframe with velocity, mass, angle, altitude, horizontal
            distance, radius and dedz as a function of time

        Returns
        -------
        outcome : Dict
            dictionary with details of the impact event, which should contain
            the key ``outcome`` (which should contain one of the following strings:
            ``Airburst``, ``Cratering`` or ``Airburst and cratering``), as well as
            the following keys:
            ``burst_peak_dedz``, ``burst_altitude``, ``burst_distance``, ``burst_energy``
        """

        outcome = {'outcome': 'Unknown',
                   'burst_peak_dedz': 0.,
                   'burst_altitude': 0.,
                   'burst_distance': 0.,
                   'burst_energy': 0.}

        # get dedz column as a series
        dedz = result.iloc[:, -1]
        outcome['burst_peak_dedz'] = dedz.max()
        # get the index of max dedz
        max_index = dedz.idxmax()
        outcome['burst_distance'] = result.loc[max_index, 'distance']
        burst_altitude = result.loc[max_index, 'altitude']
        outcome['burst_altitude'] = burst_altitude
        burst_mass = result.loc[max_index, 'mass']
        burst_velocity = result.loc[max_index, 'velocity']
        init_mass = result.loc[0, 'mass']
        init_velocity = result.loc[0, 'velocity']

        init_KE = 1/2 * init_mass * init_velocity**2 / (4.184e12)
        residual_KE = 1/2 * burst_mass * burst_velocity**2 / (4.184e12)
        KE_loss = init_KE - residual_KE

        if burst_altitude > 5000:
            outcome['outcome'] = 'Airburst'
            outcome['burst_energy'] = KE_loss
        else:
            if KE_loss > residual_KE:
                outcome['burst_energy'] = KE_loss
            else:
                outcome['burst_energy'] = residual_KE
            if max_index == dedz.size - 1:
                outcome['outcome'] = 'Cratering'
            else:
                outcome['outcome'] = 'Airburst and cratering'
        return outcome
