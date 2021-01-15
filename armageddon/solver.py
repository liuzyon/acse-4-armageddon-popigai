#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 01:13:54 2021

@author: su
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Planet():
    """
    The class called Planet is initialised with constants appropriate
    for the given target planet, including the atmospheric density profile
    and other constants
    """

    def __init__(self, atmos_func='exponential', atmos_filename=None,
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
        self.df = None

        # set function to define atmoshperic density
        if atmos_func == 'exponential':
            self.rhoa = lambda z: rho0 * np.exp(-z / H)
        elif atmos_func == 'tabular':
            self.df = pd.read_csv(atmos_filename,
                                  sep=' ', skiprows=6,
                                  names=['Altitude', 'Density', 'Height'])
            self.rhoa = lambda z: self.cal_rho_a(z)
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
        # change to radian if input is angle
        if not radians:
            theta0 = np.radians(angle)

        mass = 4 / 3 * density * np.pi * radius**3
        t0 = 0
        vmtzxr0 = np.array([velocity, mass, theta0, init_altitude, 0, radius])
        vmtzxrs_Rk4, t_all = self.Rk4(
            self.f, vmtzxr0, t0, dt, strength, density)
        # analytic
        # vmtzxrs_Rk4, t_all = self.Rk4(self.f_analy, vmtzxr0, t0, dt, strength, density)

        # change to angle if input is angle
        if not radians:
            vmtzxrs_Rk4_angle = np.degrees(vmtzxrs_Rk4[:-1, 2])

        return pd.DataFrame({'velocity': vmtzxrs_Rk4[:-1, 0],
                             'mass': vmtzxrs_Rk4[:-1, 1],
                             'angle': vmtzxrs_Rk4_angle,
                             'altitude': vmtzxrs_Rk4[:-1, 3],
                             'distance': vmtzxrs_Rk4[:-1, 4],
                             'radius': vmtzxrs_Rk4[:-1, 5],
                             'time': t_all[:-1]})

    def choose_dt0(self, dt, dp=2):
        """
        Choose dt0 that dt is a multiple of it
        but not larger than pre-set value

        Parameters
        ----------
        dt: float
            the output time-step
        dp: int
            the number of decimal place for pre-set dt0

        Return
        ------
        dt0: float
            the RK4 time-step

        Examples
        --------
        >>> choose_dt0(0.015, dp=2)
        >>> 0.0075

        >>> choose_dt0(0.09, dp=2)
        >>> 0.01
        """
        # evaluate pre-set dt0
        dt0 = 10**-dp

        # use dt as dt0 if dt is smaller than dt0
        if dt < dt0:
            return dt

        # convert dt to integer
        dt_copy = dt * 10**dp
        while int(dt_copy) != dt_copy:
            dp += 1
            dt_copy *= 10

        # check dt is a multiple of dt0
        if int(dt_copy) % int(dt0 * 10**dp) == 0:
            return dt0

        # half dt until smaller than dt0
        while dt > dt0:
            dt /= 2
        return dt

    def Rk4(self, f, y0, t0, dt, strength, density):
        """
        An RK4 implementation on the ODE system

        Parameters
        ----------
        f: function
            the ODE system that evaluate (velocity, mass, angle, altitude, distance, radius)
        y0: array-like
            the initial value of (velocity, mass, angle, altitude, distance, radius)
        t0: float
            the initial value of time
        dt: float
            the OUTPUT time-step
        strength, density: the same as solve_atmospheric_entry()

        Returns
        -------
        y: array-like
            the value of (velocity, mass, angle, altitude, distance, radius) at each OUTPUT time-step
        t: array-like
            all corresponding OUTPUT time-steps
        """
        # put initial value in
        y_all = [y0]
        t_all = [t0]

        # choose dt0 based on input dt
        dt0 = self.choose_dt0(dt)

        # find ratio and initial counter
        result = int(dt / dt0)
        count = 0

        # termination conditions (OR):
        #     mass < 0
        #     altitude < 0
        #     angle <= 0
        while y0[1] >= 0 and y0[2] > 0 and y0[3] >= 0:
            # increment counter
            count += 1

            # evaluate four stages
            k1 = dt0 * f(y0, strength, density)
            k2 = dt0 * f(y0 + 0.5 * k1, strength, density)
            k3 = dt0 * f(y0 + 0.5 * k2, strength, density)
            k4 = dt0 * f(y0 + k3, strength, density)

            # update y0 and t0
            y0 = (y0 + (1. / 6.) * (k1 + 2 * k2 + 2 * k3 + k4)).copy()
            t0 += dt0

            # write y0 and t0 to return if count is multiple of result
            if count % result == 0:
                y_all.append(y0)
                t_all.append(t0)

        return np.array(y_all), np.array(t_all)

    def f(self, vmtzxrs, strength, density):
        """
        compute ODE system and d(radius)/dt

        Parameters
        ----------
        vmtzxrs: array-like
            the value of (velocity, mass, angle, altitude, distance, radius) in a list
        strength, density: the same as solve_atmospheric_entry()

        Return
        ------
        f: array-like
            the value of time derivative for each of (velocity, mass, angle, altitude, distance, radius) in a list
        """
        f = np.zeros_like(vmtzxrs)

        # unpack input value
        v, m, theta, z, x, r = vmtzxrs

        # pre-compute values
        rhoa = self.rhoa(z)
        rhoa_A = rhoa * np.pi * r ** 2
        m_2 = 2 * m

        # evaluate ODE system
        f[0] = -self.Cd * rhoa_A * v**2 / m_2 + self.g * np.sin(theta)
        f[1] = -self.Ch * rhoa_A * v**3 / (2 * self.Q)
        f[2] = (self.g * np.cos(theta) / v
                - self.Cl * rhoa_A * v / m_2
                - v * np.cos(theta) / (self.Rp + z))
        f[3] = -v * np.sin(theta)
        f[4] = v * np.cos(theta) / (1 + z / self.Rp)

        # compute d(radius)/dt
        if rhoa * v**2 < strength:
            f[5] = 0
        else:
            f[5] = np.sqrt(3.5 * self.alpha * rhoa / density) * v
        return f

    def f_analy(self, vmtzxrs, strength, density):
        """
        compute ODE system and d(radius)/dt with simplifying assumptions

        Parameters and Return are the same as f()
        """
        f = np.zeros_like(vmtzxrs)

        v, m, theta, z, x, r = vmtzxrs

        f[0] = -self.Cd * self.rhoa(z) * np.pi * r ** 2 * v**2 / (2 * m)
        f[1] = 0
        f[2] = 0
        f[3] = -v * np.sin(theta)
        f[4] = v * np.cos(theta)
        f[5] = 0

        return f

    def cal_rho_a(self, z):
        """
        evaluation atmosphere density from tabular

        Parameters
        ----------
        z: number
            the evaluation altitude

        Return
        ------
        rho_a: float
            the evaluated atmosphere density

        >>> cal_rho_a(0)
        >>> 1.225
        >>> cal_rho_a(100000)
        >>> 4.365910282319956e-07
        """
        z_i, rho_i, H_i = self.df.iloc[min(8600, int(z / 10.))]
        return rho_i * np.exp((z_i-z) / H_i)

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
        result['dedz'] = abs(((1 / 2) * result['mass'] * result['velocity']
                              ** 2).diff() / (result['altitude'] / 1000).diff()) / (4.184e12)
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
        dedz = result.loc[:, 'dedz']
        if dedz.empty:
            return outcome
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
        init_KE = 0.5 * init_mass * init_velocity**2 / 4.184e12
        residual_KE = 0.5 * burst_mass * burst_velocity**2 / 4.184e12
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

    def extension2(self):
        # find best fitted parameter
        data = pd.read_csv('./data/ChelyabinskEnergyAltitude.csv')
        data.rename(columns={'Height (km)': 'altitude', 'Energy Per Unit Length (kt Km^-1)': 'dedz'}, inplace=True)
        data_peak_dedz = data['dedz'].max()
        max_index = data['dedz'].idxmax()
        data_burst_altitude = data.loc[max_index, 'altitude'] * 1e3  # unit m
        print('exact burst altitude and peak dedz:')
        print([data_burst_altitude, data_peak_dedz])

        # guess a initial strength and raduis
        strength = 1e6  # N/m^2
        radius = 20
        # specified parameters
        velocity = 1.92e4  # m/s
        density = 3300  # kg/m 3
        angle = 18.3  # degrees

        raduis_all = []
        strength_all = []
        altitude_all = []
        dedz_all = []
        err1_all = []
        err2_all = []
        err_all = []
        raduis_all.append(radius)
        strength_all.append(strength)
        error2 = 1000
        while error2 > 1:
            # reduce error2
            result = self.solve_atmospheric_entry(
                radius, velocity, density, strength, angle,
                init_altitude=43e3, dt=0.05, radians=False)
            result = self.calculate_energy(result)
            X = result['altitude'] / 1000  # km
            Y = result['dedz']
            altitude_all.append(X)
            dedz_all.append(Y)
            outcome = self.analyse_outcome(result)
            error1 = np.abs(outcome['burst_altitude'] - data_burst_altitude) / 1000
            error2 = np.abs(outcome['burst_peak_dedz'] - data_peak_dedz)
            error = error1 + error2
            err1_all.append(error1)
            err2_all.append(error2)
            err_all.append(error)
            print('error1:')
            print(error1)
            print('error2:')
            print(error2)
            print('present fitted burst_altitude and burst peak dedz:')
            print([outcome['burst_altitude'], outcome['burst_peak_dedz']])

            # narrow error2 also will change error1
            if outcome['burst_peak_dedz'] > data_peak_dedz:
                if error2 > 20:
                    # when error2 is large, step of radius is large too
                    radius -= 0.5
                else:
                    radius -= 0.1
                    strength -= 200
            else:
                if error2 > 20:
                    # when error2 is large, step of radius is large too
                    radius += 0.5
                else:
                    radius += 0.1
                    strength += 200
            raduis_all.append(radius)
            strength_all.append(strength)
        print('min error:')
        min_error = np.min(err_all)
        print(min_error)
        min_error_index = err_all.index(min_error)

        final_radius = raduis_all[min_error_index]
        final_strength = strength_all[min_error_index]
        print('radius')  # final best fitted parameter
        print(final_radius)
        print('strength')  # final best fitted parameter
        print(final_strength)
        final_altitude = altitude_all[min_error_index]
        final_dedz = dedz_all[min_error_index]
        # plot exact energy deposition curve vis fitted curve
        plt.plot(final_altitude, final_dedz, '.', label='fitted energy deposition curve')
        plt.plot(data['altitude'], data['dedz'], '*', label='exact energy deposition curve')
        plt.xlabel('altitude(m)')
        plt.ylabel('Energy Per Unit Length (kt Km^-1)')
        plt.legend()
        plt.show()
