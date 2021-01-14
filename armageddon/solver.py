#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 01:13:54 2021

@author: su
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import data

class Planet():
    """
    The class called Planet is initialised with constants appropriate
    for the given target planet, including the atmospheric density profile
    and other constants
    """

    def __init__(self, atmos_func='exponential', atmos_filename=None,
                 Cd=1., Ch=0.1, Q=1e7, Cl=0, alpha=0.3, Rp=6371e3,
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
            self.rhoa = lambda z: rho0 * np.exp(-z / H)
        elif atmos_func == 'tabular':
            self.rhoa = lambda z: self.extension1(z)
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
        # 角度制转弧度制
        if radians:
            theta0 = angle
        else:
            theta0 = angle * np.pi / 180

        mass = 4 / 3 * density * np.pi * radius**3
        t0 = 0
        vmtzxr0 = np.array([velocity, mass, theta0, init_altitude, 0, radius])
        vmtzxrs_Rk4, t_all = self.Rk4(
            self.f, vmtzxr0, t0, dt, strength, density)
        # analytic
        # vmtzxrs_Rk4, t_all = self.Rk4(self.f_analy, vmtzxr0, t0, dt, strength, density)

        return pd.DataFrame({'velocity': vmtzxrs_Rk4[:-1, 0],
                             'mass': vmtzxrs_Rk4[:-1, 1],
                             'angle': vmtzxrs_Rk4[:-1, 2]*180/np.pi,
                             'altitude': vmtzxrs_Rk4[:-1, 3],
                             'distance': vmtzxrs_Rk4[:-1, 4],
                             'radius': vmtzxrs_Rk4[:-1, 5],
                             'time': t_all[:-1]})

    def min_max_fact(self, a, b):
        """
        Return the largest factor of a that is <= b

        >>> min_max_fact(0.015, 0.01)
        >>> 0.005

        >>> min_max_fact(0.09, 0.01)
        >>> 0.01
        """
        # convert input to integer
        i = 0
        while int(a) != a and int(b) != b:
            i += 1
            a *= 10
            b *= 10

        if a % b == 0:
            return b * 10**-i
        b *= 10**-i

        s = int(np.sqrt(a))
        f_h, f_l = [], []
        while s >= 1:
            if a % s == 0:
                f_l.append(s)
                f_h.append(a / s)
            s -= 1

        f = np.array(f_h[::-1] + f_l) * 10**-i
        return f[f < b][0]

    def Rk4(self, f, y0, t0, dt, strength, density):
        y = np.array(y0)
        t = np.array(t0)
        y_all = [y0]
        t_all = [t0]
        dt0 = self.min_max_fact(dt, 1e-2)
        result = dt / dt0
        count = 0

        while y[1] >= 0 and y[3] >= 0 and y[2] > 0:
            count += 1
            k1 = dt0 * f(t, y, strength, density)
            k2 = dt0 * f(t + 0.5 * dt0, y + 0.5 * k1, strength, density)
            k3 = dt0 * f(t + 0.5 * dt0, y + 0.5 * k2, strength, density)
            k4 = dt0 * f(t + dt0, y + k3, strength, density)
            y = y + (1. / 6.) * (k1 + 2 * k2 + 2 * k3 + k4)
            t = t + dt0

            if count % result == 0.0:
                y_all.append(y)
                t_all.append(t)

        return np.array(y_all), np.array(t_all)

    def f(self, t, vmtzxrs, strength, density):
        f = np.zeros_like(vmtzxrs)
        v, m, theta, z, x, r = vmtzxrs
        A = np.pi * r ** 2
        rhoa = self.rhoa(z)
        f[0] = (-self.Cd * rhoa * A * (v ** 2)) / \
            (2 * m) + self.g * np.sin(theta)
        f[1] = -self.Ch * self.rhoa(z) * A * v ** 3 / (2 * self.Q)
        f[2] = (self.g * np.cos(theta) / v) - (self.Cl * rhoa * A *
                                               v / (2 * m)) - (v * np.cos(theta) / (self.Rp + z))
        f[3] = -v * np.sin(theta)
        f[4] = v * np.cos(theta) / (1 + z / self.Rp)
        if rhoa * v**2 < strength:
            f[5] = 0
        else:
            f[5] = np.sqrt(7 / 2 * self.alpha * rhoa / density) * v
        return f

    def f_analy(self, t, vmtzxrs, strength, density):
        f = np.zeros_like(vmtzxrs)

        v, m, theta, z, x, r = vmtzxrs
        A = np.pi * r ** 2

        f[0] = -self.Cd * self.rhoa(z) * A * v**2 / (2 * m)
        f[1] = 0
        f[2] = 0
        f[3] = -v * np.sin(theta)
        f[4] = v * np.cos(theta)
        f[5] = 0

        return f

    #using the method of bisection to seach dataframe quickly
    def extension1(self, s):
        df = pd.read_csv('data/AltitudeDensityTable.csv', sep=' ', skiprows=6,
                         names=['Altitude', 'Density', 'Height'])
        data = df['Altitude'].tolist()
        low = 0
        high = len(data)

        while low < high:
            mid = int((low + high) / 2)
            if data[mid] < s:
                low = mid + 1
            else:
                high = mid

        if s <= data[0]:
            x = -1
            y = -1
            z = -1
        elif s >= data[-1]:
            t = len(data) - 1
            x, y, z = df.iloc[-1]
        elif s > data[high]:
            x, y, z = df.iloc[high]
        else:
            x, y, z = df.iloc[high - 1]

        return y * np.exp((x - s) / z)


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
            ``burst_peak_dedz``, ``burst_altitude``, ``burst_distance``,
             ``burst_energy``
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
        init_KE = 1 / 2 * init_mass * init_velocity ** 2 / (4.184e12)
        residual_KE = 1 / 2 * burst_mass * burst_velocity ** 2 / (4.184e12)
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

    def plot(self, result):
        # 我们用result来自calculate_energy
        # 需要在jupyter notebook上展示给客户。
        result = result.copy()
        t_list = result['time'].tolist()
        v_list = result['velocity'].tolist()
        m_list = result['mass'].tolist()
        an_list = result['angle'].tolist()
        al_list = result['altitude'].tolist()
        d_list = result['distance'].tolist()
        r_list = result['radius'].tolist()
        e_list = result['dedz'].tolist()

        fig, axs = plt.subplots(7, 1, figsize=(8, 12))
        fig.tight_layout(w_pad=5, h_pad=5)
        axs[0].plot(t_list, v_list, 'b', label='velocity')
        axs[0].set_xlabel(r'$t$', fontsize=14)
        axs[0].set_ylabel(r'$value$', fontsize=14)
        axs[0].set_title('plot of asteroid speed changes ', fontsize=14)
        axs[0].grid(True)
        axs[0].legend(loc='best', fontsize=14)

        axs[1].plot(t_list, m_list, 'b', label='mass')
        axs[1].set_xlabel(r'$t$', fontsize=14)
        axs[1].set_ylabel(r'$value$', fontsize=14)
        axs[1].set_title('plot of asteroid mass changes', fontsize=14)
        axs[1].grid(True)
        axs[1].legend(loc='best', fontsize=14)

        axs[2].plot(t_list, an_list, 'b', label='angle')
        axs[2].set_xlabel(r'$t$', fontsize=14)
        axs[2].set_ylabel(r'$value$', fontsize=14)
        axs[2].set_title('plot of asteroid angle changes', fontsize=14)
        axs[2].grid(True)
        axs[2].legend(loc='best', fontsize=14)

        axs[3].plot(t_list, al_list, 'b', label='altitude')
        axs[3].set_xlabel(r'$t$', fontsize=14)
        axs[3].set_ylabel(r'$value$', fontsize=14)
        axs[3].set_title('plot of asteroid altitude changes', fontsize=14)
        axs[3].grid(True)
        axs[3].legend(loc='best', fontsize=14)

        axs[4].plot(t_list, d_list, 'b', label='distance')
        axs[4].set_xlabel(r'$t$', fontsize=14)
        axs[4].set_ylabel(r'$value$', fontsize=14)
        axs[4].set_title('plot of asteroid distance changes', fontsize=14)
        axs[4].grid(True)
        axs[4].legend(loc='best', fontsize=14)

        axs[5].plot(t_list, r_list, 'b', label='radius')
        axs[5].set_xlabel(r'$t$', fontsize=14)
        axs[5].set_ylabel(r'$value$', fontsize=14)
        axs[5].set_title('plot of asteroid radius changes', fontsize=14)
        axs[5].grid(True)
        axs[5].legend(loc='best', fontsize=14)

        axs[6].plot(t_list, e_list, 'b', label='energy')
        axs[6].set_xlabel(r'$t$', fontsize=14)
        axs[6].set_ylabel(r'$value$', fontsize=14)
        axs[6].set_title('plot of asteroid energy changes', fontsize=14)
        axs[6].grid(True)
        axs[6].legend(loc='best', fontsize=14)

        plt.show()
