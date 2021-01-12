"""Module dealing with postcode information."""

from datetime import date
import os
import sys

import numpy as np
import pandas as pd


def great_circle_distance(latlon1, latlon2):
    """
    Calculate the great circle distance (in metres) between pairs of
    points specified as latitude and longitude on a spherical Earth
    (with radius 6371 km).

    Parameters
    ----------

    latlon1: arraylike
        latitudes and longitudes of first point (as [n, 2] array for n points) unit:degree
    latlon2: arraylike
        latitudes and longitudes of second point (as [m, 2] array for m points) unit:degree

    Returns
    -------

    numpy.ndarray
        Distance in metres between each pair of points (as an n x m array)

    Examples
    --------

    >>> import numpy
    >>> fmt = lambda x: numpy.format_float_scientific(x, precision=3)
    >>> with numpy.printoptions(formatter={'all': fmt}): print(great_circle_distance([[54.0, 0.0], [55, 0.0]], [55, 1.0]))
    [[1.286e+05]
     [6.378e+04]]
    """
    # deal with naked list
    if not isinstance(latlon1[0], list):
        latlon = []
        latlon.append(latlon1)
        latlon1 = latlon
    if not isinstance(latlon2[0], list):
        latlon = []
        latlon.append(latlon2)
        latlon2 = latlon

    EARTH_RADIUS = 6371000  # unit is meter
    # convert to array
    latlon1 = np.array(latlon1)
    latlon2 = np.array(latlon2)
    # convert degree to radians
    latlon1 = latlon1 * np.pi / 180
    latlon2 = latlon2 * np.pi / 180
    distance = np.empty((len(latlon1), len(latlon2)), float)
    for i, item_1 in enumerate(latlon1):
        for j, item_2 in enumerate(latlon2):
            if item_1[0] == item_2[0]:
                # same latitude
                distance[i][j] = np.abs(item_1[1] - item_2[1]) * EARTH_RADIUS * np.cos(item_1[0])
            elif item_1[1] == item_2[1]:
                # same longitude
                distance[i][j] = np.abs(item_1[0] - item_2[0]) * EARTH_RADIUS
            else:
                # law of spherical cosines(through %timeit to check, this method is fastest)
                distance[i][j] = np.arccos(
                    np.sin(item_1[0]) * np.sin(item_2[0]) + np.cos(item_1[0]) * np.cos(item_2[0]) *
                    np.cos(np.abs(item_1[1] - item_2[1]))) * EARTH_RADIUS
    return distance

class PostcodeLocator(object):
    """Class to interact with a postcode database file."""

    def __init__(self, postcode_file='./resources/full_postcodes.csv',
                 census_file='./resources/population_by_postcode_sector.csv',
                 norm=great_circle_distance):
        """
        Parameters
        ----------

        postcode_file : str, optional
            Filename of a .csv file containing geographic
            location data for postcodes.

        census_file :  str, optional
            Filename of a .csv file containing census data by postcode sector.

        norm : function
            Python function defining the distance between points in latitude-longitude space.

        """
        self.postcode_file = postcode_file
        self.census_file = census_file
        self.norm = norm

    def get_postcodes_by_radius(self, X, radii, sector=False):
        """
        Return (unit or sector) postcodes within specific distances of 
        input location.

        Parameters
        ----------
        X : arraylike
            Latitude-longitude pair of centre location(一个圈中心坐标) 
        radii : arraylike
            array of radial distances from X
        sector : bool, optional
            if true return postcode sectors, otherwise postcode units

        Returns
        -------
        list of lists
            Contains the lists of postcodes closer than the elements of radii to the location X.


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), [0.13e3])
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), [0.4e3, 0.2e3], True)                                                                 
        """
        res = []

        # Read the file, Data type conversion and prepare data.
        postcodes_df = pd.read_csv(self.postcode_file, usecols=[0])
        coordinates_df = pd.read_csv(self.postcode_file, usecols=[2, 3])
        
        postcodes_array = postcodes_df.values   # units array
        coodinates_array = coordinates_df.values    # unit coordinates array

        # Calculate the distance for all postcodes.
        distances = self.norm(coodinates_array, X)


        # Iterate each radius value in list, each iteration check coordinates of all postcodes.
        # If it's within the circle, add the postcode to the postcodes_ra list for the current radius value.
        for ra in radii:
            postcodes_ra = postcodes_array[distances[:, 0] < ra]
            res.append(postcodes_ra)

        return res

    def get_population_of_postcode(self, postcodes, sector=False):
        """
        Return populations of a list of postcode units or sectors.

        Parameters
        ----------
        postcodes : list of lists
            list of postcode units or postcode sectors
        sector : bool, optional
            if true return populations for postcode sectors, otherwise postcode units

        Returns
        -------
        list of lists
            Contains the populations of input postcode units or sectors


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_population_of_postcode([['SW7 2AZ','SW7 2BT','SW7 2BU','SW7 2DD']])
        >>> locator.get_population_of_postcode([['SW7  2']], True)
        [[2283]]
        """
        census = pd.read_csv(self.census_file, usecols=[1, 4])
        postcodes = pd.read_csv(self.postcode_file, usecols=[0])

        res = []
        for level in postcodes:
            level_list = []
            for postcode in level:
                if sector:
                    row_select = census[census['geography'] == postcode]
                    population = row_select.iloc[0]['Variable: All usual residents; measures: Value'] if row_select.shape[0] > 0 else 0
                else:
                    row_select = census[census['geography'] == postcode[0:5]]
                    if row_select.shape[0] > 0:
                        population_in_sector = row_select.iloc[0]['Variable: All usual residents; measures: Value']
                        units_in_sector = postcodes['Postcode'].str.contains(postcode[0:5])
                        num_units = units_in_sector.shape[0]
                        population = population_in_sector / num_units
                    else:
                        population = 0
                level_list.append(population)
            res.append(level_list)
        return res
