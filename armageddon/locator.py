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
        latitudes and longitudes of first point (as [n, 2] array for n points)
    latlon2: arraylike
        latitudes and longitudes of second point (as [m, 2] array for m points)

    Returns
    -------

    numpy.ndarray
        Distance in metres between each pair of points (as an n x m array)

    Examples
    --------

    >>> import numpy
    >>> fmt = lambda x: numpy.format_float_scientific(x, precision=3)}
    >>> with numpy.printoptions(formatter={'all', fmt}):
        print(great_circle_distance([[54.0, 0.0], [55, 0.0]], [55, 1.0]))
    [1.286e+05 6.378e+04]
    """

    distance = np.empty((len(latlon1), len(latlon2)), float)
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
        postcodes_df = pd.read_csv(self.postcode_file, usecols=[0])
        coordinates_df = pd.read_csv(self.postcode_file, usecols=[2,3])
        postcodes_df['Postcode'] = postcodes_df['Postcode'].astype(str)
        coordinates_df['Latitude'] = coordinates_df['Latitude'].astype(float)
        coordinates_df['Longitude'] = coordinates_df['Longitude'].astype(float)

        postcodes_array = postcodes_df.values # 各个unit的数组
        coodinates_array = coordinates_df.values #各个unit对应的坐标

        distances = self.norm(coodinates_array, X)

        #对每一个半径进行迭代，每次迭代遍历所有邮编区域，看是否在圈中，如果在，加入当前半径对应的postcodes_ra里
        postcodes_ra = []
        for ra in radii:
            postcodes_ra = postcodes_array[distances < ra]

            # for i in len(distances):
            #     if distances(i) < ra :
            #         postcodes_ra.append(postcodes_array[i])

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
        """
        census = pd.read_csv(self.census_file, usecols=[1, 4])
        census['geography'] = census['geography'].astype(str)
        census['Variable: All usual residents; measures: Value'] = census['Variable: All usual residents; measures: Value'].astype(int)
        # postcodes 里是list的list
        all_list = []
        for element in postcodes:
            if sector:
                rows_select = census.loc[census['geography'].astype(int).isin(element)]
                element_list = rows_select['Variable: All usual residents; measures: Value'].values.tolist()
                all_list.append(element_list)
            else:
                pass
        return all_list
