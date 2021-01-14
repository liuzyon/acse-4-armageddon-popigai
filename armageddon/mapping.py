import folium
import numpy as np
import pandas as pd


def plot_circle(lat, lon, radius, level, map=None, **kwargs):
    """
    Plot a circle on a map (creating a new folium map instance if necessary).

    Parameters
    ----------

    lat: float
        latitude of circle to plot (degrees)
    lon: float
        longitude of circle to plot (degrees)
    radius: float
        radius of circle to plot (m)
    map: folium.Map
        existing map object

    Returns
    -------

    Folium map object

    Examples
    --------

    # >>> import folium
    # >>> armageddon.plot_circle(52.79, -2.95, 1e3, map=None)
    """

    if not map:
        map = folium.Map(location=[lat, lon], control_scale=True)

    color_dict = {4: 'red', 3: 'purple', 2: 'blue', 1: 'green'}
    folium.Circle([lat, lon], radius, color=color_dict[level], fill_color=color_dict[level], fill=True, fillOpacity=0.8, **kwargs).add_to(map)
    return map


def plot_results(burst_lat, burst_lon, blast_lat, blast_lon, radius_list, postcodes, population, sector=False):
    """
    Plot all the results on a map.

    Parameters
    ----------

    burst_lat: float
        latitude of entry point (degrees)
    burst_lon: float
        longitude of entry point (degrees)
    blast_lat: float
        latitude of surface zero point (degrees)
    blast_lon: float
        longitude of surface zero point (degrees)
    radius_list: list
        list of radii of circle to plot (m)
    postcodes: list
        list of postcodes for levels
    population: list
        list populations for postcodes

    Returns
    -------

    Folium map object

    """
    map = folium.Map(location=[blast_lat, blast_lon], control_scale=True)

    # according to entry point and surface zero point, plot the flight path.
    folium.PolyLine([
                    [burst_lat, burst_lon],
                    [blast_lat, blast_lon]
                    ], color='black').add_to(map)
    radius_list.sort(reverse=True)

    # for each radius level, plot the corresponding circle. (plot order: from low level to high level)
    for i in range(len(radius_list)):
        map = plot_circle(blast_lat, blast_lon, radius_list[i], i+1, map)

    # read the full_postcodes.csv
    postcodes_df = pd.read_csv('./armageddon/resources/full_postcodes.csv')

    if not sector:
        # plot units
        for i in range(len(postcodes)):
            for j in range(len(postcodes[i])):
                # unit postcode input no extra space, look up coordinate in full_postcodes.csv directly.
                row_select = postcodes_df[postcodes_df['Postcode'] == postcodes[i][j]]
                marker_lat = row_select.iloc[0]['Latitude']
                marker_lon = row_select.iloc[0]['Longitude']
                # plot the unit markers.
                folium.Marker([marker_lat, marker_lon],
                              popup='Unit:<br>' + postcodes[i][j] + '<br><br>All usual residents:<br>' + str(population[i][j]),
                              tooltip="Click me!").add_to(map)
    else:
        # plot sector
        for i in range(len(postcodes)):
            for j in range(len(postcodes[i])):
                # sector postcode input has one extra space, eliminate it and look in full_postcodes.csv.
                sector_postcode = postcodes[i][j][:4] + postcodes[i][j][5:]
                units_in_sector = postcodes_df[postcodes_df['Postcode'].str.contains(sector_postcode)]
                # calculate the sector coordinate(average of units coordinates) as the sector marker on the map
                marker_lat = units_in_sector['Latitude'].mean()
                marker_lon = units_in_sector['Longitude'].mean()
                # plot the sector markers.
                folium.Marker([marker_lat, marker_lon],
                              popup='Sector:<br>' + postcodes[i][j] + '<br><br>All usual residents:<br>' + str(population[i][j]),
                              tooltip="Click me!").add_to(map)
    # save the map.
    map.save("index.html")
    return map
