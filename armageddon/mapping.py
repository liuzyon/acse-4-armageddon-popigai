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


def plot_results(entry_lat, entry_lon, blast_lat, blast_lon, radius_list, postcodes=None, population=None, sector=False):
    """
    Plot all the results on a map.

    Parameters
    ----------

    entry_lat: float
        latitude of entry point (degrees)
    entry_lon: float
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
    map.add_child(folium.LatLngPopup())

    # according to entry point and surface zero point, plot the flight path.
    folium.PolyLine([
                    [entry_lat, entry_lon],
                    [blast_lat, blast_lon]
                    ], color='black').add_to(map)
    radius_list.sort(reverse=True)

    # plot the entry point marker.
    folium.Marker([entry_lat, entry_lon],
                  popup='Entry Point: (' + str(entry_lat) + ', ' + str(entry_lon) + ')',
                  icon=folium.Icon(color="green", icon="info-sign"),
                  tooltip="Click me!").add_to(map)

    # plot the blast point marker.
    folium.Marker([blast_lat, blast_lon],
                  popup='Surface Zero Location: (' + str(blast_lat) + ', ' + str(blast_lon) + ')',
                  icon=folium.Icon(color="red", icon="info-sign"),
                  tooltip="Click me!").add_to(map)

    # plot the blast point marker.

    # for each radius level, plot the corresponding circle. (plot order: from low level to high level)
    for i in range(len(radius_list)):
        map = plot_circle(blast_lat, blast_lon, radius_list[i], i+1, map)

    # plot the units or sectors marker if need.
    if postcodes is not None:
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
                    if population is not None:
                        folium.Marker([marker_lat, marker_lon],
                                      popup='Unit:<br>' + postcodes[i][j] + '<br><br>All usual residents:<br>' + str(population[i][j]),
                                      tooltip="Click me!").add_to(map)
                    else:
                        folium.Marker([marker_lat, marker_lon],
                                      popup='Unit:<br>' + postcodes[i][j],
                                      tooltip="Click me!").add_to(map)
        else:
            # plot sector
            for i in range(len(postcodes)):
                for j in range(len(postcodes[i])):
                    # sector postcode input may has one extra space, eliminate it and look in full_postcodes.csv.
                    if len(postcodes[i][j]) == 6:
                        postcodes[i][j] = postcodes[i][j][:4] + postcodes[i][j][5:]
                    units_in_sector = postcodes_df[postcodes_df['Postcode'].str.contains(postcodes[i][j])]
                    # calculate the sector coordinate(average of units coordinates) as the sector marker on the map
                    marker_lat = units_in_sector['Latitude'].mean()
                    marker_lon = units_in_sector['Longitude'].mean()
                    # plot the sector markers.
                    if population is not None:
                        folium.Marker([marker_lat, marker_lon],
                                      popup='Sector:<br>' + postcodes[i][j] + '<br><br>All usual residents:<br>' + str(population[i][j]),
                                      tooltip="Click me!").add_to(map)
                    else:
                        folium.Marker([marker_lat, marker_lon],
                                      popup='Sector:<br>' + postcodes[i][j],
                                      tooltip="Click me!").add_to(map)
    # save the map.
    map.save("index.html")
    return map
