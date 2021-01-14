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

    >>> import folium
    >>> armageddon.plot_circle(52.79, -2.95, 1e3, map=None)
    """

    if not map:
        map = folium.Map(location=[lat, lon], control_scale=True)
    
    color_dict = {4: 'red', 3: 'purple', 2: 'blue', 1: 'green'}
    folium.Circle([lat, lon], radius, color=color_dict[level], fill_color=color_dict[level], fill=True, fillOpacity=0.8, **kwargs).add_to(map)
    return map


def plot_results(burst_lat, burst_lon, blast_lat, blast_lon, radius_list, postcodes, population, sector=False):
    map = folium.Map(location=[blast_lat, blast_lon], control_scale=True)
    folium.PolyLine([
                    [burst_lat, burst_lon],
                    [blast_lat, blast_lon]
                    ], color='black').add_to(map)
    radius_list.sort(reverse=True)
    for i in range(len(radius_list)):
        map = plot_circle(blast_lat, blast_lon, radius_list[i], i+1, map)
    # postcodes, population corresponding
    postcodes_df = pd.read_csv('./armageddon/resources/full_postcodes.csv')
    
    if not sector:
        # plot units
        for i in range(len(postcodes)):
            for j in range(len(postcodes[i])):
                # 这里unit postcode无空格
                row_select = postcodes_df[postcodes_df['Postcode'] == postcodes[i][j]]
                marker_lat = row_select.iloc[0]['Latitude']
                marker_lon = row_select.iloc[0]['Longitude']
                folium.Marker([marker_lat, marker_lon],
                              popup='Unit: ' + postcodes[i][j] + '<br>All usual residents: ' + str(population[i][j]),
                              tooltip="Click me!").add_to(map)
    else:
        # plot sector
        for i in range(len(postcodes)):
            for j in range(len(postcodes[i])):
                # 这里从上一个方法结果加空格后的sector postcode
                sector_postcode = postcodes[i][j][:4] + postcodes[i][j][5:]
                units_in_sector = postcodes_df[postcodes_df['Postcode'].str.contains(sector_postcode)]
                # 算sector平均坐标
                marker_lat = units_in_sector['Latitude'].mean()
                marker_lon = units_in_sector['Longitude'].mean()
                folium.Marker([marker_lat, marker_lon],
                              popup='Sector: ' + postcodes[i][j] + '<br>All usual residents: ' + str(population[i][j]),
                              tooltip="Click me!").add_to(map)

    map.save("index.html")
    return map
