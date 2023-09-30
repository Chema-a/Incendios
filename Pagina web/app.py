from flask import Flask, render_template
import folium
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, Point
from shapely.ops import nearest_points
from folium.features import DivIcon

app = Flask(__name__)

# Carga del shapefile
gdf = gpd.read_file('Primavera.shp')

# Definición de la cuadrícula
xmin, ymin, xmax, ymax = -103.6858, 20.5430, -103.4552, 20.7269
grid_size = 0.01
cols = list(np.arange(xmin, xmax, grid_size))
rows = list(np.arange(ymin, ymax, grid_size))
rows.reverse()

polygons = []
for x in cols:
    for y in rows:
        polygons.append(Polygon([(x, y), (x + grid_size, y), (x + grid_size, y - grid_size), (x, y - grid_size)]))

grid = gpd.GeoDataFrame({'geometry': polygons})
intersected = gpd.overlay(gdf, grid, how='intersection')

# Función para encontrar el punto más cercano en una capa GeoDataFrame
def closest_point(point, gdf):
    geom_union = gdf.unary_union
    nearest = nearest_points(point, geom_union)
    match = gdf[gdf.geometry == nearest[1]]
    return match.index[0]

@app.route('/kauil')
def index():
    return render_template('index.html')

@app.route('/kauil/load_wildfires')
def load_wildfires():
    # Replace this with your actual code to load wildfires
    # For demonstration purposes, let's add a few sample wildfire points
    wildfires = [(20.68, -103.6), (20.66, -103.58), (20.71, -103.55)]
    
    m = folium.Map(location=[20.6, -103.57], zoom_start=12)

    folium.GeoJson(gdf).add_to(m)

    closest_points = []
    for p in wildfires:
        closest_points.append(closest_point(Point(p[1], p[0]), intersected))

    for index, row in intersected.iterrows():
        tooltip = f"{index}"
        y_centrum, x_centrum = row.geometry.centroid.y, row.geometry.centroid.x
        if index in closest_points:
            marker = folium.Marker([y_centrum, x_centrum], icon=folium.Icon(color='red', icon='fire'), tooltip=tooltip)
        else:
            marker = folium.Marker([y_centrum, x_centrum], tooltip=tooltip)
        marker.add_to(m)

    m.save('templates/map.html')
    return render_template('map.html')

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')