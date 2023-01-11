import geopandas as gpd 
import descartes
import matplotlib.pyplot as plt
import csv
import nasa_wildfires as fires
import json

def plot_points():
    lat=[]
    lon=[]
    with open('2021-Incendios.csv', mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            lat.append(float(row["ï»¿latitude"]))
            lon.append(float(row["longitude"]))
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            print(f'\t{row["ï»¿latitude"]} {row["longitude"]}')


            line_count += 1
        print(f'Processed {line_count} lines.')
        return lat,lon

def active_fire():
    data = fires.get_modis(region="central-america")  
    # Transform json input to python objects
    #input_dict = json.loads(data)

    # Filter python objects with list comprehensions
    print(data[3])
    output_dict = [x for x in data if x[2] <= -103.6858]


    print(output_dict)  



gdf = gpd.read_file('D:\Escritorio\Modular\Primavera.shp')
lat, lon = plot_points()
fig, axes = plt.subplots(figsize=(10,8))
gdf.plot(ax = axes, color = "white", edgecolor = "black")
active_fire()
axes.scatter(lon, lat)
plt.show()


