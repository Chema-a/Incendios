from flask import Flask, render_template,jsonify
import folium
import geopandas as gpd
from shapely.geometry import Polygon, Point
from shapely.ops import nearest_points
import numpy as np
import math
import mpmath
import ee
import requests
import folium
import branca.colormap as cm
import pdb
from time import sleep
from gevent.pywsgi import WSGIServer
import asyncio
import aiohttp
import json
import nasa_wildfires as fires

global intersected_n
app = Flask(__name__)
mpmath.mp.dps = 60
coordinatesPrimavera = [[-103.6858, 20.7269],[-103.4552,20.5430]]
# Download a regional GeoJSON of hotspots detected by the MODIS satellite in a recent 7-day period.
def WildfireHotspots():
    data = fires.get_modis(region="central-america", time_range="7d")  
    lat = []
    lon = []
    # Detectar los incendios en el rango de 7 días dentro del rango de la primavera
    for i in range(len(data["features"])):
        if data["features"][i]["geometry"]["coordinates"][0] >= coordinatesPrimavera[0][0] and data["features"][i]["geometry"]["coordinates"][0] <= coordinatesPrimavera[1][0] and data["features"][i]["geometry"]["coordinates"][1] <= coordinatesPrimavera[0][1] and data["features"][i]["geometry"]["coordinates"][1] >= coordinatesPrimavera[1][1]:
            lon.append(data["features"][i]["geometry"]["coordinates"][0])
            lat.append(data["features"][i]["geometry"]["coordinates"][1])
    
    heat_points = list(zip(lon, lat))  # Unificar las coordenadas en tuplas (x, y)
    return heat_points
# Función para cargar el shapefile
def load_shapefile(file_path):
    try:
        gdf = gpd.read_file(file_path)
        return gdf
    except Exception as e:
        raise Exception(f"Error al cargar el shapefile: {str(e)}")

# Función para definir la cuadrícula
def create_grid(xmin, ymin, xmax, ymax, grid_size):
    cols = list(range(int(xmin), int(xmax), grid_size))
    rows = list(range(int(ymin), int(ymax), -grid_size))
    return cols, rows

def distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def closest_point(p, points):
    min_distance = float('inf')
    for i, row in points.iterrows():
        dist = distance(p[0], p[1], row.geometry.centroid.x, row.geometry.centroid.y)
        if dist < min_distance:
            closest = (row.geometry.centroid.y, row.geometry.centroid.x)
            min_distance = dist
    return closest

def get_intersected_data():
    # Carga del shapefile
    gdf = load_shapefile('Primavera.shp')

        # Definición de la cuadrícula
    xmin, ymin, xmax, ymax = -103.6858, 20.5430, -103.4552, 20.7269

    grid_size = 0.01
    cols = list(np.arange(xmin, xmax, grid_size))
    rows = list(np.arange(ymin, ymax, grid_size))  


    ############################################## Cambiado a -grid_size para que disminuya en lugar de aumentar
    
    
    rows.reverse()

    polygons = []
    for x in cols:
        for y in rows:
            polygons.append(Polygon([(x,y), (x+grid_size,y), (x+grid_size,y-grid_size), (x,y-grid_size)]))


    grid = gpd.GeoDataFrame({'geometry': polygons})
    intersected = gpd.overlay(gdf, grid, how='intersection')

    return intersected

def replace_zero_with_average(arr):
    # Calcula el promedio de todos los registros del arreglo
    total_sum = sum(arr)
    non_zero_count = len(arr) - arr.count(0)
    average = total_sum / non_zero_count if non_zero_count > 0 else 0

    # Recorre el arreglo y reemplaza los valores iguales a 0 por el promedio
    for i in range(len(arr)):
        if arr[i] == 0:
            arr[i] = average

    return arr

def GetSimpleFireSpread(fuelload, fueldepth, windspeed, slope, fuelmoisture, fuelsav):
    # Parameters
    maxval= 0
    if fuelload > 0:
        wo = fuelload # Ovendry fuel loading in (lb/ft^2). Amount of primary prod???
        fd = fueldepth # Fuel depth (ft)
        if fueldepth == 0:
            fd = 0.0000000000000001

        wv = windspeed * 88# Wind velocity at midflame height (ft/minute) = 88 * mph
        fpsa = fuelsav  # Fuel Particle surface area to volume ratio (1/ft)
        mf = fuelmoisture  # Fuel particle moisture content
        h = 8000  # Fuel particle low heat content
        pp = 32.  # Ovendry particle density
        st = 0.0555  # Fuel particle mineral content
        se = 0.010  # Fuel Particle effective mineral content
        mois_ext = 0.12  # Moisture content of extinction or 0.3 if dead
        #calculate slope as degrees
        slope_rad = math.atan(slope)
        slope_degrees = slope_rad / 0.0174533 #radians
        tan_slope = math.tan(slope_rad) #  in radians
        # Betas Packing ratio
        Beta_op = 3.348 * math.pow(fpsa, -0.8189)  # Optimum packing ratio
        ODBD = wo / fd  # Ovendry bulk density
        Beta = ODBD / pp #Packing ratio
        #Beta = 0.00158
        Beta_rel = Beta / Beta_op
        # Reaction Intensity
        WN = wo / (1 + st)  # Net fuel loading
        #A = 1 / (4.774 * pow(fpsa, 0.1) - 7.27)  # Unknown const
        A =  133.0 / math.pow(fpsa, 0.7913) #updated A
        T_max = math.pow(fpsa,1.5) * math.pow(495.0 + 0.0594 * math.pow(fpsa, 1.5),-1.0)  # Maximum reaction velocity
        #T_max = (fpsa*math.sqrt(fpsa)) / (495.0 + 0.0594 * fpsa * math.sqrt(fpsa))
        T = T_max * math.pow((Beta / Beta_op),A) * math.exp(A * (1 - Beta / Beta_op))  # Optimum reaction velocity
        # moisture dampning coefficient
        NM = 1. - 2.59 * (mf / mois_ext) + 5.11 * math.pow(mf / mois_ext, 2.) - 3.52 * math.pow(mf / mois_ext,3.)  # Moisture damping coeff.
        # mineral dampning
        NS = 0.174 * math.pow(se, -0.19)  # Mineral damping coefficient
        #print(T, WN, h, NM, NS)
        RI = T * WN * h * NM * NS
        #RI = 874
        # Propogating flux ratio
        PFR = mpmath.power(192.0 + 0.2595 * fpsa, -1) * mpmath.exp((0.792 + 0.681 * mpmath.sqrt(fpsa)) * (Beta + 0.1))# Propogating flux ratio
        ## Wind Coefficient
        B = 0.02526 * math.pow(fpsa, 0.54)
        C = 7.47 * math.exp(-0.1333 * math.pow(fpsa, 0.55))
        E = 0.715 * math.exp(-3.59 * 10**-4 * fpsa)
        #WC = C * wv**B * math.pow(Beta / Beta_op, -E) #wind coefficient
        if wv > (0.9 * RI): #important - don't know source. Matches BEHAVE
            wv = 0.9 * RI
        WC = (C * wv ** B) * math.pow((Beta / Beta_op), (-E))
        #WC= WC*0.74
        #Slope  coefficient
        SC = 5.275*(Beta**-0.3)*tan_slope**2
        #Heat sink

        EHN = math.exp(-138. / fpsa)  # Effective Heating Number = f(surface are volume ratio)
        QIG = 250. + 1116. * mf  # Heat of preignition= f(moisture content)
        # rate of spread (ft per minute)
        #RI = BTU/ft^2
        numerator = (RI * PFR * (1 + WC + SC))
        denominator = (ODBD * EHN * QIG)
        R = numerator / denominator #WC and SC will be zero at slope = wind = 0
        RT = 384.0/fpsa
        HA = RI*RT
        #fireline intensity as described by Albini via USDA Forest Service RMRS-GTR-371. 2018
        FI = (384.0/fpsa)*RI*(R) ##Uses Reaction Intensity in BTU / ft/ min
        #FI = HA*R
        if (RI <= 0):
            return (maxval, maxval, maxval)
        return (R, RI, FI)
    else:
        return (maxval, maxval, maxval)


@app.route('/kauil')
def index():
    return render_template('index.html')

@app.route('/kauil/load_wildfires')
def load_wildfires():
    # return 'TODO VERDE'
    try:
        # Carga del shapefile
        gdf = load_shapefile('Primavera.shp')
        heat_points = WildfireHotspots()
        # heat_points.append(WildfireHotspots())
        # heat_points = [(-103.47079999999988,20.598000000000006)]  # Agrega las coordenadas de los puntos de calor aquí
        # Download a regional GeoJSON of hotspots detected by the MODIS satellite in a recent 7-day period.
        # Devuelve los datos de coordenadas como JSON
        print(heat_points)
        return jsonify(heat_points)

    except Exception as e:
        # Captura cualquier excepción y registra el error
        app.logger.error(f"Error en la carga de incendios: {str(e)}")


async def fetch_weather_data(session, url, params):
    async with session.get(url, params=params) as response:
        if response.status == 200:
            return await response.json()
        else:
            return None

@app.route('/kauil/load_DataAndShowSpreadRate')
async def load_DataAndShowSpreadRate():
    intersected = get_intersected_data()
    feature_collection = await getWeatherDataAsync(intersected)
    # print(feature_collection)
        # Define el nombre del archivo JSON
    # json_filename = 'feature_collection.json'
    
    # # Lee el contenido del archivo JSON
    # with open(json_filename, 'r') as json_file:
    #     feature_collection = json.load(json_file)
    return jsonify(feature_collection)

async def getWeatherDataAsync(intersected):
    async with aiohttp.ClientSession() as session:
        ee.Initialize()
        tasks = []
        coordenadas = []
        fuel_depth=[]
        fuelload = []
            # Establecer los parámetros de la consulta
        fecha_inicio = '2023-01-01'
        fecha_fin = '2023-12-31'
         # Consulta de datos de vegetación
        coleccion = ee.ImageCollection('MODIS/006/MOD13A2').filterDate(fecha_inicio, fecha_fin)

        # Load the Global Forest Change dataset
        dataset = ee.Image('UMD/hansen/global_forest_change_2019_v1_7')

        # Select the tree cover band
        treeCover = dataset.select('treecover2000')
        for i, row in intersected.iterrows():
            y_centrum, x_centrum = row.geometry.centroid.y, row.geometry.centroid.x
            coordenadas.append([x_centrum, y_centrum])

            region = ee.Geometry.Point(coordenadas[i])
            roi = ee.Geometry.Point(coordenadas[i])
            imagen = coleccion.select('NDVI').mean().clip(region)

            # Obtener una región de interés como una imagen en forma de píxeles
            datos_imagen = imagen.reduceRegion(reducer=ee.Reducer.mean(), geometry=region, scale=20)

            # Obtener el valor promedio del NDVI en la región de interés
            valor_ndvi = datos_imagen.get('NDVI')
            
            
            # Extraer el valor del NDVI
            valor_ndvi = ee.Number(valor_ndvi)

            # Estimar la masa de vegetación en lb/ft^2
            slope_ = 0.5  # Pendiente de la relación lineal
            intercept = 0.2  # Término de intercepción de la relación lineal

            vegetation_mass = valor_ndvi.multiply(slope_).add(intercept)            

            fuelload.append(vegetation_mass.getInfo())
            
            
            # Calculate the mean tree cover within the region of interest
            treeCover_mean = treeCover.reduceRegion(ee.Reducer.mean(), roi, 30)

            # Get the tree cover value
            treeCoverValue = ee.Number(treeCover_mean.get('treecover2000'))

            # Convert tree cover from percentage to feet
            conversion_factor = 43.560  # 1 acre = 43,560 square feet
            treeCoverValue_feet = treeCoverValue.multiply(conversion_factor)
            

            fuel_depth.append(treeCoverValue_feet.getInfo())
            params = {
                "lat": y_centrum,
                "lon": x_centrum,
                "appid": "88ea1e1a48e480178b9c4ab177c2475f",
                "units": "metric"
            }

            url = "https://api.openweathermap.org/data/2.5/weather"
            task = fetch_weather_data(session, url, params)
            tasks.append(task)
        weather_data = await asyncio.gather(*tasks)
        wind_speed = []
        temp = []
        slope = []
        fuel_moisture = []
        for data in weather_data:
            if data is not None:
    # Extraer la velocidad del viento
                wind_speed.append(data["wind"]["speed"])
                temp.append(data["main"]["temp"])
                # return 'ENTRA'
    # SLOPE
                # Extraer la altitud y la presión atmosférica para calcular la pendiente de la respuesta JSON
                altitude = None
                atmospheric_pressure = data["main"]["pressure"]  # En hPa
                try:
                    altitude = data["main"]["sea_level"]  # En metros
                    slope.append(0.102 * altitude / atmospheric_pressure)  # En grados
                except KeyError:
                    slope.append(0)

    # FUEL MOISTURE
                # Calcular la humedad de combustible utilizando la fórmula de Deeming
                temp_calculate = data["main"]["temp"] - 273.15  # Convertir de Kelvin a Celsius
                relative_moisture = data["main"]["humidity"]
                precipitation = data.get("rain", {}).get("1h", 0)  # Obtener la precipitación de los últimos 60 minutos (si está disponible)
                wind_speed_calculate = data["wind"]["speed"]

                fuel_moisture.append(20 + 280 / (temp_calculate + 273) - relative_moisture + 0.1 * wind_speed_calculate + 0.2 * precipitation)
            else:
                print('La solicitud falló con el código de estado:')
        fuel_depth = replace_zero_with_average(fuel_depth)
        slope = replace_zero_with_average(slope)
        # Coordenadas de los puntos de calor de referencia
        # heat_points = [(-103.5306,20.7154), (-103.5487,20.7034), (-103.5388,20.7049), (-103.5413,20.6959), (-103.5314,20.6974) ,(-103.5217,20.6988) ,(-103.5298, 20.6883), (-103.52, 20.6898 )]  # Agrega las coordenadas de los puntos de calor aquí

        # Create a list of Features
        features = []
        for i, row in intersected.iterrows():
            fuelload_val = fuelload[i]
            fuel_depth_val = fuel_depth[i]
            fuelsav = 3500
            slope_val = slope[i] # Pendiente en grados
            wind_speed_val = wind_speed[i] * 2.237 # Velocidad del viento en m/s
            temperature_val = temp[i]
            fuel_moisture_val = fuel_moisture[i] / 1000 # Humedad de combustible
            heat_points = WildfireHotspots()
            # Calcular el índice de propagación del fuego para cada punto de calor de referencia
            point_indices = []
            if heat_points:
                for heat_point in heat_points:
                    heat_x, heat_y = heat_point
                    distance_to_heat = math.sqrt((row.geometry.centroid.x - heat_x)**2 + (row.geometry.centroid.y - heat_y)**2)
                    # Calcular la influencia basada en la intensidad del fuego en los heat points (puedes ajustar el factor)
                    fire_intensity = 100  # Aquí puedes usar la intensidad real del fuego en el punto de calor
                    if distance_to_heat != 0:
                        influence_factor = fire_intensity / distance_to_heat
                    else:
                        pass
                    point_indices.append(GetSimpleFireSpread(fuelload_val, fuel_depth_val, wind_speed_val, slope_val, fuel_moisture_val, fuelsav)[0] * influence_factor)
            else:
                point_indices.append(GetSimpleFireSpread(fuelload_val, fuel_depth_val, wind_speed_val, slope_val, fuel_moisture_val, fuelsav)[0] * 1)
            # Calcular el índice de propagación del fuego para el cuadrante actual
            I_val = sum(point_indices)
            I_val = float(I_val)
            # Create a GeoJSON Feature for each point
            feature = {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [row.geometry.centroid.x, row.geometry.centroid.y]
                },
                "properties": {
                    "I": I_val  # Replace with your calculated 'I' value
                }
            }
            features.append(feature)

    feature_collection = {
                "type": "FeatureCollection",
                "features": features
            }
    return feature_collection
    

if __name__ == '__main__':
    # http_server = WSGIServer(('', 5000), app)
    # %http_server.serve_forever()
    app.run(debug=True, host='0.0.0.0')
