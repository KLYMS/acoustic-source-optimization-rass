import numpy as np
import pandas as pd
import csv
from scipy.interpolate import CubicSpline
from bisect import bisect_left
from typing import List
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

# Constants
Gamma = 1.4
SpeedOfSound = 343
R_air = 287
Troposphere_Lapse_Rate = -6.5
Tropopause_Height = 18000
Stratosphere_Lapse_Rate = 2.5
Stratopause_Height = 48000

# Flags
Ideal_Conditions = True
With_temperature = False
With_wind = False
With_spline = False

def Temperature_(z):
    if With_temperature:
        if z <= Tropopause_Height:
            return 300 + (Troposphere_Lapse_Rate / 1000) * z
        elif z > Tropopause_Height and z <= Stratopause_Height:
            return Temperature_(Tropopause_Height) + (Stratosphere_Lapse_Rate / 1000) * (z - Tropopause_Height)
        else:
            return Temperature_(Stratopause_Height)
    else:
        return 292.804

def Speed_Of_Sound(temperature):
    if With_temperature:
        return np.sqrt(Gamma* R_air * temperature)
    else:
        return SpeedOfSound
    
def avg_data(height_data: List[float], data: List[float], target_height: float) -> float:
    if target_height in height_data:
        index = height_data.index(target_height)
        return data[index]
    else:
        insert_pos = bisect_left(height_data, target_height)
        if insert_pos == 0:
            return data[0]
        elif insert_pos == len(height_data):
            return data[-1]
        else:
            h1 = height_data[insert_pos - 1]
            h2 = height_data[insert_pos]
            t1 = data[insert_pos - 1]
            t2 = data[insert_pos]
            slope_ = (t2 - t1) / (h2 - h1)
            interpolated_data = t1 + ((target_height - h1) * slope_)
            return interpolated_data

def forward_calc(index: int, height_data: List[float], temp_data: List[float]) -> float:
    h1 = height_data[index]
    h2 = height_data[index + 1]
    t1 = temp_data[index]
    t2 = temp_data[index + 1]
    return (t2 - t1) / (h2 - h1)

def backward_calc(index: int, height_data: List[float], temp_data: List[float]) -> float:
    h1 = height_data[index - 1]
    h2 = height_data[index]
    t1 = temp_data[index - 1]
    t2 = temp_data[index]
    return (t2 - t1) / (h2 - h1)

def slope(height_data: List[float], data: List[float], target_height: float) -> float:
    if len(height_data) != len(data):
        raise ValueError("height_data and temp_data must have the same size")

    if target_height in height_data:
        index = height_data.index(target_height)
        if index < len(height_data) - 1:
            return forward_calc(index, height_data, data)
        elif index == len(height_data) - 1:
            return backward_calc(index, height_data, data)
        else:
            raise IndexError("Target height is more than the last element in height_data, slope calculation is not possible.")
    else:
        insert_pos = bisect_left(height_data, target_height)
        if insert_pos >= len(height_data) - 1:
            index = len(height_data) - 1
            return backward_calc(index, height_data, data)
        return backward_calc(insert_pos, height_data, data)

def ray_tracing_eqns(y_, phi):
    x, y, z, theta = y_
    if Ideal_Conditions:
        if With_temperature:
            temp = Temperature_(z)
            if temp <= 0:
                temp =  292.804

            lapse = 0
            if (z <= Tropopause_Height):
                lapse = Troposphere_Lapse_Rate
            elif (z >= Tropopause_Height and z <= Stratopause_Height):
                lapse = Stratosphere_Lapse_Rate

            dc_dz =  np.sqrt(Gamma * R_air / temp) * 0.5 * (lapse / 1000)

        else:
            temp = 292.804
            dc_dz = 0

        v_x = 0
        v_y = 0
        dvx_dz = 0
        dvy_dz = 0

    else:
        if With_temperature:
            if With_spline:
                temp = temperature_spline(z)
                if temp <= 0:
                    temp = 292.804
                dc_dz = np.sqrt(Gamma * R_air / temp) * 0.5 * temperature_spline.derivative()(z)
            else:
                temp = avg_data(height_data, temperature_data, z)
                if temp <= 0:
                    temp =  292.804
                dc_dz = np.sqrt(Gamma * R_air / temp) * 0.5 * slope(height_data, temperature_data, z)
        else:
            temp = 292.804
            dc_dz = 0

        if With_wind:
            if With_spline:
                v_x = wind_x_spline(z)
                v_y = wind_y_spline(z)
                dvx_dz = wind_x_spline.derivative()(z)
                dvy_dz = wind_y_spline.derivative()(z)
            else:
                v_x = avg_data(height_data, wind_x_data, z)
                v_y = avg_data(height_data, wind_y_data, z)
                dvx_dz = slope(height_data, wind_x_data, z)
                dvy_dz = slope(height_data, wind_y_data, z)
        else:
            v_x = 0
            v_y = 0
            dvx_dz = 0
            dvy_dz = 0

    c_s = Speed_Of_Sound(temp)

    # PDE of Acoustic Wave Propagation 
    dx_dt = v_x + c_s * np.sin(theta) * np.cos(phi)
    dy_dt = v_y + c_s * np.sin(theta) * np.sin(phi)
    dz_dt = c_s * np.cos(theta)

    dtheta_dt = np.sin(theta) * (
                dc_dz + np.sin(theta) * (np.cos(phi) * dvx_dz + np.sin(phi) * dvy_dz)
        )

    return np.array([dx_dt, dy_dt, dz_dt, dtheta_dt])

def Rk4Solver(y, phi, h):
    k1 =  ray_tracing_eqns(y, phi)

    y_temp = y + (k1 * h) / 2
    k2 =  ray_tracing_eqns(y_temp, phi)

    y_temp = y + (k2 * h) / 2
    k3 =  ray_tracing_eqns(y_temp, phi)

    y_temp = y + k3 * h
    k4 =  ray_tracing_eqns(y_temp, phi)

    y_value = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return y_value

def integrate_(y, km, h, phi):
    t0 = 0
    t_end = (km * 1000) / SpeedOfSound
    steps = int((t_end - t0) / h)

    y_values = np.zeros((steps + 1, len(y)))
    y_values[0] = y
    t_values = np.zeros(steps + 1)
    t_values[0] = t0
    phi_values = np.zeros(steps + 1)
    phi_values[0] = phi
    for i in range(1, steps + 1):
        t_values[i] = t_values[i-1] + h
        y_values[i] = Rk4Solver(y_values[i-1], phi, h)
        phi_values[i] = phi

    return y_values, t_values, phi_values

def drange(x, y, jump):
    while x < y:
        yield x
        x += jump

def process_combination(args):
    arr, phi, theta, km, h = args
    x_val = arr[0]
    y_val = arr[1]
    z_val = arr[2]
    y = [x_val, y_val, z_val, theta]
    y_values, t_values, phi_values = integrate_(y, km, h, phi)
    values = np.column_stack([y_values, t_values, phi_values])

    return values

# read the data
if not Ideal_Conditions:
    print('Reading the Radiosonde data')
    data = pd.read_csv('/home/murali/Documents/rass/data/radiosonde_grouped_data.csv')

    height_data = data['Height'].values
    temperature_data = data['Temp'].values
    wind_speed_data = data['WS'].values
    wind_direction_data = data['WD'].values
    wind_x_data = wind_speed_data * np.cos(np.deg2rad(wind_direction_data))
    wind_y_data = wind_speed_data * np.sin(np.deg2rad(wind_direction_data))
    
    if With_spline:
        print("Creating Splines from the data")
        temperature_spline = CubicSpline(height_data, temperature_data)
        wind_speed_spline = CubicSpline(height_data, wind_speed_data)
        wind_direction_spline = CubicSpline(height_data, wind_direction_data)
        wind_x_spline = CubicSpline(height_data, wind_x_data)
        wind_y_spline = CubicSpline(height_data, wind_y_data)
    
    km = height_data.max() / 1000
else:
    km = 30

h = 2

theta_min = -np.pi / 2
theta_max = np.pi / 2
theta_step = np.pi / 350

phi_min = np.deg2rad(0)
phi_max = np.deg2rad(180)
phi_step = np.deg2rad(1)

corner_val = np.array([
    [0, 0, 0],        # center
    # [-6500, 6500, 0],     # upper left
    # [65, 65, 0],      # upper right
    # [-65, -65, 0],    # lower left
    # [65, -65, 0]      # lower right
])

edge_val = np.array([
    [0, 0, 0],       # center
    [0, 65, 0],      # upper edge
    [65, 0, 0],      # right edge
    [0, -65, 0],     # lower edge
    [-65, 0, 0]      # left edge
])

fields = ['x', 'y', 'z', 'theta', 'time', 'phi']
w_filename = ''

if With_temperature and With_wind:
    w_filename = "wTemp_wWind_4.csv"
elif With_temperature and not Ideal_Conditions:
    w_filename = "wTemp.csv"
elif With_temperature and Ideal_Conditions:
    w_filename = "ideal_wTemp.csv"
else:
    w_filename = "ideal_case.csv"

start = time.time()

combinations = [(arr, phi, theta, km, h) for arr in corner_val for phi in drange(phi_min, phi_max, phi_step) for theta in drange(theta_min, theta_max, theta_step)]

with open(w_filename, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fields)
    
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_combination, comb) for comb in combinations]
        
        for future in as_completed(futures):
            result = future.result()
            writer.writerows(result)

end = time.time()
print(f"Time Taken {end - start}")
