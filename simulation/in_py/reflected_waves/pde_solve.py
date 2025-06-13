import numpy as np
import pandas as pd
import csv
from scipy.interpolate import CubicSpline


Gamma = 1.4
SpeedOfSound = 343
R_air = 287
Troposphere_Lapse_Rate = -6.5
Tropopause_Height = 18000
Stratosphere_Lapse_Rate = 2.5
Stratopause_Height = 48000

With_temperature = True
With_wind = True
Ideal_Conditions = False


def Temperature_(z):
    if With_temperature:
        if z <= Tropopause_Height:
            return 300 + (Troposphere_Lapse_Rate / 1000) * z
        elif z > Tropopause_Height and z <= Stratopause_Height:
            return Temperature_(Tropopause_Height) + (Stratosphere_Lapse_Rate / 1000) * (z - Tropopause_Height)
        else:
            return Temperature_(Stratopause_Height)
    else:
        return 0

def Speed_Of_Sound(temperature):
    return np.sqrt(Gamma* R_air * temperature)

def ray_tracing_eqns(y_, phi):
    x, y, z, theta = y_

    if not Ideal_Conditions:
        
        if With_temperature:
            temp = temperature_spline(z)
            dc_dz = np.sqrt(Gamma*R_air/temp) * 0.5 * temperature_spline.derivative()(z)
        else:
            temp = 292.804
            dc_dz = 0
            
        if With_wind:   
            wind_speed = wind_speed_spline(z)
            wind_dir = wind_direction_spline(z)
            dvx_dz = wind_speed_spline.derivative()(z) * np.cos(np.deg2rad(wind_dir)) -(
                    np.sin(np.deg2rad(wind_dir)) * wind_direction_spline.derivative()(z) * wind_speed
            ) 
            dvy_dz = wind_speed_spline.derivative()(z) * np.sin(np.deg2rad(wind_dir)) + (
                    np.cos(np.deg2rad(wind_dir)) * wind_direction_spline.derivative()(z) * wind_speed
            )
        else:
            wind_speed = 0
            wind_dir = 0
            dvx_dz = 0
            dvy_dz = 0
            
    else:
        if With_temperature:
            temp = Temperature_(z)
            
            lapse = 0
            if (z<= Tropopause_Height):
                lapse = Troposphere_Lapse_Rate
            elif (z >= Tropopause_Height and z <= Stratopause_Height):
                lapse = Stratosphere_Lapse_Rate
                
            dc_dz =  np.sqrt(Gamma*R_air/temp) * 0.5 * (lapse/1000)
            
        else:
            temp = 292.804
            dc_dz = 0
        
        wind_speed = 0
        wind_dir = 0
        dvx_dz = 0
        dvy_dz = 0

    c_s = Speed_Of_Sound(temp)
    v_x = wind_speed * np.cos(np.deg2rad(wind_dir))
    v_y = wind_speed * np.sin(np.deg2rad(wind_dir))
    

    dx_dt = v_x + c_s * np.sin(theta) * np.cos(phi)
    dy_dt = v_y + c_s * np.sin(theta) * np.sin(phi)
    dz_dt = c_s * np.cos(theta)
    
    dtheta_dt = np.sin(theta) * (
                dc_dz + np.sin(theta) * (np.cos(phi) * dvx_dz + np.sin(phi) * dvy_dz)
        )
    
    return np.array([dx_dt, dy_dt, dz_dt, dtheta_dt])

def Rk4Solver(y, phi, h):
    
    k1 = h * ray_tracing_eqns(y,phi)
     
    y_temp = y + k1/2
    k2 = h * ray_tracing_eqns(y_temp,phi)
    
    y_temp = y + k2/2
    k3 = h * ray_tracing_eqns(y_temp,phi)
    
    y_temp = y + k3
    k4 = h * ray_tracing_eqns(y_temp,phi)
    
    y_value = y + (h/6) * (k1 + 2*k2 + 2 * k3 + k4)

    return y_value

def integrate_(y,km,h, phi):

    t0 = 0
    t_end = (km * 1000) / SpeedOfSound
    steps = int((t_end - t0)/ h)

    y_values = np.zeros((steps + 1, len(y)))
    y_values[0] = y
    t_values = np.zeros(steps + 1)
    t_values[0] = t0
    for i in range(1,steps+1):
        t_values[i] = t_values[i-1] + h
        y_values[i] = Rk4Solver(y_values[i-1], phi, h)

    return y_values , t_values

# def integrate_(y, km, h, phi): # pts do not go beyond specified km
#     distance_limit = km * 1000 
#     total_distance = 0

#     y_values = [y]
#     t_values = [0]
    
#     while total_distance < distance_limit:
#         y_next = Rk4Solver(y_values[-1], phi, h)
#         y_values.append(y_next)
#         t_values.append(t_values[-1] + h)

#         dx = y_next[0] - y_values[-2][0]
#         dy = y_next[1] - y_values[-2][1]
#         dz = y_next[2] - y_values[-2][2]
#         step_distance = np.sqrt(dx**2 + dy**2 + dz**2)
#         total_distance += step_distance

#     y_values = np.array(y_values)
#     t_values = np.array(t_values)

#     return y_values, t_values


def drange(x, y, jump):
  while x < y:
    yield x
    x += jump

if not Ideal_Conditions:
    print('Reading the Radiosonde data')
    data = pd.read_csv('/home/murali/Documents/rass/data/radiosonde_grouped_data.csv')

    height = data['Height'].values
    temperature = data['Temp'].values
    wind_speed = data['WS'].values
    wind_direction = data['WD'].values

    print("Creating Splines from the data")
    temperature_spline = CubicSpline(height, temperature,bc_type='not-a-knot')
    wind_speed_spline = CubicSpline(height, wind_speed)
    wind_direction_spline = CubicSpline(height, wind_direction)

km = 30
h = 1.5
theta_min = - np.pi/2
theta_max = np.pi/2
theta_step = np.pi/350

phi = np.deg2rad(0)

fields = ['x', 'y', 'z', 'theta', 'time']
w_filename = ''

if(With_temperature and With_wind):
    w_filename = "wTemp_wWind.csv"
elif(With_temperature):
    w_filename = "wTemp.csv"
elif(With_wind):
    w_filename = "wWind.csv"
else:
    w_filename = "ideal_case.csv"

with open(w_filename, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fields)

    for theta in drange(theta_min,theta_max,theta_step):
        y = [0,0,0,theta]
        y_values , t_values = integrate_(y,km,h,phi)
        values = np.column_stack([y_values,t_values])
        writer.writerows(values)




