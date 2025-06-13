import numpy as np
import pandas as pd
import csv
from bisect import bisect_left
from typing import List

Gamma = 1.4
SpeedOfSound = 343
R_gas = 287
Troposphere_Lapse_Rate = -6.5
Tropopause_Height = 18000
Stratosphere_Lapse_Rate = 2.5
Stratopause_Height = 48000
Ref_Wind_speed = 250.5
Ref_height = 28000
With_temperature = True
With_wind = True

class ODESystem:
    def __init__(self, phi, height_data, wind_data, wind_direction_data, temp_data):
        self.phi = phi
        self.height_data = height_data
        self.wind_data = wind_data
        self.wind_direction_data = wind_direction_data
        self.temp_data = temp_data

    def __call__(self, t, y):
        x_coord = y[0]
        y_coord = y[1]
        z_coord = y[2]
        theta = y[3]

        dx_dt = self.v_x(z_coord) + self.c_s(z_coord) * np.sin(theta) * np.cos(self.phi)
        dy_dt = self.v_y(z_coord) + self.c_s(z_coord) * np.sin(theta) * np.sin(self.phi)
        dz_dt = self.v_z(z_coord) + self.c_s(z_coord) * np.cos(theta)
        dtheta_dt = np.sin(theta) * (self.partial_cs_z(z_coord) + np.cos(self.phi) * self.partial_vx_z(z_coord) + np.sin(self.phi) * self.partial_vy_z(z_coord))

        return [dx_dt, dy_dt, dz_dt, dtheta_dt]

    def temperature(self, z):
        if With_temperature:
            if z <= Tropopause_Height:
                return 300 + (Troposphere_Lapse_Rate / 1000) * z
            elif z > Tropopause_Height and z <= Stratopause_Height:
                return self.temperature(Tropopause_Height) + (Stratosphere_Lapse_Rate / 1000) * (z - Tropopause_Height)
            else:
                return self.temperature(Stratopause_Height)
        else:
            return 0

    def wind_s_d(self,z):
        ws =  avg_data(self.height_data,self.wind_data,z)
        wd =  np.deg2rad(avg_data(self.height_data, self.wind_direction_data,z))
        return ws , wd

    def v_x(self, z):
        if With_wind:
            ws , wd = self.wind_s_d(z)
            return ws * np.cos(self.phi - wd)
        else:
            return 0
        # return 0

    def v_y(self, z):
        if With_wind:
            ws , wd = self.wind_s_d(z)
            return ws * np.sin(self.phi - wd)
        else:
            return 0

    def v_z(self, z):
        if With_wind:
            return 0
        else:
            return 0

    def c_s(self, z):
        if With_temperature:
            return np.sqrt(Gamma * R_gas * avg_data(self.height_data,self.temp_data,z))
        else:
            return SpeedOfSound

    def partial_cs_z(self, z):
        if With_temperature:
            return 0.5 * np.sqrt(Gamma * R_gas / avg_data(self.height_data,self.temp_data,z)) * slope(self.height_data,self.temp_data,z)
        else:
            return 0

    def partial_vx_z(self, z):
        if With_wind:
            wind_speed , wind_dir = self.wind_s_d(z)
            val = slope(self.height_data,self.wind_data,z) * np.cos(wind_dir) -(
                    np.sin(wind_dir) * slope(self.height_data,self.wind_data,z) * wind_speed
                    )
            return val
        else:
            return 0

    def partial_vy_z(self, z):
        if With_wind:
            wind_speed , wind_dir = self.wind_s_d(z)
            val = slope(self.height_data,self.wind_data,z) * np.sin(wind_dir) +(
                    np.cos(wind_dir) * slope(self.height_data,self.wind_data,z) * wind_speed
                    )
            return val
        else:
            return 0


class RK4Solver:
    def __init__(self, system, dt):
        self.system = system
        self.dt = dt

    def step(self, t, y):
        k1 = np.array(self.system(t, y))
        k2 = np.array(self.system(t + self.dt / 2, y + (k1 * self.dt) / 2))
        k3 = np.array(self.system(t + self.dt / 2, y + (k2 * self.dt) / 2))
        k4 = np.array(self.system(t + self.dt, y + k3 * self.dt))
        return y + (self.dt/6) * (k1 + 2 * k2 + 2 * k3 + k4)


def avg_data( height_data: List[float], data: List[float], target_height: float) -> float:
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


def integrate_(y,km,h, phi, height_data, wind_data, wind_direction_data, temp_data):

    t0 = 0
    t_end = (km * 1000) / SpeedOfSound
    steps = int((t_end - t0)/ h)
    
    ode_system = ODESystem(phi, height_data, wind_data, wind_direction_data, temp_data)
    rk4_solver = RK4Solver(ode_system, h)

    y_values = np.zeros((steps + 1, len(y)))
    y_values[0] = y
    t_values = np.zeros(steps + 1)
    t_values[0] = t0
    for i in range(1,steps+1):
        t_values[i] = t_values[i-1] + h
        y_values[i] = rk4_solver.step(t_values[i-1],y_values[i-1])

    return y_values , t_values


def drange(x, y, jump):
  while x < y:
    yield x
    x += jump

km = 30
h = 2
theta_min = - np.pi/2
theta_max = np.pi/2
theta_step = np.pi/500

# phi_min = 0 # 0 deg
# phi_max = np.pi # 180 deg
# phi_step = np.pi

phi = np.deg2rad(0)

data = pd.read_csv('/home/murali/Documents/rass/data/radiosonde_grouped_data.csv')
height_data = data['Height'].to_numpy()
wind_data = data['WS'].to_numpy()
wind_direction_data = data['WD'].to_numpy()
temp_data = data['Temp'].to_numpy()

fields = ['x', 'y', 'z', 'theta' , 'time']
w_filename = ''
if(With_temperature and With_wind):
    w_filename = "wTemp_wWind_data.csv"
elif(With_temperature):
    w_filename = "wTemp_data.csv"
elif(With_wind):
    w_filename = "wWind_data.csv"
else:
    w_filename = "ideal_case_data.csv"

with open(w_filename, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fields)
    
    for theta in drange(theta_min,theta_max,theta_step):
        y = [0,0,0,theta]
        y_values , t_vlaues = integrate_(y,km,h,phi, height_data, wind_data, wind_direction_data, temp_data)
        values = np.column_stack([y_values,t_vlaues])
        writer.writerows(values)






