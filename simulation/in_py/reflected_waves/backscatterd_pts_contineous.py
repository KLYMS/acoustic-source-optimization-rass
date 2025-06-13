import numpy as np
import pandas as pd
import time
import csv
from scipy.interpolate import RectBivariateSpline, griddata

start = time.time()

def compute_normal_vector(interp_surface, x, y):
    dz_dx = interp_surface(x, y, dx=1, dy=0)[0, 0]
    dz_dy = interp_surface(x, y, dx=0, dy=1)[0, 0]
    normal_vector = np.array([-dz_dx, -dz_dy, 1.0])
    normal_vector /= np.linalg.norm(normal_vector)
    return normal_vector

# Load time steps
time_steps = np.loadtxt('/home/murali/Documents/rass/data/sim_data/dived_data/py_data/time_step.txt')

# Parameters
grid_size = 200
rounding_pic = 0

results_file = '/home/murali/Documents/rass/data/sim_data/dived_data/py_data/backscatter_coordinates_contineous.csv'

with open(results_file, mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=['x_val', 'y_val','z_val', 'x_land', 'y_land', 'time_step'])
    writer.writeheader()
    
    for idx, t_step in enumerate(time_steps):
        if t_step != 0:
            print(f"time step: {t_step}, {idx}/{len(time_steps)}")
            df = pd.read_csv(f'/home/murali/Documents/rass/data/sim_data/dived_data/py_data/wTemp_wWind_0_{t_step}.csv')
            
            df_ = df.to_numpy()
            df_[:, 5] = np.round(df_[:, 5], decimals=rounding_pic)
            df_[:, 3] = np.round(df_[:, 3], decimals=rounding_pic)
            
            unique_phi = np.unique(df_[:, 5])
            unique_theta = np.unique(df_[:, 3])
            
            for phi in unique_phi:
                for theta in unique_theta:
                    
                    print('Phi values:', phi)
                    print('Theta values:', theta)
                    
                    mask = (df_[:, 3] == theta) & (df_[:, 5] == phi)
                    points = df_[mask]
                    
                    if len(points) == 0:
                        continue
                    
                    x = points[:, 0]
                    y = points[:, 1]
                    z = points[:, 2]
                    
                    print(f"x: {x}, y: {y}, z: {z}")
                    
                    grid_x, grid_y = np.meshgrid(
                        np.linspace(np.min(df_[:, 0]), np.max(df_[:, 0]), grid_size),
                        np.linspace(np.min(df_[:, 1]), np.max(df_[:, 1]), grid_size)
                    )
                    grid_z = griddata((df_[:, 0], df_[:, 1]), df_[:, 2], (grid_x, grid_y), method='cubic')
                    
                    interp_surface = RectBivariateSpline(
                        np.linspace(np.min(df_[:, 0]), np.max(df_[:, 0]), grid_size),
                        np.linspace(np.min(df_[:, 1]), np.max(df_[:, 1]), grid_size),
                        grid_z
                    )
                    
                    for x_val, y_val, z_val in zip(x, y, z):
                        normal_vector = compute_normal_vector(interp_surface, x_val, y_val)
                        nx, ny = normal_vector[:2]
                        
                        t = -interp_surface(x_val, y_val)[0, 0] / (nx * interp_surface(x_val, y_val, dx=1, dy=0)[0, 0] + ny * interp_surface(x_val, y_val, dx=0, dy=1)[0, 0])
                        
                        x_zero = x_val + t * nx
                        y_zero = y_val + t * ny
                        
                        if (np.abs(x_zero) < 130) & (np.abs(y_zero) < 130):
                            result = {
                                'x_val': x_val,
                                'y_val': y_val,
                                'z_val': z_val,
                                'x_land': x_zero,
                                'y_land': y_zero,
                                'time_step': t_step
                            }
                            writer.writerow(result)
                            print(f"time_step: {t_step}, Coordinates (x, y) where backscatter at z = 0: ({x_zero}, {y_zero})")

end = time.time()
print(f"time taken: {end-start}")
