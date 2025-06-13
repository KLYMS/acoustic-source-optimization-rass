import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline, griddata
import csv
import time

start = time.time()

SpeedofSound = 343

def wrap_index(idx, total_length):
    return idx % total_length

def compute_normal_vector(interp_surface, x, y):
    dz_dx = interp_surface(x, y, dx=1, dy=0)[0, 0]
    dz_dy = interp_surface(x, y, dx=0, dy=1)[0, 0]
    normal_vector = np.array([-dz_dx, -dz_dy, 1.0])
    normal_vector /= np.linalg.norm(normal_vector)
    return normal_vector

grid_size = 5

results_file = '/home/murali/Documents/rass/data/sim_data/dived_data/py_data/alglib_CS_i_0/backscatter_coordinates_i_0_d_G5.csv'
time_steps = np.loadtxt('/home/murali/Documents/rass/data/sim_data/dived_data/py_data/alglib_CS_i_0/time_step.txt')

with open(results_file, mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=['x_val', 'y_val','z_val', 'x_land', 'y_land', 'time_step'])
    writer.writeheader()

    for idx, t_step in enumerate(time_steps):
        if t_step != 0:
            print(f"time step: {t_step}, {idx + 1}/{len(time_steps)}")
            df = pd.read_csv(f'/home/murali/Documents/rass/data/sim_data/dived_data/py_data/alglib_CS_i_0/wTemp_wWind_0_{t_step}.csv')
            
            df = df[(df['z'] >= 0) & (df['x'] < SpeedofSound * t_step) & (df['y'] < SpeedofSound * t_step)].reset_index(drop=True)
        
            for row_idx, row in df.iterrows():
                phi = row['phi']
                theta = row['theta']
                
                points = []
                
                points.append({
                    'phi': phi,
                    'theta': theta,
                    'x': row['x'],
                    'y': row['y'],
                    'z': row['z']
                })
            
                for d_phi in [-1, 0, 1]:
                    for d_theta in [-1, 0, 1]:
                        if d_phi == 0 and d_theta == 0:
                            continue 
            
                        idx_phi = wrap_index(row_idx + d_phi, len(df))
                        idx_theta = wrap_index(row_idx + d_theta, len(df))
            
                        points.append({
                            'phi': df.loc[idx_phi, 'phi'],
                            'theta': df.loc[idx_theta, 'theta'],
                            'x': df.loc[idx_phi, 'x'],
                            'y': df.loc[idx_theta, 'y'],
                            'z': df.loc[idx_theta, 'z']
                        })
                        
                if len(points) < 9:
                    print(f"Not enough points found for phi={phi} and theta={theta}")
                    continue
            
                points_df = pd.DataFrame(points)
                
                grid_x, grid_y = np.meshgrid(np.linspace(points_df['x'].min(), points_df['x'].max(), grid_size), 
                                            np.linspace(points_df['y'].min(), points_df['y'].max(), grid_size))
                grid_z = griddata((points_df['x'], points_df['y']), points_df['z'], (grid_x, grid_y), method='cubic')
            
                interp_surface = RectBivariateSpline(np.linspace(points_df['x'].min(), points_df['x'].max(), grid_size), 
                                                    np.linspace(points_df['y'].min(), points_df['y'].max(), grid_size), 
                                                    grid_z)
                
                center_x, center_y, center_z = points_df.loc[0, 'x'], points_df.loc[0, 'y'], points_df.loc[0, 'z']
                normal_vector = compute_normal_vector(interp_surface, center_x, center_y)
                nx, ny = normal_vector[:2]
                
                t = -interp_surface(center_x, center_y)[0, 0] / (nx * interp_surface(center_x, center_y, dx=1, dy=0)[0, 0] + ny * interp_surface(center_x, center_y, dx=0, dy=1)[0, 0])
                
                x_zero = center_x + t * nx
                y_zero = center_y + t * ny
                
                if (np.abs(x_zero) < 130) & (np.abs(y_zero) < 130):
                    result = {'x_val': center_x, 'y_val': center_y, 'z_val': center_z, 'x_land': x_zero, 'y_land': y_zero, 'time_step': t_step}
                    writer.writerow(result)
                    print(f"time_step: {t_step}, Coordinates (x, y) where backscatter at z = 0: ({x_zero}, {y_zero})")

end = time.time()
print(f"time taken: {end - start}")
