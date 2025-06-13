import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline, griddata
from concurrent.futures import ProcessPoolExecutor
import logging
import time

# Set up logging
logging.basicConfig(filename='backscattring_progress.log',level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Function to handle wrap-around index
def wrap_index(idx, total_length):
    return idx % total_length

def compute_normal_vector(interp_surface, x, y):
    dz_dx = interp_surface(x, y, dx=1, dy=0)[0, 0]
    dz_dy = interp_surface(x, y, dx=0, dy=1)[0, 0]
    normal_vector = np.array([-dz_dx, -dz_dy, 1.0])
    normal_vector /= np.linalg.norm(normal_vector)
    return normal_vector

def process_time_step(t_step, grid_size, total_steps, step_idx):
    results = []
    logging.info(f"Processing time step: {t_step} ({step_idx + 1}/{total_steps})")
    df = pd.read_csv(f'/home/murali/Documents/rass/data/sim_data/dived_data/wTemp_wWind_0_{t_step}.csv')
    total_operations = 0
    
    for idx, row in df.iterrows():
        phi = row['phi']
        theta = row['theta']

        points = [{
            'phi': phi,
            'theta': theta,
            'x': row['x'],
            'y': row['y'],
            'z': row['z']
        }]

        for d_phi in [-1, 0, 1]:
            for d_theta in [-1, 0, 1]:
                if d_phi == 0 and d_theta == 0:
                    continue

                idx_phi = wrap_index(idx + d_phi, len(df))
                idx_theta = wrap_index(idx + d_theta, len(df))

                points.append({
                    'phi': df.loc[idx_phi, 'phi'],
                    'theta': df.loc[idx_theta, 'theta'],
                    'x': df.loc[idx_phi, 'x'],
                    'y': df.loc[idx_theta, 'y'],
                    'z': df.loc[idx_theta, 'z']
                })

        if len(points) < 9:
            logging.info(f'Not enough points found for phi={phi} and theta={theta} ------ num of pts: {len(points)}')
            continue

        points_df = pd.DataFrame(points)

        grid_x, grid_y = np.meshgrid(
            np.linspace(points_df['x'].min(), points_df['x'].max(), grid_size),
            np.linspace(points_df['y'].min(), points_df['y'].max(), grid_size)
        )
        grid_z = griddata((points_df['x'], points_df['y']), points_df['z'], (grid_x, grid_y), method='cubic')

        interp_surface = RectBivariateSpline(
            np.linspace(points_df['x'].min(), points_df['x'].max(), grid_size),
            np.linspace(points_df['y'].min(), points_df['y'].max(), grid_size),
            grid_z
        )

        center_x, center_y , center_z = points_df.loc[0, 'x'], points_df.loc[0, 'y'] , points_df.loc[0,'z']
        normal_vector = compute_normal_vector(interp_surface, center_x, center_y)
        nx, ny = normal_vector[:2]

        t = -interp_surface(center_x, center_y)[0, 0] / (
            nx * interp_surface(center_x, center_y, dx=1, dy=0)[0, 0] + ny * interp_surface(center_x, center_y, dx=0, dy=1)[0, 0])

        x_zero = center_x + t * nx
        y_zero = center_y + t * ny

        if np.abs(x_zero) < 130 and np.abs(y_zero) < 130:
            z_zero = interp_surface(x_zero, y_zero)[0, 0]
            results.append({'x_val': center_x, 'y_val': center_y, 'z_val': center_z, 'x_land': x_zero, 'y_land': y_zero, 'time_step': t_step, 'phi': phi, 'theta': theta})
            logging.info(f"time_step: {t_step}, Coordinates (x, y) where backscatter at z = 0: ({x_zero}, {y_zero})")

        total_operations += 1
    
    logging.info(f"Finished processing time step: {t_step}, Total operations: {total_operations}")
    return results

def main():
    
    grid_size = 1000
    time_steps = np.loadtxt('/home/murali/Documents/rass/data/sim_data/dived_data/time_step.txt')
    all_results = []

    total_steps = len(time_steps)

    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_time_step, t_step, grid_size, total_steps, idx): t_step for idx, t_step in enumerate(time_steps) if t_step != 0}
        
        for future in futures:
            t_step = futures[future]
            try:
                result = future.result()
                all_results.extend(result)
                logging.info(f"Completed time step: {t_step}, Results obtained: {len(result)}")
            except Exception as e:
                logging.error(f"Error processing time step {t_step}: {e}")

    results_df = pd.DataFrame(all_results)
    results_df.to_csv('/home/murali/Documents/rass/simulation/in_py/backscatter_coordinates.csv', index=False)
    logging.info("Finished writing results to CSV")

if __name__ == "__main__":
    main()
