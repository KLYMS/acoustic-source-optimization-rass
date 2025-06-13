import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
import numdifftools as nd
from concurrent.futures import ProcessPoolExecutor, as_completed

# Define the chunk size
chunk_size = 100000

def process_chunk(chunk):
    group_t = chunk.groupby('time')

    def surface_function(x, y):
        return rbf(x, y)

    def partial_derivative(func, var=0, point=[]):
        args = point[:]
        def wraps(x):
            args[var] = x
            return func(*args)
        return nd.Derivative(wraps)(point[var])

    results = []
    for t in chunk['time'].unique():
        if t != 0:
            df_ = group_t.get_group(t).to_numpy()
            unique_phi = np.unique(df_[:, 5])
            unique_theta = np.unique(df_[:, 3])
            
            for phi in unique_phi:
                for theta in unique_theta:
                    phi_idx = np.where(unique_phi == phi)[0][0]
                    theta_idx = np.where(unique_theta == theta)[0][0]
    
                    phi_values = unique_phi[max(0, phi_idx-1):min(phi_idx+2, len(unique_phi))]
                    theta_values = unique_theta[max(0, theta_idx-1):min(theta_idx+2, len(unique_theta))]
    
                    phi_mask = np.isin(df_[:, 5], phi_values)
                    theta_mask = np.isin(df_[:, 3], theta_values)
                    mask = phi_mask & theta_mask
                    
                    points = df_[mask]
                    
                    if points.shape[0] < 1:
                        continue
    
                    x = points[:, 0]
                    y = points[:, 1]
                    z = points[:, 2]
    
                    if len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
                        # Skip if x or y do not have enough unique values to create a surface
                        continue
    
                    rbf = Rbf(x, y, z, function='thin_plate')
    
                    specific_x, specific_y = 2.0, 3.0  
    
                    dz_dx = partial_derivative(surface_function, 0, [specific_x, specific_y])
                    dz_dy = partial_derivative(surface_function, 1, [specific_x, specific_y])
    
                    # Normal vector
                    normal = np.array([-dz_dx, -dz_dy, 1])
                    normal /= np.linalg.norm(normal)
    
                    results.append((specific_x, specific_y, normal))
    
    return results

if __name__ == '__main__':
    # Process data in chunks
    results = []
    with ProcessPoolExecutor() as executor:
        futures = []
        for chunk in pd.read_csv('/home/murali/Documents/rass/data/sim_data/wTemp_wWind_0_cpu_log.csv', chunksize=chunk_size):
            futures.append(executor.submit(process_chunk, chunk))
        
        for future in as_completed(futures):
            result = future.result()
            results.extend(result)

    # Print results
    for specific_x, specific_y, normal in results:
        print("Normal vector at ({}, {}): {}".format(specific_x, specific_y, normal))
