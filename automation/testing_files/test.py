import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from sklearn.cluster import DBSCAN
import time
import argparse
import os
import re

def extract_file_info(wf_file_path):
    file_name_pattern = r'/([^/]+)\_wf.csv$'
    date_pattern = r'F(\d{4})(\d{2})(\d{2})'
    
    file_name_match = re.search(file_name_pattern, wf_file_path)
    date_match = re.search(date_pattern, wf_file_path)
    
    file_name = file_name_match.group(1) if file_name_match else "No match found"
    year, month, day = map(int, date_match.groups()) if date_match else (None, None, None)
    
    return file_name, year, month, day

def calculate_land_points(df_, dist_axis, altitude_axis):
    land_pts = []
    wave_eq = CubicSpline(df_[dist_axis], df_[altitude_axis])
    derivative = wave_eq.derivative()
    
    for i in range(df_.shape[0]):
        x0 = df_.loc[i, dist_axis]
        dz_dx = derivative(x0)
        
        if dz_dx != 0:
            n_slope = np.rad2deg(np.arctan(-1 / dz_dx))
            if np.abs(n_slope) >= 70:
                land_x = x0 + wave_eq(x0) * dz_dx
                land_pts.append(land_x)
    return land_pts

def calculate_wave_points(dfa, dist_axis, altitude_axis, center):
    points_l, points_w = [], []
    dfa_centered = dfa.copy()
    dfa_centered[dist_axis] -= center
    group_t = dfa_centered.groupby('time')
    
    for time_slice, df_ in group_t:
        if time_slice != 0:
            df_ = df_.sort_values(dist_axis).drop_duplicates(subset=dist_axis).reset_index(drop=True)
            wave_eq = CubicSpline(df_[dist_axis], df_[altitude_axis])
            derivative = wave_eq.derivative()
            
            for x0 in df_[dist_axis]:
                dz_dx = derivative(x0)
                if dz_dx != 0:
                    n_slope = np.rad2deg(np.arctan(-1 / dz_dx))
                    if np.abs(n_slope) >= 70:
                        land_x = x0 + wave_eq(x0) * dz_dx
                        altitude = df_.loc[df_[dist_axis] == x0, altitude_axis].values[0]
                        points_w.append((x0, altitude))
                        
                        if np.abs(land_x) < 65:
                            points_l.append((x0, altitude))

    return np.unique(points_w, axis=0), np.unique(points_l, axis=0)

def find_nearest_location(locations, current_loc):
    distances = np.abs(np.array(locations) - current_loc)
    return locations[np.argmin(distances)]

def count_points_in_window(start, end, sorted_array):
    return np.sum((sorted_array >= start) & (sorted_array <= end))

def identify_top_windows(points, window_size, top_no_of_windows):
    top_windows, midpt_windows = [], []
    sorted_points = np.sort(points)
    
    for start_point in sorted_points:
        end_point = start_point + window_size
        point_count = count_points_in_window(start_point, end_point, sorted_points)
        
        if len(top_windows) < top_no_of_windows or point_count > top_windows[-1][2]:
            top_windows.append((start_point, end_point, point_count))
            midpt_windows.append((start_point + end_point) / 2)

    sorted_top_windows = sorted(zip(top_windows, midpt_windows), key=lambda x: x[0][2], reverse=True)
    return sorted_top_windows[:top_no_of_windows]

def calculate_angle(x, y):
    return np.degrees(np.arctan2(y, x))

def filter_points_in_cone(points, angle_start, angle_width):
    angle_end = angle_start + angle_width
    return [(x, y) for x, y in points if angle_start <= calculate_angle(x, y) <= angle_end]

def generate_angle_range(start, stop, step):
    while start < stop:
        yield start
        start += step

def main(wf_file_path, dist_axis):
    start_time = time.time()
    
    file_name, year, month, day = extract_file_info(wf_file_path)
    print(f"File: {file_name}\nDate: {year}-{month}-{day}")
    
    dfa = pd.read_csv(wf_file_path)
    altitude_axis = 'z'
    
    if dist_axis == 'x':
        dfa_filtered = dfa[(dfa['phi'] == 0) & (dfa[altitude_axis] >= 0) & (dfa['theta'] <= np.pi/2) & (dfa['theta'] >= -np.pi/2)].reset_index(drop=True)
        title = 'N-S'
    else:
        dfa_filtered = dfa[(dfa['phi'] != 0) & (dfa[altitude_axis] >= 0) & (dfa['theta'] <= np.pi/2) & (dfa['theta'] >= -np.pi/2)].reset_index(drop=True)
        title = 'E-W'
    
    window_size, top_no_of_windows = 130, 25
    initial_loc, acoustic_source_loc = 0, []
    
    group_t = dfa_filtered.groupby('time')
    for time_slice, df_ in group_t:
        if time_slice != 0:
            df_ = df_.sort_values(dist_axis).drop_duplicates(subset=dist_axis).reset_index(drop=True)
            land_pts = calculate_land_points(df_, dist_axis, altitude_axis)
            top_windows, midpt_windows = identify_top_windows(land_pts, window_size, top_no_of_windows)
            
            if top_windows:
                nearest_loc = find_nearest_location(midpt_windows, initial_loc)
                acoustic_source_loc.append(nearest_loc)
    
    acoustic_source_loc = np.array(acoustic_source_loc).reshape(-1, 1)
    db = DBSCAN(eps=5, min_samples=3).fit(acoustic_source_loc)
    
    cluster_centers = [np.mean(acoustic_source_loc[db.labels_ == label]) for label in set(db.labels_) if label != -1]
    print(f'\nCluster centers: {cluster_centers}\n')
    
    angle_width, angle_start, angle_end, angle_step = 3, 68, 113 - angle_width, 1
    points_w, points_l = zip(*[calculate_wave_points(dfa_filtered, dist_axis, altitude_axis, center) for center in cluster_centers])
    
    best_angle, max_percent_pts = 0, 0
    for angle in generate_angle_range(angle_start, angle_end, angle_step):
        cone_points_w = [point for pw in points_w for point in filter_points_in_cone(pw, angle, angle_width)]
        cone_points_l = [point for pl in points_l for point in filter_points_in_cone(pl, angle, angle_width)]
        
        percent_pts = (len(cone_points_l) / len(cone_points_w)) * 100 if cone_points_w else 0
        if percent_pts > max_percent_pts:
            best_angle, max_percent_pts = angle, percent_pts
    
    print(f"Percentage inside the cone: {max_percent_pts}, angle between: {best_angle} - {best_angle + angle_width}")
    
    points_w_flat = np.vstack(points_w)
    points_l_flat = np.vstack(points_l)
    
    plt.scatter(points_w_flat[:, 0], points_w_flat[:, 1], s=0.9, label='Wave Points')
    plt.scatter(points_l_flat[:, 0], points_l_flat[:, 1], s=0.9, color='#A05544', label='Land Points')
    
    dist = 35 * 1000
    plt.xlim((-dist, dist))
    plt.ylim((0, dist))
    plt.xlabel("Distance from source (m)")
    plt.ylabel("Altitude (m)")
    plt.title(title)
    plt.grid()
    
    theta1, theta2 = np.radians(best_angle), np.radians(best_angle + angle_width)
    plt.plot([0, np.cos(theta1) * dist], [0, np.sin(theta1) * dist], 'g--', label='Cone-section')
    plt.plot([0, np.cos(theta2) * dist], [0, np.sin(theta2) * dist], 'g--')
    
    plt.legend()
    save_path = '/home/murali/Documents/rass/automation/gen_data/'
    plt.savefig(f"{save_path}{file_name}_{title}.png")
    
    acoustic_sources_file = f"{save_path}acoustic_sources.csv"
    if os.path.exists(acoustic_sources_file):
        df = pd.read_csv(acoustic_sources_file)
    else:
        df = pd.DataFrame(columns=['Cluster_Centers', 'Axis', 'date'])
    
    df_centers = pd.DataFrame(cluster_centers, columns=['Cluster_Centers'])
    df_centers['Axis'] = dist_axis
    df_centers['date'] = pd.Timestamp(f"{year}-{month}-{day}")
    df = pd.concat([df, df_centers])
    df.to_csv(f"{save_path}{file_name}_{title}.csv", index=False)
    
    print(f"Time taken: {time.time() - start_time}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('file_path', type=str, help='Path to the CSV file')
    parser.add_argument('dist_axis', type=str, help='Distance axis (x or y)')
    
    args = parser.parse_args()
    main(args.file_path, args.dist_axis)
