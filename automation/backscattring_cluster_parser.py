import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from sklearn.cluster import DBSCAN
import time, os, re, argparse
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

def land_points(df_, dist_axis, altitude_axis):
    land_pts = []
    wave_eq = CubicSpline(df_[dist_axis], df_[altitude_axis])

    for i in range(df_.shape[0]):
        x0 = df_.loc[i, dist_axis]
        dz_dx = wave_eq.derivative()(x0)
        if dz_dx != 0:
            n_slope = np.rad2deg(np.arctan(-1 / dz_dx))

            if np.abs(n_slope) >= 70:
                land_x = x0 + wave_eq(x0) * dz_dx
                land_pts.append(land_x)
    return land_pts

def wave_points(dfa, dist_axis, altitude_axis, center):
    points_l = [] 
    points_w = []
    
    dfa_n = dfa.copy()
    dfa_n[dist_axis] = dfa_n[dist_axis] - center
    group_t = dfa_n.groupby('time')
    
    for time_slice in dfa_n['time'].unique():
        if time_slice != 0:  
            df_ = group_t.get_group(time_slice).sort_values(dist_axis)
            df_ = df_.drop_duplicates(subset=dist_axis).reset_index(drop=True)
            wave_eq = CubicSpline(df_[dist_axis], df_[altitude_axis])
            
            for x0 in df_[dist_axis].to_numpy():
                dz_dx = wave_eq.derivative()(x0)
                
                if dz_dx != 0:
                    n_slope = np.rad2deg(np.arctan(-1 / dz_dx))

                    if np.abs(n_slope) >= 70:
                        land_x = x0 + wave_eq(x0) * dz_dx
                        points_w.append((x0, df_.loc[df_[dist_axis] == x0, altitude_axis].values[0]))
                        
                        if np.abs(land_x) < 65:
                            points_l.append((x0, df_.loc[df_[dist_axis] == x0, altitude_axis].values[0]))

    points_w = list(set(points_w))
    points_l = list(set(points_l))

    points_w = np.array(points_w)
    points_l = np.array(points_l)
    
    return points_w, points_l

def find_nearest_location(locations, current_loc):
    distances = np.abs(np.array(locations) - current_loc)
    nearest_index = np.argmin(distances)
    return locations[nearest_index]

def count_points_in_window(start, end, sorted_array):
    return np.sum((sorted_array >= start) & (sorted_array <= end))

def window_points(xi, window_size, top_no_of_windows):
    top_windows, midpt_windows = [], []
    x_l_sorted = np.sort(xi)
    
    for i in range(len(x_l_sorted)):
        window_start = x_l_sorted[i]
        window_end = window_start + window_size
        count_points = count_points_in_window(window_start, window_end, x_l_sorted)
    
        if len(top_windows) < top_no_of_windows or count_points > top_windows[-1][2]:
            top_windows.append((window_start, window_end, count_points))
            midpt_windows.append((window_start + window_end) / 2)

    top_windows, midpt_windows = zip(*sorted(zip(top_windows, midpt_windows), key=lambda x: x[0][2], reverse=True))
    top_windows = top_windows[:top_no_of_windows]
    midpt_windows = midpt_windows[:top_no_of_windows]

    return top_windows, midpt_windows

def angle_with_x_axis(x, y):
    return np.degrees(np.arctan2(y, x))

def points_in_cone(points, angle_start, angle_width):
    angle_end = angle_start + angle_width
    points_in_cone = []
    
    for point in points:
        x, y = point
        angle = angle_with_x_axis(x, y)
        
        if angle < 0:
            angle += 360
            
        if angle_start <= angle <= angle_end:
            points_in_cone.append(point)
    
    return points_in_cone

def drange(start, stop, step):
    while start < stop:
        yield start
        start += step
        

def main(wf_file_path, dist_axis):
    
    pattern_filename = r'/([^/]+)\_wf.csv$'
    match_filename = re.search(pattern_filename, wf_file_path)
    
    pattern_date = r'F(\d{4})(\d{2})(\d{2})'
    match_date = re.search(pattern_date, wf_file_path)
    
    if match_filename:
        file_name_ = match_filename.group(1)
        # print(file_name_)
    else:
        print("No match found")
        
    if match_date:
        year = match_date.group(1)
        month = match_date.group(2)
        day = match_date.group(3)
        
        year , month , day = int(year), int(month), int(day)
        
        print(f"{year}-{month}-{day}")
        # print(f"Year: {year}, Month: {month}, day: {day}")
    else:
        print("No match found")

    dfa = pd.read_csv(wf_file_path)
    altitude_axis = 'z'

    if dist_axis == 'x':
        dfa_ew = dfa[(dfa['phi'] == 0) & (dfa[altitude_axis] >= 0) & (dfa['theta'] <= np.pi/2) & (dfa['theta'] >= -np.pi/2)].reset_index(drop=True)
        title = 'N-S'
    else:
        dfa_ew = dfa[(dfa['phi'] != 0) & (dfa[altitude_axis] >= 0) & (dfa['theta'] <= np.pi/2) & (dfa['theta'] >= -np.pi/2)].reset_index(drop=True)
        title = 'E-W'
        
    window_size = 130
    top_no_of_windows = 25

    initial_loc = 0
    acoustic_source_loc = []

    group_t = dfa_ew.groupby('time')
    for time_slice in dfa['time'].unique():
        if time_slice != 0:
            df_ = group_t.get_group(time_slice).sort_values(dist_axis).drop_duplicates(subset=dist_axis).reset_index(drop=True)
            land_pts = land_points(df_, dist_axis, altitude_axis)
            
            if land_pts is not None:
                top_windows, midpt_windows = window_points(land_pts, window_size, top_no_of_windows)
                
            if top_windows:
                    nearest_loc = find_nearest_location(midpt_windows, initial_loc)
                    acoustic_source_loc.append(nearest_loc)

    #  may adjust this value to obtain accurate outcomes
    
    eps = 5   
    min_samples = 3
    
    cluster_centers = []

    acoustic_source_loc = np.array(acoustic_source_loc).reshape(-1, 1)
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(acoustic_source_loc) 
    labels = db.labels_
    unique_labels = set(labels)
    unique_labels.discard(-1)

    for k in unique_labels:
        class_member_mask = (labels == k)
        xy = acoustic_source_loc[class_member_mask]
        cluster_center = np.mean(xy)
        cluster_centers.append(cluster_center)

    print(f'\nCluster centers: {cluster_centers}\n')

    angle_width = 3
    angle_start = 68
    angle_end = 113 - angle_width
    angle_step = 1

    points_w, points_l = [], []
    for center in cluster_centers:
        w, l = wave_points(dfa_ew, dist_axis, altitude_axis, center)
        points_l.append(l)
        points_w.append(w)

    angle_f, percent_pts_f = 0, 0

    for angle in drange(angle_start, angle_end, angle_step):
        points_within_cone_w = []
        points_within_cone_l = []
        
        for pw in points_w:
            points_within_cone_w.extend(points_in_cone(pw, angle, angle_width))
        
        for pl in points_l:
            points_within_cone_l.extend(points_in_cone(pl, angle, angle_width))
        
        percent_pts = (len(points_within_cone_l) / len(points_within_cone_w)) * 100 if len(points_within_cone_w) > 0 else 0
        if percent_pts_f <= percent_pts:
            angle_f = angle
            percent_pts_f = percent_pts

    print(f"Percentage inside the cone: {percent_pts_f}, angle between: {angle_f} - {angle_f + angle_width}")

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

    theta1 = np.radians(angle_f)
    theta2 = np.radians(angle_f + angle_width)

    plt.plot([0, np.cos(theta1) * dist], [0, np.sin(theta1) * dist], 'g--', label='Cone-section')
    plt.plot([0, np.cos(theta2) * dist], [0, np.sin(theta2) * dist], 'g--')

    plt.legend()
    save_path = '/home/murali/Documents/rass/automation/gen_data/'
    plt.savefig(save_path + file_name_ + "_" + title + ".png")

    acoustic_sources_file = f"{save_path}acoustic_sources.csv"
    if os.path.exists(acoustic_sources_file):
        df = pd.read_csv(acoustic_sources_file)
    else:
        df = pd.DataFrame(columns=['cluster_Centers', 'axis', 'date', 'persent_pts', 'angle_start'])
    
    df_centers = pd.DataFrame(cluster_centers, columns=['cluster_Centers'])
    df_centers['axis'] = dist_axis
    df_centers['date'] = pd.Timestamp(f"{year}-{month}-{day}")
    df_centers['persent_pts'] = percent_pts_f
    df_centers['angle_start'] = angle_f
    df = pd.concat([df, df_centers])
    df.to_csv(f"{save_path}acoustic_sources.csv", index=False)

if __name__ == "__main__":
    start_time = time.time()
    
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('file_path', type=str, help='Path to the CSV file')
    parser.add_argument('dist_axis', type=str, help='Distance axis (x or y)')
    
    # arser = argparse.ArgumentParser(description='Process some integers.')
    # parser.add_argument('--filepath', type=str, help='file path of the wavefronts data')
    # parser.add_argument('--dist_axis', type=str, help='Horizontal axis')
    # parser.add_argument('--top_windows', default=25, type=int, help='considering the top number of windows')
    # parser.add_argument('--window_size', default=130, type=float, help='length of the antenna array')
    # parser.add_argument('--step_coverage_angle', default=3, type=float, help='max coverage angle of the antenna array at a time')
    # parser.add_argument('--max_phase_coverage', default=40, type=float, help='max phase coverage of the antenna array')
    # parser.add_argument('--angle_step_size', default=1, type=float, help='step size of the cone itterator')
    # parser.add_argument('--eps', default=5, type=float, help='eps of the DBSCAN clustring')
    # parser.add_argument('--min_samples', default=3, type=int, help='min samples of the DBCSAN clustring')
    
    args = parser.parse_args()
    main(args.file_path, args.dist_axis)
    
    print(f"Time taken: {time.time() - start_time}")
