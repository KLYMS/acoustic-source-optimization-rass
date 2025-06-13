#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <Eigen/Dense>
#include <algorithm>
#include "/home/murali/Documents/rass/data/cpp_lib/alglib-cpp/src/interpolation.h"
#include "/home/murali/Documents/rass/data/cpp_lib/dbscan/dbscan.hpp"

using namespace std;
using namespace std::chrono;
using namespace Eigen;
using namespace alglib;

struct Point{
    double x ,y ;
};

// Tuning Parameters
double rho = 2;  // Smoothing parameter
alglib::ae_int_t grid_size = 250; // Grid size parameter

double offset_dist = 0; // offset dist for filtering
const double SpeedOfSound = 343;

int top_no_of_windows = 25;
double window_size = 130; 
Point initial_loc = {0,0}

// Parameters for DBSCAN
double eps = 5.0;  // The radius to search for neighbors
int min_pts = 3;   // The minimum number of points to form a cluster

const int file_no = 0; 

const bool data_filtering = true;
const bool xy_filt = true;

void log_message(const string &message);
void readCSV(vector<vector<double>>& columns, const string& filename);
void readTimeSlices(vector<int>& time_slices, const string& filename);
void filterData(vector<vector<double>>& data, int time_slice, double offset_dist);
void fitSplineSurface(alglib::spline2dinterpolant& s, const vector<double>& x, const vector<double>& y, const vector<double>& z, double rho, alglib::ae_int_t grid_size);
void findNormalAtPoint(const alglib::spline2dinterpolant& s, double x_eval, double y_eval, double& x_land, double& y_land);
int countPointsInWindow(const vector<Point>& points, double centerX, double centerY, double windowSize);
vector<Point> double_window_points(const vector<Point>& points, double window_size, int top_no_of_windows);
double calculateDistance(const Point& p1, const Point& p2);
Point findNearestWindow(const vector<Point>& top_windows, const Point& current_loc);


int main(){
    auto start = high_resolution_clock::now();

    string time_data = "/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_" + to_string(file_no) + "/time_step_cpp_alglib_CS_" + to_string(file_no) + ".txt";
    string file_path = "/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_" + to_string(file_no) + "/wTemp_wWind_cpp_alglib_CS_" + to_string(file_no) + "_";

    vector<int> time_slices;
    readTimeSlices(time_slices, time_data);

    ofstream output_file("back_coord.csv");
    output_file << "x_val,y_val,z_val,x_land,y_land,time_step,phi,theta\n";

    vector<Point> acoustic_source_loc;

    for (size_t idx = 0; idx < time_slices.size(); ++idx) {
        if (time_slices[idx] != 0) {
            log_message("Processing time step: " + to_string(time_slices[idx]) + " (" + to_string(idx + 1) + "/" + to_string(time_slices.size()) + ")");

            vector<vector<double>> data;
            vector<vector<Point>> land_points;

            readCSV(data, file_path + to_string(time_slices[idx]) + ".csv");

            if (data_filtering) filterData(data, time_slices[idx], offset_dist);

            vector<double> x = data[0];
            vector<double> y = data[1];
            vector<double> z = data[2];

            log_message("Done Data Filtering\nNow fitting the Surface having grid size: " + to_string(grid_size) + ", rho: " + to_string(rho));

            alglib::spline2dinterpolant s;

            try {
                fitSplineSurface(s, x, y, z, rho, grid_size);
            } catch (const alglib::ap_error& e) {
                log_message("ALGLIB error during spline fitting: " + string(e.msg));
                continue;
            }
            log_message("Done Surface Spline Fitting\nNow Evalating the Surface ");

            for (size_t idx_x = 0; idx_x < x.size(); ++idx_x) {
                double x_land, y_land;
                double x_eval = x[idx_x];
                double y_eval = y[idx_x];
                double z_eval = z[idx_x];

                findNormalAtPoint(s, x_eval, y_eval, x_land, y_land);
                Point p = {x_land, y_land};
                land_points.push_back(p);
            }
            
            vector<Point> top_windows;
            top_windows = double_window_points(land_points, window_size, top_no_of_windows);
            Point nearest_loc;
            nearest_loc = findNearestWindow(top_windows, initial_loc);
            acoustic_source_loc.push_back(nearest_loc);
        }
    }

    // Convert your Point data to point2 data
    vector<point2> data2D;
    for (const auto& p : acoustic_source_loc) {
        data2D.push_back({p.x, p.y});
    }

    // Perform DBSCAN clustering
    auto clusters2D = dbscan( span<const point2>(data2D.data(), data2D.size()), eps, min_pts);

    // Output the clusters and calculate centers
    std::cout << "2D Clusters:\n";
    for (const auto& cluster : clusters2D) {
        float sum_x = 0.0f;
        float sum_y = 0.0f;
        for (const auto& idx : cluster) {
            sum_x += data2D[idx].x;
            sum_y += data2D[idx].y;
            std::cout << "(" << data2D[idx].x << ", " << data2D[idx].y << ") ";
        }
        std::cout << "\n";
        float center_x = sum_x / cluster.size();
        float center_y = sum_y / cluster.size();
        std::cout << "Cluster Center: (" << center_x << ", " << center_y << ")\n";
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end - start);    
    cout << "Time Taken: " << to_string(duration.count()) << " seconds" << endl;
}

// Function to log messages
void log_message(const string &message) {
    // ofstream log_file("/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_" + to_string(file_no) + "/backscattering_progress_" + to_string(file_no) + "_G250-r2.log", ios_base::app);
    ofstream log_file("back_coord.log");
    log_file << message << endl;
}

void readCSV(vector<vector<double>>& columns, const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        log_message("Error opening file: " + filename);
        return;
    }

    string line;
    bool header = true;
    while (getline(file, line)) {
        if (header) {
            header = false;
            continue;  // Skip the header row
        }

        stringstream lineStream(line);
        string cell;
        vector<double> parsedRow;
        while (getline(lineStream, cell, ',')) {
            try {
                parsedRow.push_back(stod(cell));
            } catch (const invalid_argument& e) {
                log_message("Invalid input: " + cell + " in file: " + filename);
                parsedRow.push_back(0.0); // Default to 0.0 for invalid input
            } catch (const out_of_range& e) {
                log_message("Input out of range: " + cell + " in file: " + filename);
                parsedRow.push_back(0.0); // Default to 0.0 for out of range input
            }
        }
        if (columns.empty()) {
            columns.resize(parsedRow.size());
        }
        for (size_t i = 0; i < parsedRow.size(); ++i) {
            columns[i].push_back(parsedRow[i]);
        }
    }
    file.close();
}

void readTimeSlices(vector<int>& time_slices, const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        log_message("Error opening file: " + filename);
        return;
    }

    string line;
    while (getline(file, line)) {
        if (!line.empty()) {
            try {
                time_slices.push_back(stod(line));
            } catch (const invalid_argument& e) {
                log_message("Invalid argument: " + line);
            } catch (const out_of_range& e) {
                log_message("Out of range error: " + line);
            }
        }
    }
    file.close();
}

void filterData(vector<vector<double>>& data, int time_slice, double offset_dist) {
    vector<vector<double>> data_filtered;
    double max_dist = SpeedOfSound * time_slice + offset_dist;

    for (size_t i = 0; i < data[2].size(); ++i) {
        if (data[2][i] > 0 && (!xy_filt || (data[0][i] < max_dist && data[1][i] < max_dist))) {
            vector<double> row;
            for (size_t j = 0; j < data.size(); ++j) {
                row.push_back(data[j][i]);
            }
            data_filtered.push_back(row);
        }
    }

    // Transpose the filtered data back to the original format
    vector<vector<double>> transposed_data(data.size());
    for (const auto& row : data_filtered) {
        for (size_t j = 0; j < row.size(); ++j) {
            transposed_data[j].push_back(row[j]);
        }
    }

    data = move(transposed_data);
}


void fitSplineSurface(alglib::spline2dinterpolant& s, const vector<double>& x, const vector<double>& y, const vector<double>& z, double rho, alglib::ae_int_t grid_size) {
    int n = x.size();

    // Check if we have enough data points
    if (n < 3) {  // Need at least 3 points for spline fitting
        log_message("Not enough data points for spline fitting: " + to_string(n));
        throw alglib::ap_error("Not enough data points");
    }

    // ALGLIB objects
    alglib::real_2d_array xy;
    xy.setlength(n, 3); 
    for (int i = 0; i < n; ++i) {
        xy(i, 0) = x[i];
        xy(i, 1) = y[i];
        xy(i, 2) = z[i];
    }

    alglib::real_1d_array z_array;
    z_array.setcontent(n, &z[0]);

    // Fit spline to data
    alglib::spline2dbuilder builder;
    alglib::spline2dfitreport rep;
    alglib::spline2dbuildercreate(1, builder); 
    alglib::spline2dbuildersetpoints(builder, xy, n);
    alglib::spline2dbuildersetgrid(builder, grid_size, grid_size); 
    alglib::spline2dbuildersetalgoblocklls(builder, rho); // Use BlockLLS algorithm with smoothing parameter rho
    alglib::spline2dfit(builder, s, rep);

   log_message("Spline fitting completed. RMS error: " + to_string(rep.rmserror));
}

void findNormalAtPoint(const alglib::spline2dinterpolant& s, double x_eval, double y_eval, double& x_land, double& y_land) {
    double dzdx, dzdy, z_eval;
    alglib::spline2ddiff(s, x_eval, y_eval, z_eval, dzdx, dzdy);

    // Gradient vector (∂z/∂x, ∂z/∂y, -1)
    Vector3d gradient(dzdx, dzdy, -1.0);

    // Normal vector is the normalized gradient
    Vector3d normal = gradient.normalized();

    // log_message("Normal vector at (" + to_string(x_eval) + ", " + to_string(y_eval) + ", " + to_string(z_eval) + ") is (" + to_string(normal(0)) + ", " + to_string(normal(1)) + ", " + to_string(normal(2)) + ")");

    // Parametric form of the normal line: (x, y, z) = (x_eval, y_eval, z_eval) + t * normal
    // Solve for t when z = 0:
    double t = -z_eval / normal(2);

    // Find (x, y) coordinates at z = 0
    x_land = x_eval + t * normal(0);
    y_land = y_eval + t * normal(1);
}

int countPointsInWindow(const vector<Point>& points, double centerX, double centerY, double windowSize){
    int half_size =  windowSize / 2 ;

    int count = 0;
    for (const auto& p : points){
        if (abs(p.x - centerX) <= half_size && abs(p.y - centerY) <= half_size){
            count++;
        }
    }
    return count;
}

vector<Point> double_window_points(const vector<Point>& points, double window_size, int top_no_of_windows) {
    vector<pair<Point, int>> point_counts;

    // Count points within the window for each point
    for (const auto& p : points) {
        int pts_inside_window = countPointsInWindow(points, p.x, p.y, window_size);
        point_counts.emplace_back(p, pts_inside_window);
    }

    // Sort the points based on the count in descending order
    sort(point_counts.begin(), point_counts.end(), [](const pair<Point, int>& a, const pair<Point, int>& b) {
        return b.second < a.second; // Sort in descending order of count
    });

    // Select the top points based on the count
    vector<Point> top_pts;
    for (int i = 0; i < min(top_no_of_windows, static_cast<int>(point_counts.size())); ++i) {
        top_pts.push_back(point_counts[i].first);
    }

    return top_pts;
}

// Function to calculate the Euclidean distance between two points
double calculateDistance(const Point& p1, const Point& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Function to find the nearest window (or location) from a list of potential windows
Point findNearestWindow(const vector<Point>& top_windows, const Point& current_loc) {
    double min_distance = numeric_limits<double>::max();
    Point nearest_window;
    
    for (const auto& window : top_windows) {
        double distance = calculateDistance(window, current_loc);
        if (distance < min_distance) {
            min_distance = distance;
            nearest_window = window;
        }
    }

    return nearest_window;
}

