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

using namespace std;
using namespace std::chrono;
using namespace Eigen;
using namespace alglib;

// Constants and global variables
double rho = 10;  // Smoothing parameter
alglib::ae_int_t grid_size = 200; // Grid size parameter
double offset_dist = 0; // offset dist for filtering
const double SpeedOfSound = 343;
const bool data_filtering = true;
const bool xy_filt = true;

void log_message(const string &message);
void readCSV(vector<vector<double>>& columns, const string& filename);
void readTimeSlices(vector<int>& time_slices, const string& filename);
void filterData(vector<vector<double>>& data, int time_slice, double offset_dist);
void fitBicubicSurface(bicubicspline2dinterpolant& s, const vector<double>& x, const vector<double>& y, const vector<double>& z, alglib::ae_int_t grid_size);
void findNormalAtPoint(const bicubicspline2dinterpolant& s, double x_eval, double y_eval, double& x_land, double& y_land);

int main() {
    auto start = high_resolution_clock::now();

    string time_data = "/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/time_step_cpp_alglib_CS_i_0.txt";
    string file_path = "/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/wTemp_wWind_cpp_alglib_CS_i_0_";

    vector<int> time_slices;
    readTimeSlices(time_slices, time_data);

    ofstream output_file("backscatter_coordinates_i_0_G200-r10.csv");
    output_file << "x_val,y_val,z_val,x_land,y_land,time_step,phi,theta\n";

    for (size_t idx = 0; idx < time_slices.size(); ++idx) {
        if (time_slices[idx] != 0) {
            log_message("Processing time step: " + to_string(time_slices[idx]) + " (" + to_string(idx + 1) + "/" + to_string(time_slices.size()) + ")");

            vector<vector<double>> data;
            readCSV(data, file_path + to_string(time_slices[idx]) + ".csv");

            if (data_filtering) filterData(data, time_slices[idx], offset_dist);

            bicubicspline2dinterpolant s;

            try {
                fitBicubicSurface(s, data[0], data[1], data[2], grid_size);
            } catch (const alglib::ap_error& e) {
                log_message("ALGLIB error during bicubic fitting: " + string(e.msg));
                continue;
            }

            for (size_t idx_x = 0; idx_x < data[0].size(); ++idx_x) {
                double x_land, y_land;
                double x_eval = data[0][idx_x];
                double y_eval = data[1][idx_x];
                double z_eval = data[2][idx_x];

                findNormalAtPoint(s, x_eval, y_eval, x_land, y_land);

                if (abs(x_land) <= 130 && abs(y_land) <= 130) {
                    output_file << fixed << setprecision(6)
                        << x_eval << ","
                        << y_eval << ","
                        << z_eval << ","
                        << x_land << ","
                        << y_land << ","
                        << time_slices[idx] << ","
                        << data[5][idx_x] << ","  // phi
                        << data[3][idx_x] << "\n";  // theta

                    // log_message("time_step: " + to_string(time_slices[idx]) + ", Coordinates (x, y) where backscatter at z = 0: (" + to_string(x_land) + ", " + to_string(y_land) + ")");
                }
            }
        }
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end - start);
    log_message("Time Taken: " + to_string(duration.count()) + " seconds");

    return 0;
}

// Function to log messages
void log_message(const string &message) {
    ofstream log_file("ackscattering_progress_test.log", ios_base::app);
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

void fitBicubicSurface(const vector<double>& x, const vector<double>& y, const vector<double>& z, const double rho, const ae_int_t grid_size, spline2dinterpolant& s) {
    int n = x.size();

    // Check if we have enough data points
    if (n < 4) {  // Need at least 4 points for bicubic spline fitting
        // Handle error, not enough data points
        return;
    }

    // ALGLIB objects
    real_1d_array x_array, y_array;
    real_2d_array z_array;
    x_array.setcontent(n, &x[0]);
    y_array.setcontent(n, &y[0]);
    z_array.setcontent(n, 3, &z[0]);

    // Fit bicubic spline to data
    spline2dinterpolant s_temp;
    Spline2DBuildBicubicV(x_array, y_array, z_array, n, 3, s_temp);

    // Copy the spline interpolant to output parameter
    s = s_temp;
}

void findNormalAtPoint(const bicubicspline2dinterpolant& s, double x_eval, double y_eval, double& x_land, double& y_land) {
    double dzdx, dzdy, z_eval;
    alglib::bicubicspline2ddiff(s, x_eval, y_eval, z_eval, dzdx, dzdy);

    // Gradient vector (∂z/∂x, ∂z/∂y, -1)
    Vector3d gradient(dzdx, dzdy, -1.0);

    // Normal vector is the normalized gradient
    Vector3d normal = gradient.normalized();

    // Parametric form of the normal line: (x, y, z) = (x_eval, y_eval, z_eval) + t * normal
    // Solve for t when z = 0:
    double t = -z_eval / normal(2);

    // Find (x, y) coordinates at z = 0
    x_land = x_eval + t * normal(0);
    y_land = y_eval + t * normal(1);
}

