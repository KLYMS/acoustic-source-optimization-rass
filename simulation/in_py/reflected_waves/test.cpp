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

// Tuning Parameters
double rho = 10;  // Smoothing parameter
double offset_dist = 0; // offset dist for filtering
const double SpeedOfSound = 343;

const bool data_filtering = true;
const bool xy_filt = true;

void log_message(const string &message);
void readCSV(vector<vector<double>>& columns, const string& filename);
void readTimeSlices(vector<int>& time_slices, const string& filename);
void filterData(vector<double>& x, vector<double>& y, vector<double>& z, int time_slice, double offset_dist);
void fitSplineSurface(alglib::spline2dinterpolant& s, const vector<double>& x, const vector<double>& y, const vector<double>& z);
void findNormalAtPoint(const alglib::spline2dinterpolant& s, double x_eval, double y_eval, double& x_land, double& y_land);

int main() {
    auto start = high_resolution_clock::now();

    string time_data = "/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/time_step_cpp_alglib_CS_i_0.txt";
    string file_path = "/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/wTemp_wWind_cpp_alglib_CS_i_0_";

    vector<int> time_slices;
    readTimeSlices(time_slices, time_data);

    ofstream output_file("/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/backscatter_coordinates_i_0_G200-r10.csv");
    output_file << "x_val,y_val,z_val,x_land,y_land,time_step,phi,theta\n";

    for (size_t idx = 0; idx < time_slices.size(); ++idx) {
        if (time_slices[idx] != 0) {
            log_message("Processing time step: " + to_string(time_slices[idx]) + " (" + to_string(idx + 1) + "/" + to_string(time_slices.size()) + ")");

            vector<vector<double>> data;
            readCSV(data, file_path + to_string(time_slices[idx]) + ".csv");

            vector<double> x = data[0];
            vector<double> y = data[1];
            vector<double> z = data[2];

            if (data_filtering) filterData(x, y, z, time_slices[idx], offset_dist);

            alglib::spline2dinterpolant s;

            try {
                fitSplineSurface(s, x, y, z);
            } catch (const alglib::ap_error& e) {
                log_message("ALGLIB error during spline fitting: " + string(e.msg));
                continue;
            }

            for (size_t idx_x = 0; idx_x < x.size(); ++idx_x) {
                double x_land, y_land;
                double x_eval = x[idx_x];
                double y_eval = y[idx_x];
                double z_eval = z[idx_x];

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
                                << data[3][idx_x] << "\n"; // theta

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
    ofstream log_file("/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/backscattering_progress_i_0_G200-r10.log", ios_base::app);
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
            continue; // Skip the header row
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

// Function to filter data based on z values
void filterData(vector<double>& x, vector<double>& y, vector<double>& z, int time_slice, double offset_dist) {
    vector<double> x_filtered, y_filtered, z_filtered;
    double max_dist = SpeedOfSound * time_slice + offset_dist;

    if (!xy_filt) {
        for (size_t i = 0; i < z.size(); ++i) {
            if (z[i] > 0) {
                x_filtered.push_back(x[i]);
                y_filtered.push_back(y[i]);
                z_filtered.push_back(z[i]);
            }
        }
    } else {
        for (size_t i = 0; i < z.size(); ++i) {
            if (z[i] > 0 && abs(x[i]) < max_dist && abs(y[i]) < max_dist) {
                x_filtered.push_back(x[i]);
                y_filtered.push_back(y[i]);
                z_filtered.push_back(z[i]);
            }
        }
    }

    x = x_filtered;
    y = y_filtered;
    z = z_filtered;
}

void fitSplineSurface(alglib::spline2dinterpolant& s, const vector<double>& x, const vector<double>& y, const vector<double>& z) {
    int n = x.size();

    // Check if we have enough data points
    if (n < 2) { // Need at least 2 points for spline fitting
        log_message("Not enough data points for spline fitting: " + to_string(n));
        throw alglib::ap_error("Not enough data points");
    }

    // Create real_1d_array objects for x, y, and z values
    alglib::real_1d_array x_array, y_array, f_array;
    x_array.setlength(n);
    y_array.setlength(n);
    f_array.setlength(n);

    // Fill in x, y, and f arrays
    for (int i = 0; i < n; ++i) {
        x_array[i] = x[i];
        y_array[i] = y[i];
        f_array[i] = z[i];
    }

    alglib::spline2dbuildbicubicv(x_array, n, y_array, n, f_array, 1, s);

    log_message("Spline fitting completed.");
}

void findNormalAtPoint(const alglib::spline2dinterpolant& s, double x_eval, double y_eval, double& x_land, double& y_land) {
    double dzdx, dzdy, z_eval;
    alglib::spline2ddiff(s, x_eval, y_eval, z_eval, dzdx, dzdy);

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
