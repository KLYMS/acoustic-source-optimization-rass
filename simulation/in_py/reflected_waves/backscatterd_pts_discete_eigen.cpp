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
double rho = 3;  // Smoothing parameter
alglib::ae_int_t grid_size = 25; // Grid size parameter
double offset_dist = 1500; // Offset distance for filtering
const double SpeedOfSound = 343;

const bool data_filtering = true;
const bool xy_filt = true;

void log_message(const string &message);
void readCSV(vector<vector<double>>& columns, const string& filename);
void readTimeSlices(vector<int>& time_slices, const string& filename);
int wrap_index(int idx, int total_length);
void collectPoints(const vector<vector<double>>& data, int idx, vector<vector<double>>& points);
void filterData(vector<vector<double>>& data, int time_slice, double offset_dist);
void fitSplineSurface(alglib::spline2dinterpolant& s, const vector<double>& x, const vector<double>& y, const vector<double>& z, double rho, alglib::ae_int_t grid_size);
void findNormalAtPoint(const alglib::spline2dinterpolant& s, double x_eval, double y_eval, double& x_land, double& y_land);

int main() {
    auto start = high_resolution_clock::now();

    string time_data = "/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/time_step_cpp_alglib_CS_i_0.txt";
    string file_path = "/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/wTemp_wWind_cpp_alglib_CS_i_0_";

    vector<int> time_slices;
    readTimeSlices(time_slices, time_data);

    ofstream output_file("/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/backscatter_coordinates_i_0_d_G25-r3.csv");
    output_file << "x_val,y_val,z_val,x_land,y_land,time_step,phi,theta\n";

    for (size_t idx = 0; idx < time_slices.size(); ++idx) {
        if (time_slices[idx] != 0) {
            log_message("Processing time step: " + to_string(time_slices[idx]) + " (" + to_string(idx + 1) + "/" + to_string(time_slices.size()) + ")");

            vector<vector<double>> data;
            readCSV(data, file_path + to_string(time_slices[idx]) + ".csv");

            if (data_filtering) filterData(data, time_slices[idx], offset_dist);

            if (data.empty() || data[0].empty()){
                log_message(to_string(time_slices[idx]) + " data is empty");
                continue;
            };

            size_t data_size = data[0].size();
            for (size_t i = 0; i < data_size; ++i) {
                double phi = data[5][i];
                double theta = data[3][i];

                vector<vector<double>> points;

                points.push_back({ data[0][i], data[1][i], data[2][i], data[3][i], data[5][i]});

                collectPoints(data,i,points);

                

                if (points[0].size() < 9) {
                    log_message("Not enough points found for phi=" + to_string(phi) + " and theta=" + to_string(theta) + " ------ num of pts: " + to_string(points[0].size()));
                    continue;
                }

                alglib::spline2dinterpolant s;
                try {
                    fitSplineSurface(s, points[0], points[1], points[2], rho, grid_size);
                } catch (const alglib::ap_error& e) {
                    log_message("ALGLIB error during spline fitting: " + string(e.msg));
                    continue;
                }   

                double x_land, y_land;
                double x_eval = points[0][0];
                double y_eval = points[1][0];
                double z_eval = points[2][0];

                findNormalAtPoint(s, x_eval, y_eval, x_land, y_land);

                if (abs(x_land) <= 130 && abs(y_land) <= 130) {
                    output_file << fixed << setprecision(6)
                                << x_eval << ","
                                << y_eval << ","
                                << z_eval << ","
                                << x_land << ","
                                << y_land << ","
                                << time_slices[idx] << ","
                                << phi << ","
                                << theta << "\n";
                    log_message("Got point");
                }
            }
        }
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end - start);
    log_message("Time Taken: " + to_string(duration.count()) + " seconds");

    return 0;
}

void log_message(const string &message) {
    ofstream log_file("/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/alglib_CS_i_0/backscattering_progress_i_0_d_G25-r3.log", ios_base::app);
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

// Function to handle wrap-around index
int wrap_index(int idx, int total_length) {
    return (idx % total_length + total_length) % total_length;
}

void collectPoints(const vector<vector<double>>& data, int idx, vector<vector<double>>& points) {
    int total_length = data[0].size(); // Number of data points

    for (int d_phi = -1; d_phi <= 1; ++d_phi) {
        for (int d_theta = -1; d_theta <= 1; ++d_theta) {
            if (d_phi == 0 && d_theta == 0) {
                continue;
            }

            int idx_phi = wrap_index(idx + d_phi, total_length);
            int idx_theta = wrap_index(idx + d_theta, total_length);

            // data format: data[column_index][row_index]
            points.push_back({ data[0][idx_phi], data[1][idx_theta], data[2][idx_theta], data[3][idx_theta], data[5][idx_phi]});
        }
    }
}


void fitSplineSurface(alglib::spline2dinterpolant& s, const vector<double>& x, const vector<double>& y, const vector<double>& z, double rho, alglib::ae_int_t grid_size) {
    int n = x.size();

    if (n < 3) {
        log_message("Not enough data points for spline fitting: " + to_string(n));
        throw alglib::ap_error("Not enough data points");
    }

    alglib::real_2d_array xy;
    xy.setlength(n, 3); 
    for (int i = 0; i < n; ++i) {
        xy(i, 0) = x[i];
        xy(i, 1) = y[i];
        xy(i, 2) = z[i];
    }

    alglib::real_1d_array z_array;
    z_array.setcontent(n, &z[0]);

    alglib::spline2dbuilder builder;
    alglib::spline2dfitreport rep;
    alglib::spline2dbuildercreate(1, builder); 
    alglib::spline2dbuildersetpoints(builder, xy, n);
    alglib::spline2dbuildersetgrid(builder, grid_size, grid_size); 
    alglib::spline2dbuildersetalgoblocklls(builder, rho); 
    alglib::spline2dfit(builder, s, rep);

    log_message("Spline fitting completed. RMS error: " + to_string(rep.rmserror));
}

void findNormalAtPoint(const alglib::spline2dinterpolant& s, double x_eval, double y_eval, double& x_land, double& y_land) {
    double dzdx, dzdy, z_eval;
    alglib::spline2ddiff(s, x_eval, y_eval, z_eval, dzdx, dzdy);

    Vector3d gradient(dzdx, dzdy, -1.0);
    Vector3d normal = gradient.normalized();

    double t = -z_eval / normal(2);

    x_land = x_eval + t * normal(0);
    y_land = y_eval + t * normal(1);
}
