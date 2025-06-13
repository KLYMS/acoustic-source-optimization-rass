#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include "/home/murali/Documents/rass/data/cpp_lib/alglib-cpp/src/interpolation.h"

using namespace std;
using namespace alglib;
using namespace std::chrono;

void log_message(const string &message);
int wrap_index(int idx, int total_length);
void compute_normal_vector(rbfmodel &model, double x, double y, double &nx, double &ny, double &nz);
int read_csv(vector<vector<double>> &columns, const string &filename);
void process_time_step(int t_step, int total_steps, int step_idx, const string &file_path, ofstream &output_file, double epsilon);

int main() {
    double epsilon = 1e-8;  // Small step for finite difference
    vector<int> time_steps;

    // Read the time steps from the file
    ifstream time_step_file("/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/time_step_cpp_alglib_CS_0.txt");
    string line;
    while (getline(time_step_file, line)) {
        if (!line.empty()) {
            try {
                time_steps.push_back(stoi(line));
            } catch (const invalid_argument &e) {
                cerr << "Invalid argument: " << line << endl;
                log_message("Invalid argument: " + line);
            } catch (const out_of_range &e) {
                cerr << "Out of range error: " << line << endl;
                log_message("Out of range error: " + line);
            }
        }
    }
    time_step_file.close();

    string file_path = "/home/murali/Documents/rass/data/sim_data/dived_data/cpp_data/wTemp_wWind_cpp_alglib_CS_0_";

    auto start_time = chrono::high_resolution_clock::now();

    ofstream output_file("/home/murali/Documents/rass/simulation/in_py/reflected_waves/backscatter_coordinates.csv");
    output_file << "x_val,y_val,z_val,x_land,y_land,time_step,phi,theta\n";

    for (size_t idx = 0; idx < time_steps.size(); ++idx) {
        if (time_steps[idx] != 0) {
            process_time_step(time_steps[idx], time_steps.size(), idx, file_path, output_file, epsilon);
        }
    }

    output_file.close();

    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_time = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();

    log_message("Finished writing results to CSV");
    log_message("Total time taken: " + to_string(elapsed_time) + " seconds");

    return 0;
}

// Function to log messages
void log_message(const string &message) {
    ofstream log_file("backscattering_progress.log", ios_base::app);
    log_file << message << endl;
}

// Function to handle wrap-around index
int wrap_index(int idx, int total_length) {
    return (idx + total_length) % total_length;
}

// Function to compute normal vector
void compute_normal_vector(rbfmodel &model, double x, double y, double &nx, double &ny, double &nz) {
    double epsilon = 1e-8;  // Small step for finite difference
    double dz_dx = (rbfcalc2(model, x + epsilon, y) - rbfcalc2(model, x - epsilon, y)) / (2 * epsilon);
    double dz_dy = (rbfcalc2(model, x, y + epsilon) - rbfcalc2(model, x, y - epsilon)) / (2 * epsilon);

    nx = -dz_dx;
    ny = -dz_dy;
    nz = 1.0;
    double norm = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= norm;
    ny /= norm;
    nz /= norm;
}

// Function to read CSV file
int read_csv(vector<vector<double>> &columns, const string &filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return 1;
    }

    string line;
    getline(file, line);
    stringstream headerStream(line);
    string headerCell;

    while (getline(headerStream, headerCell, ',')) {
        columns.emplace_back();  // Add an empty vector for each column
    }

    while (getline(file, line)) {
        stringstream lineStream(line);
        string cell;
        size_t colIndex = 0;

        while (getline(lineStream, cell, ',')) {
            try {
                cell.erase(remove_if(cell.begin(), cell.end(), ::isspace), cell.end());
                double value = stod(cell);
                columns[colIndex].push_back(value);
            } catch (const invalid_argument &e) {
                cerr << "Invalid input for column " << colIndex + 1 << ": " << cell << endl;
            } catch (const out_of_range &e) {
                cerr << "Input out of range for column " << colIndex + 1 << ": " << cell << endl;
            }
            colIndex++;
        }
    }

    file.close();
    return 0;
}

// Function to process time step
void process_time_step(int t_step, int total_steps, int step_idx, const string &file_path, ofstream &output_file, double epsilon) {
    log_message("Processing time step: " + to_string(t_step) + " (" + to_string(step_idx + 1) + "/" + to_string(total_steps) + ")");
    
    vector<vector<double>> data;
    read_csv(data, file_path + to_string(t_step) + ".csv");

    int total_operations = 0;

    for (size_t idx = 0; idx < data.size(); ++idx) {
        double phi = data[idx][5];
        double theta = data[idx][3];
        vector<vector<double>> points = {data[idx]};

        for (int d_phi = -1; d_phi <= 1; ++d_phi) {
            for (int d_theta = -1; d_theta <= 1; ++d_theta) {
                if (d_phi == 0 && d_theta == 0) continue;

                int idx_phi = wrap_index(idx + d_phi, data.size());
                int idx_theta = wrap_index(idx + d_theta, data.size());
                points.push_back(data[idx_phi]);
                points.push_back(data[idx_theta]);
            }
        }

        if (points.size() < 9) {
            log_message("Not enough points found for phi=" + to_string(phi) + " and theta=" + to_string(theta) + " ------ num of pts: " + to_string(points.size()));
            continue;
        }
        
        real_2d_array pts;
        pts.setlength(points.size(), 3);
        for (size_t i = 0; i < points.size(); ++i) {
            pts[i][0] = points[i][0];
            pts[i][1] = points[i][1];
            pts[i][2] = points[i][2];
        }

        rbfmodel model;
        rbfreport rep;
        rbfcreate(2, 1, model);
        rbfsetpoints(model, pts);
        rbfsetalgothinplatespline(model);
        rbfbuildmodel(model, rep);

        double center_x = points[0][1], center_y = points[0][2], center_z = points[0][3];
        double nx, ny, nz;
        compute_normal_vector(model, center_x, center_y, nx, ny, nz);

        double t = -rbfcalc2(model, center_x, center_y) / (nx * (rbfcalc2(model, center_x + epsilon, center_y) - rbfcalc2(model, center_x, center_y)) / epsilon +
                                                           ny * (rbfcalc2(model, center_x, center_y + epsilon) - rbfcalc2(model, center_x, center_y)) / epsilon);

        double x_zero = center_x + t * nx;
        double y_zero = center_y + t * ny;

        if (abs(x_zero) < 130 && abs(y_zero) < 130) {
            output_file << fixed << setprecision(6)
                        << center_x << ","
                        << center_y << ","
                        << center_z << ","
                        << x_zero << ","
                        << y_zero << ","
                        << t_step << ","
                        << phi << ","
                        << theta << "\n";

            log_message("time_step: " + to_string(t_step) + ", Coordinates (x, y) where backscatter at z = 0: (" + to_string(x_zero) + ", " + to_string(y_zero) + ")");
        }

        total_operations++;
    }

    log_message("Finished processing time step: " + to_string(t_step) + ", Total operations: " + to_string(total_operations));
}
