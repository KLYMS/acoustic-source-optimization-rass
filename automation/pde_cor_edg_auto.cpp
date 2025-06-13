#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <chrono>
#include <regex>
#include <sstream>
#include "/home/murali/Documents/rass/data/cpp_lib/alglib-cpp/src/interpolation.h"

using namespace std;
using namespace std::chrono;
using namespace alglib;
    
// Constants
const double Gamma = 1.4;
const double SpeedOfSound = 343;
const double R_air = 287;
const double Troposphere_Lapse_Rate = -6.5;
const double Tropopause_Height = 18000;
const double Stratosphere_Lapse_Rate = 2.5;
const double Stratopause_Height = 48000;

// Flags
const bool Ideal_Conditions = false;
const bool With_temperature = true;
const bool With_wind = true;

const bool With_spline = true;
const bool akima_spline = true; // using akima spline interpolation
const bool zero_insertion = false;

// Function declarations
int read_csv(vector<vector<double>>& columns, string& filename);
double Temperature_(double z);
double Speed_Of_Sound(double temperature);
double avg_data(const vector<double>& height_data, const vector<double>& data, double target_height);
double forward_calc(int index, const vector<double>& height_data, const vector<double>& temp_data);
double backward_calc(int index, const vector<double>& height_data, const vector<double>& temp_data);
double slope(const vector<double>& height_data, const vector<double>& data, double target_height);
vector<double> ray_tracing_eqns(vector<double> y_, double phi, const vector<double>& height_data, const vector<double>& temperature_data, const vector<double>& wind_x_data, const vector<double>& wind_y_data);
vector<double> Rk4Solver(vector<double> y, double phi, double h, const vector<double>& height_data, const vector<double>& temperature_data, const vector<double>& wind_x_data, const vector<double>& wind_y_data);
tuple<vector<vector<double>>, vector<double>, vector<double>> integrate_(vector<double> y, double km, double h, double phi, const vector<double>& height_data, const vector<double>& temperature_data, const vector<double>& wind_x_data, const vector<double>& wind_y_data);
vector<double> drange(double x, double y, double jump);
void write_to_csv(const string& filename, const vector<vector<double>>& data);
void build_spline(const vector<double>& x, const vector<double>& y, spline1dinterpolant& c);

// ALGLIB spline objects
spline1dinterpolant temperature_spline, wind_x_spline, wind_y_spline;

int main(int argc, char* argv[]) {
    chrono::high_resolution_clock::time_point start = high_resolution_clock::now();
    
    // inputs default values
    string input_filepath;
    double phi_step_deg = 90;

    if (argc > 1) {
        input_filepath = argv[1];
    } else {
        cout << "No file path provided" << endl;
    }

    if (argc > 2) {
        phi_step_deg = stod(argv[2]);
    }

    string output_filepath = "/home/murali/Documents/rass/automation/gen_data/";
    string base_filename;

    regex pattern(R"(/([^/]+)_cleaned\.csv$)");
    smatch match_filename;

    if (regex_search(input_filepath, match_filename ,pattern)){
        base_filename = match_filename[1];
        // cout << base_filename << endl;
    }else { cout << "No match found" << endl; }

    // Variables
    vector<vector<double>> columns;
    vector<double> height_data, temperature_data, wind_x_data, wind_y_data;
    double km;

    if (!Ideal_Conditions) {
        
        vector<double> wind_speed_data, wind_direction_data;

        // std::cout << "Reading the Radiosonde data" << endl;
        read_csv(columns, input_filepath);
        height_data = columns[0];
        temperature_data = columns[1];
        wind_speed_data = columns[2];
        wind_direction_data = columns[3];

        wind_x_data.resize(wind_speed_data.size());
        wind_y_data.resize(wind_speed_data.size());
        
        for (int i=0; i< wind_speed_data.size();i++){
            wind_x_data[i] = wind_speed_data[i] * cos(wind_direction_data[i] * (M_PI / 180));
            wind_y_data[i] = wind_speed_data[i] * sin(wind_direction_data[i] * (M_PI / 180));
        }

        if (With_spline){
            // cout << "Creating splines from the data" << endl;

            if (height_data.size() == temperature_data.size() && 
                height_data.size() == wind_x_data.size() && 
                height_data.size() == wind_y_data.size()) {

                build_spline(height_data, temperature_data, temperature_spline);
                build_spline(height_data, wind_x_data, wind_x_spline);
                build_spline(height_data, wind_y_data, wind_y_spline);

            } else {
                cout << "Error: Data arrays are of different sizes." << endl;
                return 1;
            }
        }

        wind_speed_data.clear();
        wind_direction_data.clear();

        km = *max_element(height_data.begin(), height_data.end()) / 1000;
        // cout << "km: " << km << endl;

    }else km = 30;

    vector<vector<double>> corner_val = {
        {0,0,0},         // centre
        // {-65,65,0},   // upper left corner
        // {65,65,0},   // upper right corner
        // {-65,-65,0},   // lower left corner
        // {65,-65,0}   // lower right corner
    };

    vector<vector<double>> edge_val = {
        // {0, 65, 0},  // upper edge
        // {65, 0, 0},  // right edge
        // {0, -65, 0},  // lower edge
        // {-65, 0, 0}  // left edge
    };

    string w_filename = output_filepath + base_filename + "_wf.csv";

    ofstream csvfile(w_filename);
    csvfile << "x,y,z,theta,time,phi\n";

    double h = 1; // the time calc of the wavefront here  1 sec gap between wavefronts

    double theta_min = -M_PI / 2;
    double theta_max = M_PI / 2;
    double theta_step = M_PI / 350;

    double phi_min = 0;
    double phi_max = M_PI;
    double phi_step = phi_step_deg * M_PI / 180;

    vector<double> theta_values = drange(theta_min,theta_max,theta_step); 
    vector<double> phi_values = drange(phi_min,phi_max,phi_step);

    for (vector<double>& arr : corner_val) {
        for (const double& phi : phi_values) {
            for (const double& theta : theta_values) {
                vector<double> y = { arr[0], arr[1], arr[2], theta };
                tuple<vector<vector<double>>, vector<double>, vector<double>> result = integrate_(y, km, h, phi, height_data, temperature_data, wind_x_data, wind_y_data);
                vector<vector<double>>& y_values = get<0>(result);
                vector<double>& t_values = get<1>(result);
                vector<double>& phi_values = get<2>(result);
                for (size_t i = 0; i < y_values.size(); ++i) {
                    csvfile << y_values[i][0] << "," << y_values[i][1] << "," << y_values[i][2] << ","
                            << y_values[i][3] << "," << t_values[i] << "," << phi_values[i] << "\n";
                }
            }
        }
    }
    
    csvfile.close();

    chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
    chrono::seconds duration = chrono::duration_cast<chrono::seconds>(end - start);
    // std::cout << "Time Taken: " << duration.count() << " seconds" << endl;
    return 0;
}

void build_spline(const vector<double>& x, const vector<double>& y, spline1dinterpolant& c) {
    // Convert vectors to ALGLIB's real_1d_array
    real_1d_array x_alglib, y_alglib;

    // Debugging: Check for NaN or Inf values
    for (size_t i = 0; i < x.size(); ++i) {
        if (isnan(x[i]) || isnan(y[i]) || isinf(x[i]) || isinf(y[i])) {
            cerr << "Error: NaN or Inf value found at index " << i << endl;
            return;
        }
    }

    // Create new vectors with edge cases handled
    vector<double> x_extended = x;
    vector<double> y_extended = y;

    if (zero_insertion){
        // Add edge points with zero values if necessary
        if (x[0] > 0) {
            x_extended.insert(x_extended.begin(), 0);
            y_extended.insert(y_extended.begin(), 0);
        }
        if (x.back() < 0) {
            x_extended.push_back(0);
            y_extended.push_back(0);
        }
    }
    
    x_alglib.setcontent(x_extended.size(), x_extended.data());
    y_alglib.setcontent(y_extended.size(), y_extended.data());

    ae_int_t n = x.size();

    try {
        if (!akima_spline){

        // Build the cubic spline interpolant with first derivative boundary conditions
        // Boundary condition types and values
        ae_int_t boundltype = 1;  // First derivative boundary condition
        double boundl = 0.0;      // Value of the first derivative at the left boundary
        ae_int_t boundrtype = 1;  // First derivative boundary condition
        double boundr = 0.0;      // Value of the first derivative at the right boundary
        
        // Build the cubic spline interpolant with first derivative boundary conditions
        cout << "Bulding Cubic Spline" << endl;
        spline1dbuildcubic(x_alglib, y_alglib, n, boundltype, boundl, boundrtype, boundr, c);
        }
        else {
            // Build the Akima spline interpolant
            // cout << "Bulding akima Spline" << endl;
            spline1dbuildakima(x_alglib, y_alglib, n,  c);
        }
    }
    catch (alglib::ap_error& e) {
        cerr << "ALGLIB error: " << e.msg << endl;
    }

}

int read_csv(vector<vector<double>>& columns, string& filename){
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file" << endl;
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
            } catch (const invalid_argument& e) {
                cerr << "Invalid input for column " << colIndex + 1 << ": " << cell << endl;
            } catch (const std::out_of_range& e) {
                cerr << "Input out of range for column " << colIndex + 1 << ": " << cell << endl;
            }
            colIndex++;
        }
    }

    file.close();
    return 0;
}

double Temperature_(double z) {
    if (With_temperature) {
        if (z <= Tropopause_Height) {
            return 300 + (Troposphere_Lapse_Rate / 1000) * z;
        } else if (z > Tropopause_Height && z <= Stratopause_Height) {
            return Temperature_(Tropopause_Height) + (Stratosphere_Lapse_Rate / 1000) * (z - Tropopause_Height);
        } else {
            return Temperature_(Stratopause_Height);
        }
    } else {
        return 292.804;
    }
}

double Speed_Of_Sound(double temperature) {
    if (With_temperature) {
        return sqrt(Gamma * R_air * temperature);
    } else {
        return SpeedOfSound;
    }
}

double avg_data(const vector<double>& height_data, const vector<double>& data, double target_height) {
    if (target_height < height_data.front() || target_height > height_data.back()) {
        cerr << "Target height out of bounds" << endl;
        return 0;
    }

    auto it = lower_bound(height_data.begin(), height_data.end(), target_height);
    if (it == height_data.end() || it == height_data.begin()) {
        cerr << "Error in finding bounding heights" << endl;
        return 0;
    }

    size_t index = it - height_data.begin();
    double x0 = height_data[index - 1];
    double y0 = data[index - 1];
    double x1 = height_data[index];
    double y1 = data[index];

    return y0 + (target_height - x0) * (y1 - y0) / (x1 - x0);
}

double forward_calc(int index, const vector<double>& height_data, const vector<double>& temp_data) {
    double slope = (temp_data[index + 1] - temp_data[index]) / (height_data[index + 1] - height_data[index]);
    return temp_data[index] + slope * (height_data[index] - height_data[index]);
}

double backward_calc(int index, const vector<double>& height_data, const vector<double>& temp_data) {
    double slope = (temp_data[index] - temp_data[index - 1]) / (height_data[index] - height_data[index - 1]);
    return temp_data[index] - slope * (height_data[index] - height_data[index - 1]);
}

double slope(const vector<double>& height_data, const vector<double>& data, double target_height) {
    auto it = lower_bound(height_data.begin(), height_data.end(), target_height);
    if (it == height_data.end() || it == height_data.begin()) {
        cerr << "Error in finding bounding heights for slope calculation" << endl;
        return 0;
    }

    size_t index = it - height_data.begin();
    double x0 = height_data[index - 1];
    double y0 = data[index - 1];
    double x1 = height_data[index];
    double y1 = data[index];

    return (y1 - y0) / (x1 - x0);
}

vector<double> ray_tracing_eqns(vector<double> y_, double phi, const vector<double>& height_data, const vector<double>& temperature_data, const vector<double>& wind_x_data, const vector<double>& wind_y_data) {
    double z = y_[2];
    double theta = y_[3];
    double temp, dc_dz, v_x, v_y, dvx_dz, dvy_dz;

    if (Ideal_Conditions) {
        if (With_temperature) {
            temp = Temperature_(z);
            if (temp <= 0) temp = 292.804;

            double lapse = 0;

            if (z <= Tropopause_Height) lapse = Troposphere_Lapse_Rate;
            else if (z > Tropopause_Height && z <= Stratopause_Height) lapse = Stratosphere_Lapse_Rate;

            dc_dz = sqrt((Gamma * R_air) / temp) * 0.5 * (lapse / 1000);

        } else {
            temp = 292.804;
            dc_dz = 0;
        }

        v_x = 0;
        v_y = 0;
        dvx_dz = 0;
        dvy_dz = 0;

    } else {
        if (With_temperature) {
            if (With_spline) {
                double ds, d2s;
                spline1ddiff(temperature_spline, z, temp, ds, d2s);
                if (temp <= 0) temp = 292.804;
                dc_dz = sqrt(Gamma * R_air / temp) * 0.5 * ds;
            } else {
                temp = avg_data(height_data, temperature_data, z);
                if (temp <= 0) temp = 292.804;
                dc_dz = sqrt(Gamma * R_air / temp) * 0.5 * slope(height_data, temperature_data, z);
            }

        } else {
            temp = 292.804;
            dc_dz = 0;
        }

        if (With_wind) {
            if (With_spline) {
                double ds, d2s;
                spline1ddiff(wind_x_spline, z, v_x, ds, d2s);
                dvx_dz = ds;

                spline1ddiff(wind_y_spline, z, v_y, ds, d2s);
                dvy_dz = ds;

            } else {
                v_x = avg_data(height_data, wind_x_data, z);
                v_y = avg_data(height_data, wind_y_data, z);
                dvx_dz = slope(height_data, wind_x_data, z);
                dvy_dz = slope(height_data, wind_y_data, z);
            }

        } else {
            v_x = 0;
            v_y = 0;
            dvx_dz = 0;
            dvy_dz = 0;
        }
    }

    double c_s = Speed_Of_Sound(temp);

    // PDE of Acoustic Wave Propagation
    double dx_dt = v_x + c_s * sin(theta) * cos(phi);
    double dy_dt = v_y + c_s * sin(theta) * sin(phi);
    double dz_dt = c_s * cos(theta);
    double dtheta_dt = sin(theta) * (dc_dz + sin(theta) * (cos(phi) * dvx_dz + sin(phi) * dvy_dz));

    return { dx_dt, dy_dt, dz_dt, dtheta_dt };
}


vector<double> Rk4Solver(vector<double> y, double phi, double h, const vector<double>& height_data, const vector<double>& temperature_data, const vector<double>& wind_x_data, const vector<double>& wind_y_data) {
    vector<double> k1 = ray_tracing_eqns(y, phi, height_data, temperature_data, wind_x_data, wind_y_data);

  
    vector<double> y_temp(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y_temp[i] = y[i] + (k1[i] * h) / 2;
    }
    vector<double> k2 = ray_tracing_eqns(y_temp, phi, height_data, temperature_data, wind_x_data, wind_y_data);

    for (size_t i = 0; i < y.size(); ++i) {
        y_temp[i] = y[i] + (k2[i] * h) / 2;
    }
    vector<double> k3 = ray_tracing_eqns(y_temp, phi, height_data, temperature_data, wind_x_data, wind_y_data);

    for (size_t i = 0; i < y.size(); ++i) {
        y_temp[i] = y[i] + k3[i] * h;
    }
    vector<double> k4 = ray_tracing_eqns(y_temp, phi, height_data, temperature_data, wind_x_data, wind_y_data);

    vector<double> y_value(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y_value[i] = y[i] + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }

    y_temp.clear();

    return y_value;
}


tuple<vector<vector<double>>, vector<double>, vector<double>> integrate_(vector<double> y, double km, double h, double phi, const vector<double>& height_data, const vector<double>& temperature_data, const vector<double>& wind_x_data, const vector<double>& wind_y_data) {
    double t0 = 0;
    double t_end = (km * 1000) / SpeedOfSound;
    int steps = (t_end - t0) / h ;

    vector<vector<double>> y_values(steps + 1, vector<double>(y.size()));
    y_values[0] = y;
    vector<double> t_values(steps + 1);
    t_values[0] = t0;
    vector<double> phi_values(steps + 1);
    phi_values[0] = phi;

    for (int i = 1; i <= steps; i++) {
        t_values[i] = t_values[i - 1] + h;
        y_values[i] = Rk4Solver(y_values[i - 1], phi, h, height_data, temperature_data, wind_x_data, wind_y_data);
        phi_values[i] = phi;
    }
    
    return { y_values, t_values, phi_values };
}

vector<double> drange(double x, double y, double jump) {
    vector<double> range;
    for (double val = x; val < y; val += jump) {
        range.push_back(val);
    }
    return range;
}

void write_to_csv(const string& filename, const vector<vector<double>>& data) {
    ofstream csvfile(filename);
    for (const auto& row : data) {
        for (const auto& value : row) {
            csvfile << value << ",";
        }
        csvfile << "\n";
    }
    csvfile.close();
}
