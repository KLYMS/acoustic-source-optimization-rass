#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <regex>
#include "/home/murali/Documents/rass/data/cpp_lib/alglib-cpp/src/interpolation.h"

using namespace std;
using namespace std::chrono;
using namespace alglib;

int read_csv(vector<vector<double>>& columns, string& filename);
void processWaveformFile(const string &wf_file_path, const string &dist_axis);

int main(){
        
}

void processWaveformFile(const string &wf_file_path, const string &dist_axis) {

    // Extract date from the filename
    regex date_pattern(R"((\d{4})-(\d{2})-(\d{2}))");
    smatch match_date;
    regex_search(wf_file_path, match_date, date_pattern);
    int year = match_date.empty() ? 0 : stoi(match_date[1].str());
    int month = match_date.empty() ? 0 : stoi(match_date[2].str());
    int day = match_date.empty() ? 0 : stoi(match_date[3].str());

    cout << "Year: " << year << ", Month: " << month << ", Day: " << day << endl;

    DataFrame dfa = readCSV(wf_file_path);
    DataFrame dfa_ew = filterDataFrame(dfa, dist_axis);

    double window_size = 130;
    int top_no_of_windows = 25;

    double initial_loc = 0;
    vector<double> acoustic_source_loc;

    vector<double> unique_times = dfa_ew.time;
    sort(unique_times.begin(), unique_times.end());
    unique_times.erase(unique(unique_times.begin(), unique_times.end()), unique_times.end());

    for (auto time_slice : unique_times) {
        if (time_slice != 0) {
            DataFrame filtered_df;
            for (size_t i = 0; i < dfa_ew.time.size(); i++) {
                if (dfa_ew.time[i] == time_slice) {
                    filtered_df.time.push_back(dfa_ew.time[i]);
                    filtered_df.x.push_back(dfa_ew.x[i]);
                    filtered_df.y.push_back(dfa_ew.y[i]);
                    filtered_df.z.push_back(dfa_ew.z[i]);
                    filtered_df.phi.push_back(dfa_ew.phi[i]);
                    filtered_df.theta.push_back(dfa_ew.theta[i]);
                }
            }

            // Remove duplicates
            auto last = unique(filtered_df.x.begin(), filtered_df.x.end());
            filtered_df.x.erase(last, filtered_df.x.end());

            vector<double> land_pts = landPoints(filtered_df.x, filtered_df.z);

            if (!land_pts.empty()) {
                auto [top_windows, midpt_windows] = windowPoints(land_pts, window_size, top_no_of_windows);
                double nearest_loc = findNearestLocation(midpt_windows, initial_loc);
                acoustic_source_loc.push_back(nearest_loc);
            }
        }
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