#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <regex>
#include <chrono>
#include <fstream>
#include "/home/murali/Documents/rass/data/cpp_lib/alglib-cpp/src/interpolation.h"
#include "/home/murali/Documents/rass/data/cpp_lib/alglib-cpp/src/statistics.h"
#include "/home/murali/Documents/rass/data/cpp_lib/matplotlibcpp.h"
#include "/home/murali/Documents/rass/data/cpp_lib/csv.h"  


namespace plt = matplotlibcpp;
using namespace alglib;
using namespace std;
using namespace std::chrono;

struct DataFrame {
    vector<double> time;
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> phi;
    vector<double> theta;
};

DataFrame readCSV(const string &filepath) {
    DataFrame df;
    io::CSVReader<6> in(filepath);
    in.read_header(io::ignore_extra_column, "time", "x", "y", "z", "phi", "theta");

    double time, x, y, z, phi, theta;
    while (in.read_row(time, x, y, z, phi, theta)) {
        df.time.push_back(time);
        df.x.push_back(x);
        df.y.push_back(y);
        df.z.push_back(z);
        df.phi.push_back(phi);
        df.theta.push_back(theta);
    }
    return df;
}

DataFrame filterDataFrame(const DataFrame &df, const string &dist_axis) {
    DataFrame filtered_df;
    if (dist_axis == "x") {
        for (size_t i = 0; i < df.time.size(); i++) {
            if (df.phi[i] == 0 && df.z[i] >= 0 && df.theta[i] <= M_PI / 2 && df.theta[i] >= -M_PI / 2) {
                filtered_df.time.push_back(df.time[i]);
                filtered_df.x.push_back(df.x[i]);
                filtered_df.y.push_back(df.y[i]);
                filtered_df.z.push_back(df.z[i]);
                filtered_df.phi.push_back(df.phi[i]);
                filtered_df.theta.push_back(df.theta[i]);
            }
        }
    } else {
        for (size_t i = 0; i < df.time.size(); i++) {
            if (df.phi[i] != 0 && df.z[i] >= 0 && df.theta[i] <= M_PI / 2 && df.theta[i] >= -M_PI / 2) {
                filtered_df.time.push_back(df.time[i]);
                filtered_df.x.push_back(df.x[i]);
                filtered_df.y.push_back(df.y[i]);
                filtered_df.z.push_back(df.z[i]);
                filtered_df.phi.push_back(df.phi[i]);
                filtered_df.theta.push_back(df.theta[i]);
            }
        }
    }
    return filtered_df;
}

vector<double> landPoints(const vector<double> &dist_axis, const vector<double> &altitude_axis) {
    vector<double> land_pts;
    real_1d_array x, y;
    x.setcontent(dist_axis.size(), &dist_axis[0]);
    y.setcontent(altitude_axis.size(), &altitude_axis[0]);

    spline1dinterpolant spline;
    spline1dbuildcubic(x, y, spline);

    for (size_t i = 0; i < dist_axis.size(); i++) {
        double x0 = dist_axis[i];
        double dz_dx = spline1ddiff1(spline, x0);
        if (dz_dx != 0) {
            double n_slope = atan(-1 / dz_dx) * 180 / M_PI;
            if (abs(n_slope) >= 70) {
                double land_x = x0 + spline1dcalc(spline, x0) * dz_dx;
                land_pts.push_back(land_x);
            }
        }
    }
    return land_pts;
}

vector<pair<double, double>> wavePoints(const DataFrame &df, const string &dist_axis, const string &altitude_axis, double center) {
    vector<pair<double, double>> points_l, points_w;

    // Adjusted data
    vector<double> adj_dist_axis = (dist_axis == "x") ? df.x : df.y;
    for (auto &d : adj_dist_axis) d -= center;

    vector<double> unique_times = df.time;
    sort(unique_times.begin(), unique_times.end());
    unique_times.erase(unique(unique_times.begin(), unique_times.end()), unique_times.end());

    for (auto time_slice : unique_times) {
        if (time_slice != 0) {
            // Filter by time slice and sort by dist_axis
            DataFrame filtered_df;
            for (size_t i = 0; i < df.time.size(); i++) {
                if (df.time[i] == time_slice) {
                    filtered_df.time.push_back(df.time[i]);
                    filtered_df.x.push_back(df.x[i]);
                    filtered_df.y.push_back(df.y[i]);
                    filtered_df.z.push_back(df.z[i]);
                    filtered_df.phi.push_back(df.phi[i]);
                    filtered_df.theta.push_back(df.theta[i]);
                }
            }

            // Remove duplicates
            auto last = unique(filtered_df.x.begin(), filtered_df.x.end());
            filtered_df.x.erase(last, filtered_df.x.end());

            // Cubic spline interpolation
            real_1d_array x, y;
            x.setcontent(filtered_df.x.size(), &filtered_df.x[0]);
            y.setcontent(filtered_df.z.size(), &filtered_df.z[0]);
            spline1dinterpolant spline;
            spline1dbuildcubic(x, y, spline);

            for (auto x0 : filtered_df.x) {
                double dz_dx = spline1ddiff1(spline, x0);
                if (dz_dx != 0) {
                    double n_slope = atan(-1 / dz_dx) * 180 / M_PI;
                    if (abs(n_slope) >= 70) {
                        double land_x = x0 + spline1dcalc(spline, x0) * dz_dx;
                        points_w.push_back({x0, spline1dcalc(spline, x0)});
                        if (abs(land_x) < 65) {
                            points_l.push_back({x0, spline1dcalc(spline, x0)});
                        }
                    }
                }
            }
        }
    }

    // Remove duplicates
    sort(points_w.begin(), points_w.end());
    points_w.erase(unique(points_w.begin(), points_w.end()), points_w.end());

    sort(points_l.begin(), points_l.end());
    points_l.erase(unique(points_l.begin(), points_l.end()), points_l.end());

    return {points_w, points_l};
}

double findNearestLocation(const vector<double> &locations, double current_loc) {
    auto min_it = min_element(locations.begin(), locations.end(), [current_loc](double a, double b) {
        return abs(a - current_loc) < abs(b - current_loc);
    });
    return *min_it;
}

size_t countPointsInWindow(double start, double end, const vector<double> &sorted_array) {
    return count_if(sorted_array.begin(), sorted_array.end(), [start, end](double x) {
        return x >= start && x <= end;
    });
}

pair<vector<pair<double, double>>, vector<double>> windowPoints(const vector<double> &xi, double window_size, int top_no_of_windows) {
    vector<pair<double, double>> top_windows;
    vector<double> midpt_windows;
    vector<double> x_l_sorted = xi;
    sort(x_l_sorted.begin(), x_l_sorted.end());

    for (size_t i = 0; i < x_l_sorted.size(); i++) {
        double window_start = x_l_sorted[i];
        double window_end = window_start + window_size;
        size_t count_points = countPointsInWindow(window_start, window_end, x_l_sorted);

        if (top_windows.size() < top_no_of_windows || count_points > top_windows.back().second) {
            top_windows.push_back({window_start, window_end});
            midpt_windows.push_back((window_start + window_end) / 2);
        }
    }

    // Sort and retain top windows
    vector<pair<pair<double, double>, double>> windows_with_counts;
    for (size_t i = 0; i < top_windows.size(); i++) {
        windows_with_counts.push_back({top_windows[i], countPointsInWindow(top_windows[i].first, top_windows[i].second, x_l_sorted)});
    }
    sort(windows_with_counts.begin(), windows_with_counts.end(), [](auto &a, auto &b) {
        return a.second > b.second;
    });
    if (windows_with_counts.size() > top_no_of_windows) {
        windows_with_counts.resize(top_no_of_windows);
    }

    top_windows.clear();
    midpt_windows.clear();
    for (auto &wc : windows_with_counts) {
        top_windows.push_back(wc.first);
        midpt_windows.push_back((wc.first.first + wc.first.second) / 2);
    }

    return {top_windows, midpt_windows};
}

double angleWithXAxis(double x, double y) {
    return atan2(y, x) * 180 / M_PI;
}

vector<pair<double, double>> pointsInCone(const vector<pair<double, double>> &points, double angle_start, double angle_width) {
    double angle_end = angle_start + angle_width;
    vector<pair<double, double>> points_in_cone;

    for (auto &point : points) {
        double x = point.first;
        double y = point.second;
        double angle = angleWithXAxis(x, y);
        if (angle < 0) angle += 360;
        if (angle >= angle_start && angle <= angle_end) {
            points_in_cone.push_back(point);
        }
    }
    return points_in_cone;
}

void processWaveformFile(const string &wf_file_path, const string &dist_axis) {
    auto start = high_resolution_clock::now();

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

    double eps = 5;
    int min_samples = 3;
    vector<double> cluster_centers;

    real_2d_array xy;
    xy.setlength(acoustic_source_loc.size(), 1);
    for (size_t i = 0; i < acoustic_source_loc.size(); i++) {
        xy[i][0] = acoustic_source_loc[i];
    }

    clusterizerstate s;
    clusterizercreate(s);
    clusterizersetpoints(s, xy, 2);
    clusterizersetkmeanslimits(s, 5, 0);
    clusterizerreport rep;
    kmeansgenerate(s, min_samples, eps, 5, rep);

    for (int i = 0; i < rep.cidx.length(); i++) {
        if (rep.cidx[i] != -1) {
            cluster_centers.push_back(rep.c[i][0]);
        }
    }

    cout << "Cluster centers: ";
    for (auto center : cluster_centers) {
        cout << center << " ";
    }
    cout << endl;

    double angle_width = 3;
    double angle_start = 68;
    double angle_end = 113 - angle_width;
    double angle_step = 1;

    vector<vector<pair<double, double>>> points_w, points_l;

    for (auto center : cluster_centers) {
        auto [w, l] = wavePoints(dfa_ew, dist_axis, "z", center);
        points_w.push_back(w);
        points_l.push_back(l);
    }

    double angle_f = 0, percent_pts_f = 0;

    for (double angle = angle_start; angle < angle_end; angle += angle_step) {
        vector<pair<double, double>> points_within_cone_w;
        vector<pair<double, double>> points_within_cone_l;

        for (auto &pw : points_w) {
            auto pw_in_cone = pointsInCone(pw, angle, angle_width);
            points_within_cone_w.insert(points_within_cone_w.end(), pw_in_cone.begin(), pw_in_cone.end());
        }

        for (auto &pl : points_l) {
            auto pl_in_cone = pointsInCone(pl, angle, angle_width);
            points_within_cone_l.insert(points_within_cone_l.end(), pl_in_cone.begin(), pl_in_cone.end());
        }

        double percent_pts = (points_within_cone_l.size() / static_cast<double>(points_within_cone_w.size())) * 100;
        if (percent_pts_f <= percent_pts) {
            angle_f = angle;
            percent_pts_f = percent_pts;
        }
    }

    cout << "Percentage inside the cone: " << percent_pts_f << ", angle between: " << angle_f << " - " << angle_f + angle_width << endl;

    vector<pair<double, double>> points_w_flat;
    vector<pair<double, double>> points_l_flat;
    for (auto &pw : points_w) {
        points_w_flat.insert(points_w_flat.end(), pw.begin(), pw.end());
    }
    for (auto &pl : points_l) {
        points_l_flat.insert(points_l_flat.end(), pl.begin(), pl.end());
    }

    vector<double> px_w, py_w, px_l, py_l;
    for (auto &p : points_w_flat) {
        px_w.push_back(p.first);
        py_w.push_back(p.second);
    }
    for (auto &p : points_l_flat) {
        px_l.push_back(p.first);
        py_l.push_back(p.second);
    }

    plt::scatter(px_w, py_w, 0.9, {{"label", "Wave Points"}});
    plt::scatter(px_l, py_l, 0.9, {{"color", "#A05544"}, {"label", "Land Points"}});

    double dist = 35 * 1000;
    plt::xlim(-dist, dist);
    plt::ylim(0, dist);
    plt::xlabel("Distance from source (m)");
    plt::ylabel("Altitude (m)");
    plt::title(dist_axis == "x" ? "N-S" : "E-W");
    plt::grid(true);

    double theta1 = angle_f * M_PI / 180;
    double theta2 = (angle_f + angle_width) * M_PI / 180;

    plt::plot({0, cos(theta1) * dist}, {0, sin(theta1) * dist}, "g--", {{"label", "Cone-section"}});
    plt::plot({0, cos(theta2) * dist}, {0, sin(theta2) * dist}, "g--");

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end - start).count();
    cout << "Time taken: " << duration << " seconds" << endl;

    plt::legend();
    plt::save("/home/murali/Documents/rass/automation/gen_data/" + to_string(year) + "-" + to_string(month) + "-" + to_string(day) + "_" + (dist_axis == "x" ? "N-S" : "E-W") + ".png");

    ofstream centers_file("/home/murali/Documents/rass/automation/gen_data/acoustic_sources.csv", ios_base::app);
    centers_file << "Cluster_Centers,Axis,date\n";
    for (auto center : cluster_centers) {
        centers_file << center << "," << (dist_axis == "x" ? "x" : "y") << "," << year << "-" << month << "-" << day << "\n";
    }
    centers_file.close();
}

int main(int argc, char **argv) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <file_path> <dist_axis>" << endl;
        return 1;
    }

    string file_path = argv[1];
    string dist_axis = argv[2];

    processWaveformFile(file_path, dist_axis);

    return 0;
}
