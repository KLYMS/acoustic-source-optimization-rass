#include<bits/stdc++.h>
#include "spline.h"

using namespace std;

const double Gamma = 1.4;
const double R_gas = 287;
const double Troposphere_Lapse_Rate = -6.5;
const double Tropopause_Height = 18000; 
const double Stratosphere_Lapse_Rate = 2.5;
const double Stratopause_Height = 48000;
/* const double Stratopause_Temperature = 265; */
const double Ref_Wind_speed = 250.5;
const double Ref_height = 28000;
const double Alpha= 0.143;
bool With_temperature = true;
bool With_wind = true;
const double SpeedOfSound = 343;
const double phi = 0;


int read_csv(vector<vector<double>>& columns, string& filename);
vector<vector<double>> integrate(vector<double> y0,double km, double h);

tk::spline height_temp, height_wind;

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


/* struct Spline { */
/*     double a, b, c, d, x; */
/* }; */

/* class CubicSpline { */
/* public: */
/*     CubicSpline(const std::vector<double>& x, const std::vector<double>& y) { */
/*         int n = x.size(); */
/*         if (n < 2) throw std::invalid_argument("There must be at least two points."); */

/*         splines.resize(n); */

/*         for (int i = 0; i < n; ++i) { */
/*             splines[i].x = x[i]; */
/*             splines[i].a = y[i]; */
/*         } */

/*         vector<double> alpha(n - 1), l(n), mu(n), z(n); */
/*         vector<double> h(n - 1); */

/*         for (int i = 0; i < n - 1; ++i) { */
/*             h[i] = x[i + 1] - x[i]; */
/*             if (h[i] == 0.0) throw invalid_argument("Consecutive x values must be different."); */
/*         } */

/*         for (int i = 1; i < n - 1; ++i) { */
/*             alpha[i] = 3.0 / h[i] * (splines[i + 1].a - splines[i].a) - 3.0 / h[i - 1] * (splines[i].a - splines[i - 1].a); */
/*         } */

/*         l[0] = 1.0; */
/*         mu[0] = z[0] = 0.0; */

/*         for (int i = 1; i < n - 1; ++i) { */
/*             l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]; */
/*             mu[i] = h[i] / l[i]; */
/*             z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]; */
/*         } */

/*         l[n - 1] = 1.0; */
/*         z[n - 1] = 0.0; */
/*         splines[n - 1].c = 0.0; */

/*         for (int j = n - 2; j >= 0; --j) { */
/*             splines[j].c = z[j] - mu[j] * splines[j + 1].c; */
/*             splines[j].b = (splines[j + 1].a - splines[j].a) / h[j] - h[j] * (splines[j + 1].c + 2.0 * splines[j].c) / 3.0; */
/*             splines[j].d = (splines[j + 1].c - splines[j].c) / (3.0 * h[j]); */
/*         } */
/*     } */

/*     double interpolate(double x) const { */
/*         if (splines.empty()) throw logic_error("Spline coefficients are not computed."); */

/*         const Spline* s = nullptr; */
/*         if (x <= splines.front().x) { */
/*             s = &splines.front(); */
/*         } else if (x >= splines.back().x) { */
/*             s = &splines.back(); */
/*         } else { */
/*             for (size_t i = 1; i < splines.size(); ++i) { */
/*                 if (x < splines[i].x) { */
/*                     s = &splines[i - 1]; */
/*                     break; */
/*                 } */
/*             } */
/*         } */

/*         double dx = x - s->x; */
/*         return s->a + s->b * dx + s->c * dx * dx + s->d * dx * dx * dx; */
/*     } */

/* private: */
/*     vector<Spline> splines; */
/* }; */

class ODESystem {
    public:
        ODESystem(double phi,const vector<double>& heights, const vector<double>& temps, const vector<double>& winds) 
            :  phi(phi) {        
        
                if (With_temperature) {
                    height_temp.set_points(heights, temps);
                }
                if (With_wind) {
                    height_wind.set_points(heights, winds);
                }
            }
        
        vector<double> operator()(double t, vector<double>& y){
        double x_coord = y[0];
        double y_coord = y[1];
        double z_coord = y[2];
        double theta = y[3];

        double dx_dt = v_x(z_coord) + c_s(z_coord) * sin(theta) * cos(phi);
        double dy_dt = v_y(z_coord) + c_s(z_coord) * sin(theta) * sin(phi);
        double dz_dt = v_z(z_coord) + c_s(z_coord) * cos(theta);
        double dtheta_dt = sin(theta) * (partial_cs_z(z_coord) + cos(phi) * partial_vx_z(z_coord) + sin(phi) * partial_vy_z(z_coord));

        return {dx_dt, dy_dt, dz_dt, dtheta_dt};
    }
        
    private:
        double phi;
                
        /* double temperature(double z) { */
        /*     if(With_temperature){ */
        /*         if (z <= Tropopause_Height) return 300 + (Troposphere_Lapse_Rate / 1000) * z; */ 
        /*         else if(z > Tropopause_Height && z<= Stratopause_Height ) return temperature(Tropopause_Height) + (Stratosphere_Lapse_Rate / 1000) * (z - Tropopause_Height); */ 
        /*         else return temperature(Stratopause_Height); */
        /*     }else return 0; */
        /* } */

        double v_x(double z){
            if(With_wind){
                return height_wind(z);
            }else return 0;
        }
      
        double v_y(double z){
            return 0;
        }
       
        double v_z(double z){
            return 0;
        }

        double c_s(double z){
            if(With_temperature){
                double t = height_temp(z);
                double vel = sqrt(Gamma * R_gas * t);
                return vel;
            }
            else return SpeedOfSound;
        }

        double partial_cs_z(double z){
            double n = sqrt((Gamma * R_gas)/(4 * height_temp(z)));
            return n * height_temp.deriv(1,z);
            /* if(With_temperature){ */
            /*     double lapse; */
            /*     if(z <= Tropopause_Height)  lapse = Troposphere_Lapse_Rate; */
            /*     else if(z >= Tropopause_Height && z <= Stratopause_Height) lapse = Stratosphere_Lapse_Rate; */
            /*     else lapse = 0; */
                
            /*     double num = sqrt(Gamma * R_gas) * (lapse/1000); */
            /*     double t = temperature(z); */
            /*     return  num/(2 * sqrt(t)); */
            /* } */
            /* else return 0.0; */
        }

        double partial_vx_z(double z){
            if(With_wind){
               return height_wind.deriv(1,z);
            }
            else return 0;
        }

        double partial_vy_z(double z){
            if(With_wind){
                return 0.0;
            }
            else return 0.0;
        }
};



class RK4Solver {
    public:
        RK4Solver(double h) : h(h) {}
        vector<double> step( ODESystem& system, double t, vector<double>& y){
            vector<double> k1,k2,k3,k4, y_temp;

            k1 = system(t,y);
            y_temp = y;
            for(int i=0; i< y.size();i++) y_temp[i] = y[i] + 0.5 * k1[i];
            
            k2 = system(t+0.5*h,y_temp);
            y_temp = y;
            for(int i=0; i< y.size();i++) y_temp[i] = y[i] + 0.5 * k2[i];

            k3 = system(t + 0.5 * h, y_temp);
            y_temp = y;
            for(int i=0; i< y.size();i++) y_temp[i] = y[i] + k3[i];
            k4 = system(t+h,y_temp);

            vector<double> y_next(y.size());
            for(int i=0;i<y.size();i++){
                y_next[i] = y[i] + (h/6) * (k1[i]+ 2 * (k2[i]+k3[i]) + k4[i]);
            }

            return y_next;

        }

    private:
        double h;
};

vector<vector<double>> integrate(vector<double> y0,double km, double h, vector<double>& heights, vector<double>& temps, vector<double> winds){    
    
    double t0 = 0.0;
    double tf = (km * 1000)/SpeedOfSound;
    int steps = (tf - t0)/h;
    ODESystem system(phi,heights,temps,winds);
    RK4Solver solver(h);
    
    vector<vector<double>> y_values(steps, vector<double>( y0.size() ));
    y_values[0] = y0;
    for(int i=1;i<steps;i++){
        double t = t0 + (i-1) * h;
        y_values[i] = solver.step(system,t,y_values[i-1]);
    }

    vector<double> final_values =  y_values.back();
    cout << "Final x: " << final_values[0]
              << ", y: " << final_values[1]
              << ", z: " << final_values[2]
              << ", theta: " << final_values[3] << endl;

    cout << "steps: " << steps << endl;
    cout << endl;

    return y_values;

   
}

int main(){
    double km=30;
    double h = 5;
    double theta_min = -M_PI / 2;
    double theta_max = M_PI / 2;
    double theta_step = M_PI / 300;
    
    vector<vector<double>> data;
    string filename = "wind_radio_data.csv";
    read_csv(data,filename);
    
    vector<double> heights = data[0];
    vector<double> temps = data[1];
    vector<double> winds = data[2];
   
    ofstream outputFile;
    if (With_temperature && With_wind) outputFile.open("wTemp_wWind.csv");
    else if (With_temperature) outputFile.open("ideal_wTemp.csv");
    else if (With_wind) outputFile.open("ideal_wWind.csv");
    else outputFile.open("ideal_case.csv");

    if (!outputFile.is_open()) {
        cout << "Failed to open the file." << endl;
        return 1;
    }
    
    outputFile << "x,y,z,theta" << endl;    
    
    for(double i = theta_min; i<=theta_max;i+=theta_step){
        vector<double> y = {0,0,0,i};
        vector<vector<double>> y_values = integrate(y,km,h,heights,temps,winds);        
        for (const auto& row : y_values){
            for (size_t i = 0; i < row.size(); i++) {
                outputFile << row[i];
                if (i < row.size() - 1) outputFile << ",";
            }
            outputFile << "\n";
        }
    }
    outputFile.close();

    return 0;
}


