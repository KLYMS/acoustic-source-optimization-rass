#include<bits/stdc++.h>
using namespace std;

const double Gamma = 1.4;
const double R_air = 287;
const double Troposphere_Lapse_Rate = -6.5;
const double Tropopause_Height = 18000; 
const double Stratosphere_Lapse_Rate = 2.5;
const double Stratopause_Height = 48000;
/* const double Stratopause_Temperature = 265; */
const double Ref_Wind_speed = 250.5;
const double Ref_height = 28000;
const double Alpha= 0.143;

/* wind constants */ 
const double Boundary_layer_height = 2300;
const double kappa = 0.4; // von Kármán constant
const double u_star = 0.3; // Friction velocity for boundary layer
const double z0 = 0.1; // Roughness length for boundary layer

const double u_g = 10.0; // Geostrophic wind speed for troposphere
const double D = 1000.0; // Ekman depth for troposphere
const double wphi = 0.0; // Phase shift for troposphere

const double R_dry = 287.0; // Gas constant for dry air
const double f = 1e-4; // Coriolis parameter
const double dTdy = 0.01; // Meridional temperature gradient for stratosphere

class ODESystem {
    public:
        ODESystem(double phi) 
            : phi(phi) {}        
        
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
                
        double temperature(double z) {
            if (z <= Tropopause_Height) return 300 + (Troposphere_Lapse_Rate / 1000) * z; 
            else if(z > Tropopause_Height && z<= Stratopause_Height ) return  (300 + (Troposphere_Lapse_Rate / 1000) * Tropopause_Height) + (Stratosphere_Lapse_Rate / 1000) * (z - Tropopause_Height); 
            else return  (300 + (Troposphere_Lapse_Rate / 1000) * Tropopause_Height) + (Stratosphere_Lapse_Rate / 1000) * (Stratopause_Height - Tropopause_Height);
        }

        double Wind_speed(double z){
            if(z <= Boundary_layer_height && z > 0) return (u_star * log(z/z0))/kappa;
            else if(z > Boundary_layer_height && z <= Tropopause_Height) return u_g * (1- (exp(-z/D) * cos(z/D + wphi)));
            else if(z > Tropopause_Height ) return (-R_dry * dTdy * z)/f;
            else return 0;
        }
        
        // wind velovity
        double v_x(double z){
            return Wind_speed(z);
        }

        double v_y(double z){
            return 0;
        }

        double v_z(double z){
            return 0;
        }

        // sound velocity
        double c_s(double z){
            double t = temperature(z);
            double vel = sqrt(Gamma * R_air * t);
            return vel;
        }

        // partial derivatives
        double partial_cs_z(double z){
            double lapse;
            if(z <= Tropopause_Height)  lapse = Troposphere_Lapse_Rate;
            else if(z >= Tropopause_Height && z <= Stratopause_Height) lapse = Stratosphere_Lapse_Rate;
            else lapse = 0;
                
            double num = sqrt(Gamma * R_air) * (lapse/1000);
            double t = temperature(z);
            return  num/(2 * sqrt(t));
        }

        double partial_vx_z(double z){
            if(z <= Boundary_layer_height  && z > 0) return u_star / (z * kappa);
            else if(z <= Tropopause_Height && z > Boundary_layer_height){
                double e1 = (u_g * exp(- z/D))/D;
                double e2 = cos(z/D + wphi) + sin(z/D + wphi);
                return e1 * e2;
            }
            else if(z > Tropopause_Height ) return  (-R_dry * dTdy) / f;
            else return 0;
        }

        double partial_vy_z(double z){
            return 0.0;
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

vector<vector<double>> integrate(vector<double> y0,double km, double h){    
    
    double c_s = 343;
    double t0 = 0.0;
    double tf = (km * 1000)/c_s;
    double phi = 0;
    int steps = (tf - t0)/h;
    ODESystem system(phi);
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
    double h = 0.01;
    double theta_min = -M_PI / 2;
    double theta_max = M_PI / 2;
    double theta_step = M_PI / 400;

    ofstream outputFile("ideal_wWind.csv");
    if (!outputFile.is_open()) {
        cout << "Failed to open the file." << endl;
        return 1;
    } 
    outputFile << "x,y,z,theta" << endl;

    for(double i = theta_min; i<=theta_max;i+=theta_step){
        vector<double> y = {0,0,0,i};
        vector<vector<double>> y_values = integrate(y,km,h);        
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


