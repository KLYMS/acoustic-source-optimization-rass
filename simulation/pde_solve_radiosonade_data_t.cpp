#include<bits/stdc++.h>
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
bool With_wind = false;
const double SpeedOfSound = 343;
const double phi = 0;
vector<double> height_data;
vector<double> temp_data;
 


double avg_data(double target_height){

    auto it = find(height_data.begin(), height_data.end(), target_height);
    
    if (it != height_data.end()) {
        int index = distance(height_data.begin(), it);
        return temp_data[index];
    } else {
        auto insert_pos = lower_bound(height_data.begin(), height_data.end(), target_height);
        int index = distance(height_data.begin(), insert_pos);
        double interpolated_data;
        if (index == 0) {
            interpolated_data = temp_data[0];
        } else if (index == height_data.size()) {
            interpolated_data = temp_data.back();
        } else {
            int h1 = height_data[index - 1];
            int h2 = height_data[index];
            double t1 = temp_data[index - 1];
            double t2 = temp_data[index];
            double slope_ =  (t2 - t1) / (h2 - h1);
            interpolated_data  = t1 + ((target_height - h1) * slope_);
        }
        return interpolated_data;       
        /* double average_temp; */
        /* if (index == 0) { */
        /*     average_temp = temp_data[0]; */
        /* } else if (index == height_data.size()) { */
        /*     average_temp = temp_data.back(); */
        /* } else { */
        /*     average_temp = (temp_data[index - 1] + temp_data[index]) / 2.0; */
        /* } */
        
        /* return average_temp; */
    }

}

double forward_calc(int index, const std::vector<double>& height_data, const std::vector<double>& temp_data){
        double h1 = height_data[index];
        double h2 = height_data[index + 1];
        double t1 = temp_data[index];
        double t2 = temp_data[index + 1];
        return (t2 - t1) / (h2 - h1);
}

double backward_calc(int index, const std::vector<double>& height_data, const std::vector<double>& temp_data){
        double h1 = height_data[index - 1];
        double h2 = height_data[index];
        double t1 = temp_data[index - 1];
        double t2 = temp_data[index];
        return (t2 - t1) / (h2 - h1);
}

double slope(const std::vector<double>& height_data, const std::vector<double>& temp_data, double target_height) {
    if (height_data.size() != temp_data.size()) {
        throw invalid_argument("height_data and temp_data must have the same size");
    }

    auto it = find(height_data.begin(), height_data.end(), target_height);

    if (it != height_data.end()) {
        int index = distance(height_data.begin(), it);
        if (index < height_data.size() - 1) {
            return forward_calc(index,height_data,temp_data);
        
        }else if(index == height_data.size()-1){
            return backward_calc(index, height_data, temp_data);
        
        }else {
            throw out_of_range("Target height is more then the last element in height_data, slope calculation is not possible.");
        }

    } else {
        auto insert_pos = lower_bound(height_data.begin(), height_data.end(), target_height);
        int index = distance(height_data.begin(), insert_pos);

        if (index >= height_data.size()-1) {
            index = height_data.size()-1;
            return backward_calc(index, height_data,temp_data);
        }
        return backward_calc(index,height_data,temp_data);
    }
}


class ODESystem {
    public:
        ODESystem(double phi) 
            :  phi(phi) {}        
        
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
                
        double temperature(double z){
            if(With_temperature){
                if (z <= Tropopause_Height) return 300 + (Troposphere_Lapse_Rate / 1000) * z; 
                else if(z > Tropopause_Height && z<= Stratopause_Height ) return temperature(Tropopause_Height) + (Stratosphere_Lapse_Rate / 1000) * (z - Tropopause_Height); 
                else return temperature(Stratopause_Height);
            }else return 0;
        }

        
        double v_x(double z){
            return 0;
        }
      
        double v_y(double z){
            return 0;
        }
       
        double v_z(double z){
            return 0;
        }

        double c_s(double z){
            if(With_temperature){
                double t = avg_data(z);
                double vel = sqrt(Gamma * R_gas * t);
                return vel;
            }
            else return SpeedOfSound;
        }

        double partial_cs_z(double z){
            if(With_temperature){
                double n = sqrt((Gamma * R_gas)/(4 * avg_data(z)));
                return n * slope(height_data,temp_data,z);
                /* double lapse; */
                /* if(z <= Tropopause_Height)  lapse = Troposphere_Lapse_Rate; */
                /* else if(z >= Tropopause_Height && z <= Stratopause_Height) lapse = Stratosphere_Lapse_Rate; */
                /* else lapse = 0; */
                
                /* double num = sqrt(Gamma * R_gas) * (lapse/1000); */
                /* double t = temperature(z); */
                /* return  num/(2 * sqrt(t)); */
            }
            else return 0.0;
        }

        double partial_vx_z(double z){
            if(With_wind){
                return 0;
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

vector<vector<double>> integrate(vector<double> y0,double km, double h){    
    
    double t0 = 0.0;
    double tf = (km * 1000)/SpeedOfSound;
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

int read_csv(vector<double>& column1, vector<double>& column2){

    ifstream file("temp_radio_data.csv");
    if (!file.is_open()) {
        cout << "Error opening file" << endl;
        return 1;
    }

    string line;
    getline(file, line);

    while (std::getline(file, line)) {
        stringstream lineStream(line);
        string cell;

        if (getline(lineStream, cell, ',')) {
            try {
                cell.erase(remove_if(cell.begin(), cell.end(), ::isspace), cell.end());
                double value1 = stod(cell);
                column1.push_back(value1);
            } catch (const invalid_argument& e) {
                cerr << "Invalid input for column 1: " << cell << endl;
            } catch (const out_of_range& e) {
                cerr << "Input out of range for column 1: " << cell << endl;
            }
        }
        if (getline(lineStream, cell, ',')) {
            try {
                cell.erase(remove_if(cell.begin(), cell.end(), ::isspace), cell.end());
                double value2 = stod(cell);
                column2.push_back(value2);
            } catch (const invalid_argument& e) {
                cerr << "Invalid input for column 2: " << cell << endl;
            } catch (const out_of_range& e) {
                cerr << "Input out of range for column 2: " << cell << endl;
            }
        }
    }

    file.close();

    return 0;
}




int main(){
    double km=30;
    double h = 5;
    double theta_min = -M_PI / 2;
    double theta_max = M_PI / 2;
    double theta_step = M_PI / 300;
    
    read_csv(height_data,temp_data);

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


