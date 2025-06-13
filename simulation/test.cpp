#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

int main() {
    // Define vectors to hold the data from each column
    std::vector<double> column1;
    std::vector<double> column2;
    std::vector<double> column3;

    // Open the CSV file
    std::ifstream file("radiosonade_data_w.csv");
    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }

    std::string line;
    // Skip the header line
    std::getline(file, line);

    // Read data line by line
    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string cell;

        // Read the first column
        if (std::getline(lineStream, cell, ',')) {
            try {
                // Remove any leading/trailing spaces
                cell.erase(std::remove_if(cell.begin(), cell.end(), ::isspace), cell.end());
                double value1 = std::stod(cell);
                column1.push_back(value1);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid input for column 1: " << cell << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Input out of range for column 1: " << cell << std::endl;
            }
        }

        // Read the second column
        if (std::getline(lineStream, cell, ',')) {
            try {
                // Remove any leading/trailing spaces
                cell.erase(std::remove_if(cell.begin(), cell.end(), ::isspace), cell.end());
                double value2 = std::stod(cell);
                column2.push_back(value2);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid input for column 2: " << cell << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Input out of range for column 2: " << cell << std::endl;
            }
        }

        // Read the third column
        if (std::getline(lineStream, cell, ',')) {
            try {
                // Remove any leading/trailing spaces
                cell.erase(std::remove_if(cell.begin(), cell.end(), ::isspace), cell.end());
                double value3 = std::stod(cell);
                column3.push_back(value3);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid input for column 3: " << cell << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Input out of range for column 3: " << cell << std::endl;
            }
        }
    }

    file.close();

    // Print the data to verify
    std::cout << "Column 1: ";
    for (double val : column1) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Column 2: ";
    for (double val : column2) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Column 3: ";
    for (double val : column3) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}

