#include <iostream>
#include <vector>
#include <fstream> // for file I/O
#include <cmath> // for absolute value
#include <algorithm> // for maximum
#include <chrono> // for timing
using namespace std;
using Matrix = vector<vector<double>>; // define a type alias for matrix of doubles

// Function to initialize a 2D grid with boundary conditions
Matrix initialize_grid(int N, double initial_value) {
    // create an NxN grid filled with the initial value
    Matrix grid(N, vector<double>(N, initial_value));

    // set boundary conditions
    for (int i = 0; i < N; ++i) {
        grid[0][i] = 0.0;       // top boundary
        grid[N-1][i] = 0.0;     // bottom boundary
        grid[i][0] = 0.0;       // left boundary
        grid[i][N-1] = 0.0;     // right boundary
    }

    // insert optional central block
    for (int i = N/3; i < 2*N/3; ++i) {
        for (int j = N/3; j < 2*N/3; ++j) {
            grid[i][j] = 1.0; // central block of high concentration
        }
    }

    return grid;
}

// Function to update the grid values based on advection-diffusion equation
void update_solution(Matrix &grid, const Matrix &velocity_x, const Matrix &velocity_y, double D, double dt, double dx, double dy) {
    int N = grid.size(); // number of grid points
    Matrix new_grid = grid; // create copy of current grid

    // iterate over the interior grid points
    for (int i = 1; i < (N - 1); ++i) {
        for (int j = 1; j < (N - 1); ++j) {

            // compute advection term using upwind differencing
            double advection_x;
            if (velocity_x[i][j] > 0) {
                advection_x = velocity_x[i][j] * (grid[i][j] - grid[i][j-1]) / dx;
            } else {
                advection_x = velocity_x[i][j] * (grid[i][j+1] - grid[i][j]) / dx;
            }

            double advection_y;
            if (velocity_y[i][j] > 0) {
                advection_y = velocity_y[i][j] * (grid[i][j] - grid[i-1][j]) / dy;
            } else {
                advection_y = velocity_y[i][j] * (grid[i+1][j] - grid[i][j]) / dy;
            }

            // compute diffusion term using central difference
            double diffusion = D * ((grid[i+1][j] - 2 * grid[i][j] + grid[i-1][j]) / (dx * dx) +
                                    (grid[i][j+1] - 2 * grid[i][j] + grid[i][j-1]) / (dy * dy));
            
            // store value in new grid
            new_grid[i][j] = max(0.0, grid[i][j] - dt * (advection_x + advection_y) + dt * diffusion);
        }
    }
    grid = new_grid; // update grid with new values
}

// Function to save the grid data to a CSV file
void save_grid_to_csv(const Matrix &grid, const string &filename) {
    ofstream file(filename); // open the file for writing

    // loop through each row of the grid
    for (const auto &row : grid) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j]; // write the value
            if (j < row.size() - 1) file << ","; // add a comma if not the last column
        }
        file << "\n"; // newline after each row
    }

    file.close(); // close the file
}

// Main program for the simulation
int main() {
    auto start = chrono::high_resolution_clock::now(); // start timer

    int N = 600; // grid size
    double dx = 1.0 / N; // grid spacing in x-direction
    double dy = 1.0 / N; // grid spacing in y-direction

    // initialize grid with initial values
    Matrix grid = initialize_grid(N, 0.0);
    save_grid_to_csv(grid, "data/initial.csv");

    double D = 0.01; // diffusion coefficient

    // initialize velocity fields
    Matrix velocity_x(N, vector<double>(N, 2.0)); // velocity in x-direction (right=positive)
    Matrix velocity_y(N, vector<double>(N, 0.0)); // velocity in y-direction (down=positive)

    // run simulation for set time steps
    int num_steps = 100;
    for (int t = 0; t < num_steps; ++t) {
        
        // calculate stable time step
        double max_velocity = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                max_velocity = max(max_velocity, max(fabs(velocity_x[i][j]), fabs(velocity_y[i][j])));
            }
        }

        double dt_advection = dx / max_velocity; // Courant-Friedrichs-Lewy condition
        double dt_diffusion = (dx * dx) / (4 * D); // stability condition for diffusion
        double dt = min(dt_advection, dt_diffusion) / 2; // choose smaller time step
        if (max_velocity == 0 && D == 0) {
            double dt = 0.01; // default time step
        }

        // save grid data every 1/10th of the entire simulation
        if (t % (num_steps/10) == 0) {
            save_grid_to_csv(grid, "data/output_step_" + to_string(t) + ".csv");
        }

        update_solution(grid, velocity_x, velocity_y, D, dt, dx, dy); // update grid values
    }

    auto end = chrono::high_resolution_clock::now(); // end timer
    chrono::duration<double> elapsed = end - start; // calculate runtime
    cout << "Elapsed time: " << elapsed.count() << " seconds" << endl; // output runtime

    return 0; 
}