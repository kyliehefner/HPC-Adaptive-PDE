#include <iostream>
#include <vector>
using namespace std;
using Matrix = vector<vector<double>>; // define a type alias for matrix of doubles

// Function to initialize a 2D grid with boundary conditions
Matrix initialize_grid(int N, double initial_value) {
    // create an NxN grid filled with the initial value
    Matrix grid(N, vector<double>(N, initial_value));

    // set boundary conditions (Dirichlet with fixed bottom "source")
    for (int i = 0; i < N; ++i) {
        grid[0][i] = 0.0;       // top boundary: fixed to 0
        grid[N-1][i] = 1.0;     // bottom boundary: fixed to 1
        grid[i][0] = 0.0;       // left boundary: fixed to 0
        grid[i][N-1] = 0.0;     // right boundary: fixed to 0
    }
    return grid;
}

// Function to update the grid values based on advection-diffusion equation
void update_solution(Matrix &grid, const Matrix &velocity_x, const Matrix &velocity_y, double D, double dt, double dx, double dy) {
    int N = grid.size(); // number of grid points
    Matrix new_grid = grid; // create copy of current grid

    // iterature over the interior grid points
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            // compute advection term in x and y directions
            double advection_x = velocity_x[i][j] * (grid[i][j] - grid[i-1][j]) / dx;
            double advection_y = velocity_y[i][j] * (grid[i][j] - grid[i][j-1]) / dy;

            // compute diffusion term using central difference
            double diffusion = D * ((grid[i+1][j] - 2 * grid[i][j] + grid[i-1][j]) / (dx * dx) +
                                    (grid[i][j+1] - 2 * grid[i][j] + grid[i][j-1]) / (dy * dy));
            
            // store value in new grid
            new_grid[i][j] = grid[i][j] - dt * (advection_x + advection_y) + dt * diffusion;
        }
    }
    grid = new_grid; // update grid with new values
}

// Main program for the simulation
int main() {
    int N = 10;             // grid size
    double D = 0.1;          // diffusion coefficient
    double dt = 0.01;        // time step
    double dx = 1.0 / N;     // grid spacing in x-direction
    double dy = 1.0 / N;     // grid spacing in y-direction

    // initialize grid with initial values
    Matrix grid = initialize_grid(N, 0.0);

    // initialize velocity fields
    Matrix velocity_x(N, vector<double>(N, 1.0)); // constant velocity in x-direction
    Matrix velocity_y(N, vector<double>(N, 0.0)); // no velocity in y-direction

    // run simulation for 100 time steps
    for (int t = 0; t < 100; ++t) {
        update_solution(grid, velocity_x, velocity_y, D, dt, dx, dy); // update grid values
    }

    // print final grid to console (for debugging)
    for (const auto &row : grid) {
        for (const auto &val : row) {
            cout << val << " ";
        }
        cout << "\n";
    }
    return 0; 
}