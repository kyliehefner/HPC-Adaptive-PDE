#include <iostream>
#include <vector>
using namespace std;

// Function to initialize a 2D grid with boundary conditions
vector<vector<double>> initialize_grid(int N, double initial_value) {
    // Create an NxN grid filled with the initial value
    vector<vector<double>> grid(N, vector<double>(N, initial_value));

    // Set boundary conditions (Dirichlet with fixed bottom "source")
    for (int i = 0; i < N; ++i) {
        grid[0][i] = 0.0;       // Top boundary: Fixed to 0
        grid[N-1][i] = 1.0;     // Bottom boundary: Fixed to 1
        grid[i][0] = 0.0;       // Left boundary: Fixed to 0
        grid[i][N-1] = 0.0;     // Right boundary: Fixed to 0
    }
    return grid;
}


int main() {
    initialize_grid(10, 0); // 10x10 grid with initial values set to 0
}