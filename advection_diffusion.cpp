#include <iostream>
#include <vector>
#include <fstream> // for file I/O
#include <cmath> // for absolute value
#include <algorithm> // for maximum
#include <chrono> // for timing
#include "omp.h" // openmp
using namespace std;
using Matrix = vector<vector<double>>; // define a type alias for matrix of doubles

// Node class for the quadtree
struct Node {
    double value; // the value at this node
    int level;    // refinement level of this node
    Node* children[4]; // array of pointers to child nodes

    // constructure initializing a node with value and refinement
    Node(double val, int lvl) : value(val), level(lvl) {
        for (int i = 0; i < 4; ++i) children[i] = nullptr; // initialize child nodes as nullptr
    }

    // check if this node is a leaf (no subdivisions)
    bool is_leaf() const {
        return children[0] == nullptr;
    }
};

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

    /*
    // insert optional central block
    for (int i = N/3; i < 2*N/3; ++i) {
        for (int j = N/3; j < 2*N/3; ++j) {
            grid[i][j] = 1.0; // central block of high concentration
        }
    }
    */

    int particle_size = N / 10; // Each particle occupies about 1/10th of N
    int num_particles = 5;     // Number of particles (can be adjusted)

    for (int p = 0; p < num_particles; ++p) {
        // Determine particle center coordinates
        int center_i = (p + 1) * (N / (num_particles + 1)); // Spread particles evenly
        int center_j = (p + 1) * (N / (num_particles + 1));

        // Fill a square around the center
        for (int i = max(1, center_i - particle_size / 2); i < min(N - 1, center_i + particle_size / 2); ++i) {
            for (int j = max(1, center_j - particle_size / 2); j < min(N - 1, center_j + particle_size / 2); ++j) {
                grid[i][j] = 1.0; // Particle concentration
            }
        }
    }

    return grid;
}

// Function to calculate the error at a grid cell based on gradients
double calculate_error(const Matrix &grid, int i, int j, double dx, double dy) {
    // calculate partial derivatives using central differencing
    double dudx = (grid[i+1][j] - grid[i-1][j]) / (2 * dx);
    double dudy = (grid[i][j+1] - grid[i][j-1]) / (2 * dy);

    // combine the absolute values of gradients in x and y directions
    return (fabs(dudx) + fabs(dudy));
}

// Function to refine a node by splitting it into 4 child nodes
void refine_grid(Node* node) {
    if (!node->is_leaf()) return; // if the node is already refined, do nothing

    double val = node->value;     // get the value of the current node
    int next_level = node->level + 1; // increment the refinement level

    // create 4 child nodes with the parent's value and next refinement level
    for (int i = 0; i < 4; ++i) {
        node->children[i] = new Node(val, next_level);
    }
}

// Function to coarsen a node by merging its 4 child nodes into the parent
void coarsen_grid(Node* node) {
    if (node->is_leaf()) return; // if the node is already a leaf, do nothing

    double avg_value = 0.0; // initialize the average value

    // loop through the 4 child nodes to calculate their average value
    for (int i = 0; i < 4; ++i) {
        avg_value += node->children[i]->value; // accumulate the value of each child
        delete node->children[i];             // delete the child node to free memory
        node->children[i] = nullptr;          // set the child pointer to null
    }

    // set the parent node's value to the average of its children
    node->value = avg_value / 4.0;
}

// Function to recursively check leaf nodes to find deepest level
int find_deepest_level(Node* node) {
    if (node->is_leaf()) return node->level; // Leaf node: return its level

    int max_level = node->level; // Current level
    for (int i = 0; i < 4; ++i) {
        if (node->children[i] != nullptr) {
            max_level = std::max(max_level, find_deepest_level(node->children[i]));
        }
    }
    return max_level;
}

// Function using deepest level to calculate smallest cell size
double find_min_dx(Node* root, double domain_size, int initial_grid_size) {
    int deepest_level = find_deepest_level(root); // Find the deepest level in the quadtree
    double initial_dx = domain_size / initial_grid_size; // Coarsest grid cell size
    return initial_dx / pow(2, deepest_level); // Adjust for the refinement level
}

// Function to apply AMR logic: refine or coarsen grid cells based on error
void apply_amr(Node* root, const Matrix &grid, double dx, double dy, double refine_threshold, double coarsen_threshold) {
    //#pragma omp parallel for collapse(2)
    for (int i = 1; i < grid.size() - 1; ++i) {     // loop through interior cells (avoid boundaries)
        for (int j = 1; j < grid[0].size() - 1; ++j) {
            // calculate the error at the current grid cell
            double error = calculate_error(grid, i, j, dx, dy);

            //#pragma omp critical // lock to ensure thread-safety
            //{
                // refine the cell if the error exceeds the refine threshold
                if (error > refine_threshold) {
                    refine_grid(root);
                }
                // coarsen the cell if the error is below the coarsen threshold
                else if (error < coarsen_threshold) {
                    coarsen_grid(root);
                }
            //}
        }
    }
}

// Function to update the grid values based on advection-diffusion equation
void update_solution(Matrix &grid, const Matrix &velocity_x, const Matrix &velocity_y, double D, double dt, double dx, double dy) {
    int N = grid.size(); // number of grid points
    Matrix new_grid = grid; // create copy of current grid

    // Parallelize the grid update
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < (N - 1); ++i) {
        for (int j = 1; j < (N - 1); ++j) { // iterate over the interior grid points

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

    omp_set_num_threads(4);

    // Initialize grid and parameters
    int N = 300; // grid size
    double dx = 1.0 / N, dy = 1.0 / N; // grid spacing
    Matrix grid = initialize_grid(N, 0.0); // initialize grid
    Node* root = new Node(0.0, 0); // Create a root node with initial value
    double D = 0.2; // diffusion coefficient
    Matrix velocity_x(N, vector<double>(N, 100.0)); // velocity in x-direction (right=positive)
    Matrix velocity_y(N, vector<double>(N, 0.0)); // velocity in y-direction (down=positive)
    double refine_threshold = (0.20 * 100); // threshold for refinement (20% of max error)
    double coarsen_threshold = (0.05 * 100); // threshold for coursening (5% of max error)

    // Time-stepping loop
    int num_steps = 200;
    for (int t = 0; t < num_steps; ++t) {

        // apply amr to refine/coarsen the grid (every 5 steps)
        if (t % 5 == 0) {
            apply_amr(root, grid, dx, dy, refine_threshold, coarsen_threshold);
        }

        // find smallest cell size (min_dx) after AMR
        double min_dx = find_min_dx(root, 1.0, 4);

        // calculate stable time step
        double max_velocity = 0.0;
        //#pragma omp parallel for reduction(max : max_velocity) collapse(2)
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

        // save grid data every 1/20th of the entire simulation
        if (t % (num_steps/20) == 0) {
            save_grid_to_csv(grid, "data/output_step_" + to_string(t) + ".csv");
        }

        update_solution(grid, velocity_x, velocity_y, D, dt, dx, dy); // update grid values
    }

    auto end = chrono::high_resolution_clock::now(); // end timer
    chrono::duration<double> elapsed = end - start; // calculate runtime
    cout << "Elapsed time: " << elapsed.count() << " seconds" << endl; // output runtime

    return 0; 
}