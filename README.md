# **High-Performance Adaptive Solver for the Advection-Diffusion Equation Using OpenMP and MPI**

---

### **Project Overview:**  
This project aims to develop a high-performance solver for the **Advection-Diffusion Equation** using **Adaptive Mesh Refinement (AMR)** to dynamically adjust grid resolution based on the solution’s behavior. The advection-diffusion equation models the combined effects of transport (advection) and spreading (diffusion), making it both mathematically interesting and computationally intensive. The solver will leverage **MPI** for distributed memory parallelism and **OpenMP** for shared memory parallelism to achieve scalability.

---

### **Advection-Diffusion Equation:**
The governing equation:
\[
\frac{\partial u}{\partial t} + \vec{v} \cdot \nabla u = D \nabla^2 u
\]
Where:
- \(u(x, y, t)\) is the scalar field (e.g., concentration or temperature).
- \(\vec{v} = (v_x, v_y)\) is the velocity vector (advection term).
- \(D\) is the diffusion coefficient.

---

### **Objectives:**
1. **Develop a Baseline Solver:**
   - Implement a 2D finite difference scheme for the advection-diffusion equation on a uniform grid.
2. **Incorporate Adaptive Mesh Refinement (AMR):**
   - Dynamically refine/coarsen the grid based on steep gradients or high error regions.
3. **Parallelize the Solver:**
   - Use **MPI** for domain decomposition and communication between processes.
   - Use **OpenMP** to parallelize computations within each process.
4. **Evaluate Scalability:**
   - Analyze the solver’s performance in terms of runtime, speedup, and load balancing across nodes and threads.
5. **Visualize Results:**
   - Generate visualizations of the solution evolution and the adaptive grid.

---

### **Methodology:**

#### **Step 1: Discretization**
- Use an explicit finite difference method:
  - Advection: Upwind or central differencing.
  - Diffusion: Central differencing.
- Discretized equation:
  \[
  u_{i,j}^{n+1} = u_{i,j}^n - \Delta t \vec{v} \cdot \nabla u + D \Delta t \nabla^2 u
  \]

#### **Step 2: Adaptive Mesh Refinement**
- Implement AMR with a **quadtree grid**:
  - Refine regions with steep gradients or large solution changes.
  - Coarsen regions with low variation.
  - Use a threshold-based error estimator, such as:
    \[
    \text{Error} = \left| \frac{\partial u}{\partial x} \right| + \left| \frac{\partial u}{\partial y} \right|
    \]

#### **Step 3: Parallelization**
1. **MPI (Distributed Memory):**
   - Decompose the computational domain into subdomains for each MPI process.
   - Use communication primitives to synchronize boundary conditions and adaptively refine/coarsen the grid across processes.
2. **OpenMP (Shared Memory):**
   - Parallelize time-stepping, error estimation, and grid refinement within each subdomain.

#### **Step 4: Testing and Scalability**
- Test the solver on increasing grid sizes and number of processes.
- Measure:
  - Runtime vs. grid size.
  - Load balancing across processes.

#### **Step 5: Visualization**
- Generate visualizations of:
  - Solution evolution (e.g., concentration spreading).
  - Adaptive grid changes over time.

---

### **Expected Outcomes:**
1. A fully functional 2D advection-diffusion solver with adaptive mesh refinement.
2. Scalability analysis showing performance improvements with parallelization.
3. Visualizations of the adaptive grid and solution dynamics over time.
4. A compelling demonstration of AMR’s efficiency and accuracy compared to a uniform grid solver.

---

### **Tools and Technologies:**
- **Programming Language:** C++
- **Parallelization Libraries:** OpenMP, MPI
- **Visualization Tools:** Matplotlib, Paraview, or VTK
- **Testing Environment:** HPC cluster or multi-core machine

---

### **Evaluation Criteria:**
1. **Correctness:** Accurate solution compared to analytical benchmarks.
2. **Performance:** Demonstration of runtime improvements with AMR and parallelization.
3. **Scalability:** Efficient use of compute resources across multiple nodes.
4. **Visualization:** Clear and engaging representations of adaptive grids and solution dynamics.
