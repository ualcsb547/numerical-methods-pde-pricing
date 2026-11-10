# Numerical Methods, PDEs, and Optimization Engine 📐

A specialized repository dedicated to deterministic algorithms, numerical optimization schemes, and differential equation solvers written entirely in **MATLAB** applied to scientific computing, engineering physics, and financial architectures.

## 🔬 Core Solvers & Algorithmic Blocks

### 1. Partial Differential Equations (PDEs) & Stability Frameworks
*   **Heat Equation Solver (`ejer1_calorMDF.m`):** Finite Difference Methods (FDM) applied to parabolic differential systems modeling heat transfer or asset-volatility diffusion landscapes.
*   **Crank-Nicolson Scheme (`ejer3CN.m`):** Implementation of the implicit, second-order, unconditionally stable Crank-Nicolson method, a mathematical pillar used in quantitative finance to solve complex Black-Scholes pricing grids.
*   **Wave Equation Solver (`ejerondasMDF.m`):** Numerical integration of second-order hyperbolic partial differential equations modeling dynamic wave propagation.

### 2. Finite Difference Methods (MDF) for Boundary Value Problems
*   **Dirichlet Boundaries (`mdf_dirichlet.m`):** Systems architecture engineered to approximate boundary value challenges where edge states are explicitly fixed.
*   **Neumann Boundary Conditions (`mdf_neumann_der.m` / `mdf_neumann_izq.m`):** Algorithms implemented to handle flux-based and derivative boundary values on both left and right domain margins.

### 3. Shooting Methods for Boundary Value Systems (ODEs)
*   **Linear Shooting Engine (`midisplin.m`):** Transforms Boundary Value Problems (BVPs) into independent Initial Value Problems (IVPs) via linear algebra combinations.
*   **Non-Linear Shooting Method (`midisp_nolineal.m`):** Iterative non-linear BVP root-finding solver leveraging mathematical approximations.
*   **Time-Stepping Integrators (`mieuler.m` / `mirk4.m`):** Classical Explicit Euler and High-Order Runge-Kutta 4th Order (RK4) numerical ODE solvers.

---

## 🛠️ Repository Structure
```text
└── matlab/
    ├── 📁 finite_differences/
    │   ├── mdf_dirichlet.m
    │   ├── mdf_neumann_der.m
    │   └── mdf_neumann_izq.m
    ├── 📁 boundary_value_problems/
    │   ├── midisplin.m
    │   ├── midisp_nolineal.m
    │   └── mirk4.m
    └── 📁 pde_solvers/
        ├── ejer1_calorMDF.m
        ├── ejer3CN.m
        └── ejerondasMDF.m
```
