
# GenPen: Generalized Smoothing for Linear ODEs

This repository provides MATLAB and Fortran code implementing the generalized smoothing methodology for linear ordinary differential equations (ODEs) as introduced in the aforementioned paper. The approach extends traditional smoothing techniques by incorporating ODE constraints, enabling more accurate estimation of dynamic systems from noisy data.

## Overview

The generalized smoother addresses the problem of estimating the state of a system governed by linear ODEs when observations are noisy and possibly incomplete. By integrating the ODE structure into the smoothing process, the method ensures that the estimated trajectories adhere closely to the underlying dynamics.

Key features include:

- Handling both homogeneous and non-homogeneous linear ODEs.
- Incorporation of B-spline basis functions for flexible function approximation.
- Implementation of generalized cross-validation (GCV) for optimal smoothing parameter selection.
- Support for second-order ODE simulations and real-world data applications.

## 📁 Repository Structure

The repository contains the following key files:

- `Gen_Pen.m`: Main function for generalized smoothing with homogeneous ODEs.
- `Gen_Pen_Non.m`: Extension for non-homogeneous ODEs.
- `Generate_sol_ODE.m`: Generates solutions for specified ODEs.
- `Generate_sol_homogenous_ODE.m`: Generates solutions for homogeneous ODEs.
- `Second_order_ODE_simulations.m`: Scripts for simulating second-order ODEs.
- `GCV.F`: Fortran code for computing the GCV criterion.
- `bsplinepen.m`, `bsplinepenJ.m`: Functions for B-spline penalty computations.
- `dbeta_dy.m`, `dbeta_dyNH.m`: Derivatives of beta with respect to y for homogeneous and non-homogeneous cases.
- `dc_dbeta.m`, `dc_dbetaNH.m`: Derivatives of c with respect to beta.
- `dpen.m`: Computes the penalty term.
- `min_gcv.m`, `min_gcv_non.m`: Functions to minimize the GCV criterion for homogeneous and non-homogeneous cases.
- `sintpoly3.F`: Fortran code related to polynomial integration.
- `stinv.m`: Computes the inverse of the matrix.
- `GS_Diabetic_data.m`: Application of the method to diabetic patient data.
- `Melanoma.m`: Application to melanoma dataset.
- `Simulated_Data_Non_Homo.m`: Simulation for non-homogeneous ODE data.

## 🛠️ Requirements

- MATLAB (version R2014b or later recommended)
- Fortran compiler (e.g., gfortran) for compiling `.F` files
- Basic knowledge of ODEs and numerical methods

## 🚀 Getting Started

1. **Clone the repository:**

   ```bash
   git clone https://github.com/mcareyucd/GenPen.git
   cd GenPen
   ```

2. **Compile Fortran code (if necessary):**

   Ensure that the Fortran compiler is available in your system path.

   ```bash
   gfortran -o GCV GCV.F
   gfortran -o sintpoly3 sintpoly3.F
   ```

3. **Run MATLAB scripts:**

   Open MATLAB and navigate to the `GenPen` directory. You can then run the desired scripts, for example:

   ```matlab
   Second_order_ODE_simulations
   ```

   or apply the method to real data:

   ```matlab
   GS_Diabetic_data
   Melanoma
   ```

## 📈 Applications

The generalized smoothing approach is applicable in various domains where systems are modeled by linear ODEs, including:

- Biomedical engineering (e.g., glucose-insulin dynamics)
- Pharmacokinetics/pharmacodynamics modeling
- Environmental modeling
- Engineering systems analysis

## 📖 Citation

If you utilize this code in your research, please cite the original paper:

Carey, M., Gath, E. G., & Hayes, K. (2017). A Generalized Smoother for Linear Ordinary Differential Equations. *Journal of Computational and Graphical Statistics*, 26(3), 671–681. [https://doi.org/10.1080/10618600.2016.1265526](https://doi.org/10.1080/10618600.2016.1265526)

## 📝 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
