#include <cmath>
#include <fstream>
#include <iomanip>
#include "BLAS_Wrapper.h"
#include "Burgers.h"
#include <iostream>
using namespace std;

/**
 * @brief Public Constructor: Accepts a Model instance reference as input
 * @param &m reference to Model instance
 * */
Burgers::Burgers(Model &m) {
    model = &m;

    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    /// Allocate memory to instance variables
    U = new double[Nyr*Nxr];
    V = new double[Nyr*Nxr];
    NextU = new double[Nyr*Nxr];
    NextV = new double[Nyr*Nxr];
}

/**
 * @brief Destructor: Deletes all allocated pointers in the class instance
 * */
Burgers::~Burgers() {
    /// Delete U and V
    delete[] U;
    delete[] V;
    delete[] NextU;
    delete[] NextV;
    /// model is not dynamically alloc
}

/**
 * @brief Sets initial velocity field in x,y for U0 (V0 = U0)
 * */
void Burgers::SetInitialVelocity() {
    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double x0 = model->GetX0();
    double y0 = model->GetY0();
    double dx = model->GetDx();
    double dy = model->GetDy();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    /// Compute U0;
    for (int i = 0; i < Nxr; i++) {
        for (int j = 0; j < Nyr; j++) {
            // Assumes x0 and y0 are identifying top LHS of matrix
            double y = y0 - (j+1)*dy;
            double x = x0 + (i+1)*dx;
            double r = pow(x*x+y*y, 0.5);
            // Store in column-major format
            U[i*Nyr+j] = (r <= 1.0)? 2.0*pow(1.0-r,4.0) * (4.0*r+1.0) : 0.0;
            V[i*Nyr+j] = (r <= 1.0)? 2.0*pow(1.0-r,4.0) * (4.0*r+1.0) : 0.0;
        }
    }
}

/**
 * @brief Sets velocity field in x,y for U, V
 * */
void Burgers::SetIntegratedVelocity() {
    /// Get model parameters
    int Nt = model->GetNt();
    double* temp = nullptr;
    /// Compute U, V for every step k
    for (int k = 0; k < Nt-1; k++) {
        ComputeNextVelocityState();
        temp = NextU;
        NextU = U;
        U = temp;

        temp = NextV;
        NextV = V;
        V = temp;
        cout << "step: " << k << "\n";
    }
}

/**
 * @brief Writes the velocity field for U, V into a file
 * IMPORTANT: Run SetIntegratedVelocity() first
 * */
void Burgers::WriteVelocityFile() {
    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    /// Alloc 2D pointer
    double** Vel = new double*[Nyr];
    for (int j = 0; j < Nyr; j++) {
        Vel[j] = new double[Nxr];
    }

    /// Write U, V into "data.txt"
    ofstream of;
    of.open("data.txt", ios::out | ios::trunc);
    of.precision(4); // 4 s.f.
    /// Write U velocities
    of << "U velocity field:" << endl;
    wrap(U, Nyr, Nxr, Vel);
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            if (j == 0 || i == 0 || j == Ny-1 || i == Nx-1) {
                of << 0 << ' ';
            }
            else {
                of << Vel[j-1][i-1] << ' ';
            }
        }
        of << endl;
    }
    /// Write V velocities
    of << "V velocity field:" << endl;
    wrap(V, Nyr, Nxr, Vel);
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            if (j == 0 || i == 0 || j == Ny-1 || i == Nx-1) {
                of << 0 << ' ';
            }
            else {
                of << Vel[j-1][i-1] << ' ';
            }
        }
        of << endl;
    }
    of.close();

    /// Delete 2D temp pointer
    for (int j = 0; j < Nyr; j++) {
        delete[] Vel[j];
    }
    delete[] Vel;
}

/**
 * @brief Calculates and sets energy of each velocity field per timestamp
 * */
void Burgers::SetEnergy() {
    /// Get Model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double dx = model->GetDx();
    double dy = model->GetDy();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    /// Calculate Energy
    double ddotU = F77NAME(ddot)(Nyr*Nxr, U, 1, U, 1);
    double ddotV = F77NAME(ddot)(Nyr*Nxr, V, 1, V, 1);
    E = 0.5 * (ddotU + ddotV) * dx*dy;
}

/**
 * @brief Computes linear and non-linear terms for U and V
 * */
void Burgers::ComputeNextVelocityState() {
    /// Get model parameters
    int Nyr = model->GetNy() - 2;
    int Nxr = model->GetNx() - 2;

    /// Compute first, second derivatives, & non-linear terms
    double alpha_sum = model->GetAlpha_Sum();
    double beta_dx_sum = model->GetBetaDx_Sum();
    double beta_dy_sum = model->GetBetaDy_Sum();
    double beta_dx_2 = model->GetBetaDx_2();
    double beta_dy_2 = model->GetBetaDy_2();
    double bdx = model->GetBDx();
    double bdy = model->GetBDy();

    /// Pointers to row shifts in U,V
    int iPlus, iMinus;
    double bdxU, bdyV;
    for (int i = 0; i < Nxr; i++) {
        int start = i*Nyr;
        iMinus = (i-1)*Nyr;
        iPlus = (i+1)*Nyr;
        for (int j = 0; j < Nyr; j++) {
            int curr = start + j;
            bdxU = bdx * U[curr];
            bdyV = bdy * V[curr];

            double alpha_total = alpha_sum - bdxU - bdyV;
            NextU[curr] = alpha_total * U[curr];
            NextV[curr] = alpha_total * V[curr];
            if (i < Nxr-1) {
                NextU[curr] += beta_dx_2 * U[iPlus+j];
                NextV[curr] += beta_dx_2 * V[iPlus+j];
            }
            if (i > 0) {
                double bdxU_total = bdxU + beta_dx_sum;
                NextU[curr] += bdxU_total * U[iMinus+j];
                NextV[curr] += bdxU_total * V[iMinus+j];
            }
            if (j < Nyr-1) {
                NextU[curr] += beta_dy_2 * U[curr+1];
                NextV[curr] += beta_dy_2 * V[curr+1];
            }
            if (j > 0) {
                double bdyV_total = bdyV + beta_dy_sum;
                NextU[curr] += bdyV_total * U[curr-1];
                NextV[curr] += bdyV_total * V[curr-1];
            }
        }
    }

    for (int k = 0; k < Nyr*Nxr; k++) {
        NextU[k] += U[k];
        NextV[k] += V[k];
    }
}

/**
 * @brief Wraps a column-major 1D pointer into a pre-allocated row-major 2D pointer
 * @param A 1D pointer in column-major format
 * @param Nxr Nxr
 * @param Nyr Nyr
 * @param res 2D pointer pre-allocated with memory
 * */
void Burgers::wrap(double* A, int Nyr, int Nxr, double** res) {
    for (int i = 0; i < Nyr*Nxr; i++) {
        int col = i / Nyr; // divisor result
        int row = i % Nyr; // remainder after division
        res[row][col] = A[i];
    }
}