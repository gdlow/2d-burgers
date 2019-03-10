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

}

/**
 * @brief Destructor: Deletes all allocated pointers in the class instance
 * */
Burgers::~Burgers() {
    /// Delete U and V
    delete[] U;
    delete[] V;

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

    /// Compute U, V for every step k
    for (int k = 0; k < Nt-1; k++) {
        double* NextU = GetNextU();
        double* NextV = GetNextV();
        /// Delete current pointer and point to NextVel
        delete[] U;
        delete[] V;
        U = NextU;
        V = NextV;
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
 * @brief Private helper function that computes and returns next velocity state based on previous inputs
 * */
double* Burgers::GetNextU() {
    /// Get model parameters
    int Nyr = model->GetNy()-2;
    int Nxr = model->GetNx()-2;
    int NyrNxr = Nyr*Nxr;
    double dt = model->GetDt();
    double bdx = model->GetBDx();
    double bdy = model->GetBDy();

    /// Set aliases for computation
    double* Vel = U;
    double* Other = V;

    /// Generate NextVel
    double* NextVel = new double[NyrNxr];

    /// Compute first, second derivatives, & non-linear terms
    double alpha_sum = model->GetAlpha_Sum();
    double beta_dx_sum = model->GetBetaDx_Sum();
    double beta_dy_sum = model->GetBetaDy_Sum();
    double beta_dx_2 = model->GetBetaDx_2();
    double beta_dy_2 = model->GetBetaDy_2();

    double* Vel_iMinus = nullptr;
    double* Vel_iPlus = nullptr;
    for (int i = 0; i < Nxr; i++) {
        if (i > 0) Vel_iMinus = &(Vel[(i-1)*Nyr]);
        if (i < Nxr-1) Vel_iPlus = &(Vel[(i+1)*Nyr]);
        int start = i*Nyr;
        for (int j = 0; j < Nyr; j++) {
            int curr = start + j;
            NextVel[curr] = (alpha_sum - bdx * Vel[curr] - bdy * Other[curr]) * Vel[curr];
            NextVel[curr] = (i>0)? NextVel[curr] + (bdx * Vel[curr] + beta_dx_sum) * Vel_iMinus[j] : NextVel[curr];
            NextVel[curr] = (j>0)? NextVel[curr] + (bdy * Other[curr] + beta_dy_sum) * Vel[curr-1] : NextVel[curr];
            NextVel[curr] = (i<Nxr-1)? NextVel[curr] + beta_dx_2 * Vel_iPlus[j] : NextVel[curr];
            NextVel[curr] = (j<Nyr-1)? NextVel[curr] + beta_dy_2 * Vel[curr+1] : NextVel[curr];
        }
    }

    /// Loop unrolling to improve cache performance
    int i, j;
    int unrollfactor = 7;
    int maxval = NyrNxr - unrollfactor;
    for (i = 0; i < maxval; i+= unrollfactor+1) {
        NextVel[i] *= dt;
        NextVel[i+1] *= dt;
        NextVel[i+2] *= dt;
        NextVel[i+3] *= dt;
        NextVel[i+4] *= dt;
        NextVel[i+5] *= dt;
        NextVel[i+6] *= dt;
        NextVel[i+7] *= dt;
    }
    for (j = i; j < NyrNxr; j++) {
        NextVel[j] *= dt;
    }
    F77NAME(daxpy)(NyrNxr, 1.0, Vel, 1, NextVel, 1);

    return NextVel;
}

/**
 * @brief Private helper function that computes and returns next velocity state based on previous inputs
 * */
double* Burgers::GetNextV() {
    /// Get model parameters
    int Nyr = model->GetNy()-2;
    int Nxr = model->GetNx()-2;
    int NyrNxr = Nyr*Nxr;
    double dt = model->GetDt();
    double bdx = model->GetBDx();
    double bdy = model->GetBDy();

    /// Set aliases for computation
    double* Vel = V;
    double* Other = U;

    /// Generate NextVel
    double* NextVel = new double[NyrNxr];

    /// Compute first, second derivatives, & non-linear terms
    double alpha_sum = model->GetAlpha_Sum();
    double beta_dx_sum = model->GetBetaDx_Sum();
    double beta_dy_sum = model->GetBetaDy_Sum();
    double beta_dx_2 = model->GetBetaDx_2();
    double beta_dy_2 = model->GetBetaDy_2();

    double* Vel_iMinus = nullptr;
    double* Vel_iPlus = nullptr;
    for (int i = 0; i < Nxr; i++) {
        if (i > 0) Vel_iMinus = &(Vel[(i-1)*Nyr]);
        if (i < Nxr-1) Vel_iPlus = &(Vel[(i+1)*Nyr]);
        int start = i*Nyr;
        for (int j = 0; j < Nyr; j++) {
            int curr = start + j;
            NextVel[curr] = (alpha_sum - bdy * Vel[curr] - bdx * Other[curr]) * Vel[curr];
            NextVel[curr] = (i>0)? NextVel[curr] + (bdx * Other[curr] + beta_dx_sum) * Vel_iMinus[j] : NextVel[curr];
            NextVel[curr] = (j>0)? NextVel[curr] + (bdy * Vel[curr] + beta_dy_sum) * Vel[curr-1] : NextVel[curr];
            NextVel[curr] = (i<Nxr-1)? NextVel[curr] + beta_dx_2 * Vel_iPlus[j] : NextVel[curr];
            NextVel[curr] = (j<Nyr-1)? NextVel[curr] + beta_dy_2 * Vel[curr+1] : NextVel[curr];
        }
    }

    /// Loop unrolling to improve cache performance
    int i, j;
    int unrollfactor = 7;
    int maxval = NyrNxr - unrollfactor;
    for (i = 0; i < maxval; i+= unrollfactor+1) {
        NextVel[i] *= dt;
        NextVel[i+1] *= dt;
        NextVel[i+2] *= dt;
        NextVel[i+3] *= dt;
        NextVel[i+4] *= dt;
        NextVel[i+5] *= dt;
        NextVel[i+6] *= dt;
        NextVel[i+7] *= dt;
    }
    for (j = i; j < NyrNxr; j++) {
        NextVel[j] *= dt;
    }
    F77NAME(daxpy)(NyrNxr, 1.0, Vel, 1, NextVel, 1);

    return NextVel;
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