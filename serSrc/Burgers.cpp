#include <cmath>
#include <fstream>
#include <iomanip>
#include "BLAS_Wrapper.h"
#include "Helpers.h"
#include "Burgers.h"

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

    /// Matrix coefficients
    dVel_dx_2_coeffs = new double[Nxr*Nxr];
    dVel_dy_2_coeffs = new double[Nyr*Nyr];
    dVel_dx_coeffs = new double[Nxr*Nxr];
    dVel_dy_coeffs = new double[Nyr*Nyr];

    /// Term arrays
    dVel_2 = new double[Nyr*Nxr];
    dVel = new double[Nyr*Nxr];
}

/**
 * @brief Destructor: Deletes all allocated pointers in the class instance
 * */
Burgers::~Burgers() {
    /// Delete matrix coefficients
    delete[] dVel_dx_2_coeffs;
    delete[] dVel_dy_2_coeffs;
    delete[] dVel_dx_coeffs;
    delete[] dVel_dy_coeffs;

    /// Delete U and V
    delete[] U;
    delete[] V;

    /// Delete term arrays
    delete[] dVel_2;
    delete[] dVel;

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
            double r = ComputeR(x, y);
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

    /// Set Matrix Coefficients
    SetMatrixCoefficients();

    /// Compute U, V for every step k
    for (int k = 0; k < Nt-1; k++) {
        double* NextU = NextVelocityState(U, V, true);
        double* NextV = NextVelocityState(U, V, false);
        /// Delete current pointer and point to NextVel
        delete[] U;
        delete[] V;
        U = NextU;
        V = NextV;
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
 * @param Ui U velocity per timestamp
 * @param Vi V velocity per timestamp
 * @param SELECT_U true if the computation is for U
 * */
double* Burgers::NextVelocityState(double* Ui, double* Vi, bool SELECT_U) {
    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double dt = model->GetDt();
    double bdx = model->GetBDx();
    double bdy = model->GetBDy();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    /// Set aliases for computation
    double* Vel = (SELECT_U) ? Ui : Vi;
    double* Other = (SELECT_U)? Vi : Ui;

    /// Generate NextVel
    double* NextVel = new double[Nyr*Nxr];

    /// Compute first & second derivatives
    // x
    for (int i = 0; i < Nyr; i++) {
        F77NAME(dgbmv)('N', Nxr, Nxr, 1, 1, 1.0, dVel_dx_2_coeffs, Nxr, &(Vel[i]), Nyr, 0.0, &(dVel_2[i]), Nyr);
        F77NAME(dgbmv)('N', Nxr, Nxr, 1, 0, 1.0, dVel_dx_coeffs, Nxr, &(Vel[i]), Nyr, 0.0, &(dVel[i]), Nyr);
    }
    // y
    for (int i = 0; i < Nyr*Nxr; i += Nyr) {
        F77NAME(dgbmv)('N', Nyr, Nyr, 1, 1, 1.0, dVel_dy_2_coeffs, Nyr, &(Vel[i]), 1, 1.0, &(dVel_2[i]), 1);
        F77NAME(dgbmv)('N', Nyr, Nyr, 1, 0, 1.0, dVel_dy_coeffs, Nyr, &(Vel[i]), 1, 1.0, &(dVel[i]), 1);
    }

    /// Matrix addition through all terms
    if (SELECT_U) {
        for (int i = 0; i < Nyr*Nxr; i++) {
            double Vel_Vel = bdx * Vel[i] * Vel[i];
            double Vel_Other = bdy * Vel[i] * Other[i];
            double Vel_Vel_Minus_1 = (i < Nyr)? 0 : bdx * Vel[i-Nyr] * Vel[i];
            double Vel_Other_Minus_1 = (i % Nyr == 0)? 0 : bdy * Vel[i-1] * Other[i];

            NextVel[i] = dVel_2[i] - dVel[i] -
                         (Vel_Vel + Vel_Other - Vel_Vel_Minus_1 - Vel_Other_Minus_1);
            NextVel[i] *= dt;
            NextVel[i] += Vel[i];
        }
    }
    else {
        for (int i = 0; i < Nyr*Nxr; i++) {
            double Vel_Vel = bdy * Vel[i] * Vel[i];
            double Vel_Other = bdx * Vel[i] * Other[i];
            double Vel_Vel_Minus_1 = (i % Nyr == 0)? 0 : bdy * Vel[i-1] * Vel[i];
            double Vel_Other_Minus_1 = (i < Nyr)? 0 : bdx * Vel[i-Nyr] * Other[i];

            NextVel[i] = dVel_2[i] - dVel[i] -
                         (Vel_Vel + Vel_Other - Vel_Vel_Minus_1 - Vel_Other_Minus_1);
            NextVel[i] *= dt;
            NextVel[i] += Vel[i];
        }
    }

    return NextVel;
}

/**
 * @brief Private helper function that sets matrix coefficients for differentials
 * */
void Burgers::SetMatrixCoefficients() {
    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    /// Burger constants
    double alpha_dx_2 = model->GetAlphaDx_2();
    double beta_dx_2 = model->GetBetaDx_2();
    double alpha_dy_2 = model->GetAlphaDy_2();
    double beta_dy_2 = model->GetBetaDy_2();
    double alpha_dx_1 = model->GetAlphaDx_1();
    double beta_dx_1 = model->GetBetaDx_1();
    double alpha_dy_1 = model->GetAlphaDy_1();
    double beta_dy_1 = model->GetBetaDy_1();

    /// Set coefficients (alpha along the LD)
    GenSymmBanded(alpha_dx_2, beta_dx_2, Nxr, dVel_dx_2_coeffs);
    GenSymmBanded(alpha_dy_2, beta_dy_2, Nyr, dVel_dy_2_coeffs);
    GenTrmmBanded(alpha_dx_1, beta_dx_1, Nxr, dVel_dx_coeffs);
    GenTrmmBanded(alpha_dy_1, beta_dy_1, Nyr, dVel_dy_coeffs);
}

/**
 * @brief Private helper function computing R for SetInitialVelocity()
 * @param x real x distance to origin
 * @param y real y distance to origin
 * */
double Burgers::ComputeR(double x, double y) {
    double r = pow(x*x+y*y, 0.5);
    return r;
}
