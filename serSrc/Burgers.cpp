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
    delete[] U; delete[] V;

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

    /// Compute U0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // Assumes x0 and y0 are identifying top LHS of matrix
            double y = y0 - j*dy;
            double x = x0 + i*dx;
            double r = ComputeR(x, y);
            // Store in column-major format
            U[(i-1)*Nyr+(j-1)] = (r <= 1.0)? 2.0*pow(1.0-r,4.0) * (4.0*r+1.0) : 0.0;
            V[(i-1)*Nyr+(j-1)] = (r <= 1.0)? 2.0*pow(1.0-r,4.0) * (4.0*r+1.0) : 0.0;
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
        CopyAndDelete(NextU, NextV);
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
 * @param Ui U velocity per timestamp (i.e. supply U[k])
 * @param Vi V velocity per timestamp (i.e. supply V[k])
 * @param SELECT_U true if the computation is for U
 * */
double* Burgers::NextVelocityState(double* Ui, double* Vi, bool SELECT_U) {
    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double dt = model->GetDt();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double b = model->GetB();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    /// Set aliases for computation
    double* Vel = (SELECT_U) ? Ui : Vi;
    double* Other = (SELECT_U)? Vi : Ui;

    /// Generate term arrays
    double* NextVel = new double[Nyr*Nxr];
    double* dVel_dx_2 = new double[Nyr*Nxr];
    double* dVel_dy_2 = new double[Nyr*Nxr];
    double* dVel_dx = new double[Nyr*Nxr];
    double* dVel_dy = new double[Nyr*Nxr];
    double* Vel_Vel = nullptr;
    double* Vel_Other = nullptr;
    double* Vel_Vel_Minus_1 = nullptr;
    double* Vel_Other_Minus_1 = nullptr;

    /// Compute second derivatives
    F77NAME(dsymm)('R', 'U', Nyr, Nxr, 1.0, dVel_dx_2_coeffs, Nxr, Vel, Nyr, 0.0, dVel_dx_2, Nyr);
    F77NAME(dsymm)('L', 'U', Nyr, Nxr, 1.0, dVel_dy_2_coeffs, Nyr, Vel, Nyr, 0.0, dVel_dy_2, Nyr);

    /// Compute first derivatives
    F77NAME(dcopy)(Nyr*Nxr, Vel, 1, dVel_dx, 1);
    F77NAME(dcopy)(Nyr*Nxr, Vel, 1, dVel_dy, 1);
    F77NAME(dtrmm)('R', 'U', 'N', 'N', Nyr, Nxr, -1.0, dVel_dx_coeffs, Nxr, dVel_dx, Nyr);
    F77NAME(dtrmm)('L', 'L', 'N', 'N', Nyr, Nxr, -1.0, dVel_dy_coeffs, Nyr, dVel_dy, Nyr);

    /// Compute b terms
    if (SELECT_U) {
        Vel_Vel = MatMul(Vel, Vel, Nyr, Nxr, false, false, b/dx);
        Vel_Other = MatMul(Vel, Other, Nyr, Nxr, false, false, b/dy);
        Vel_Vel_Minus_1 = MatMul(Vel, Vel, Nyr, Nxr, true, false, b/dx);
        Vel_Other_Minus_1 = MatMul(Vel, Other, Nyr, Nxr, false, true, b/dy);
    }
    else {
        Vel_Vel = MatMul(Vel, Vel, Nyr, Nxr, false, false, b/dy);
        Vel_Other = MatMul(Vel, Other, Nyr, Nxr, false, false, b/dx);
        Vel_Vel_Minus_1 = MatMul(Vel, Vel, Nyr, Nxr, false, true, b/dy);
        Vel_Other_Minus_1 = MatMul(Vel, Other, Nyr, Nxr, true, false, b/dx);
    }

    /// Matrix addition through all terms
    for (int i = 0; i < Nyr*Nxr; i++) {
        NextVel[i] = dVel_dx_2[i] + dVel_dy_2[i] + dVel_dx[i] + dVel_dy[i] -
                     (Vel_Vel[i] + Vel_Other[i] - Vel_Vel_Minus_1[i] - Vel_Other_Minus_1[i]);
        NextVel[i] *= dt;
        NextVel[i] += Vel[i];
    }

    /// Delete pointers
    delete[] dVel_dx_2;
    delete[] dVel_dy_2;
    delete[] dVel_dx;
    delete[] dVel_dy;
    delete[] Vel_Vel;
    delete[] Vel_Other;
    delete[] Vel_Vel_Minus_1;
    delete[] Vel_Other_Minus_1;

    return NextVel;
}

/**
 * @brief Copies next state velocity into previous state, and deletes the temporary pointer
 * */
void Burgers::CopyAndDelete(double* NextU, double* NextV) {
    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    for (int i = 0; i < Nyr*Nxr; i++) {
        U[i] = NextU[i];
        V[i] = NextV[i];
    }

    delete[] NextU;
    delete[] NextV;
}

/**
 * @brief Private helper function that sets matrix coefficients for differentials
 * */
void Burgers::SetMatrixCoefficients() {
    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double ax = model->GetAx();
    double ay = model->GetAy();
    double c = model->GetC();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    /// Generate and set coefficients
    GenSymm((-2.0*c)/pow(dx,2.0), c/pow(dx,2.0), Nxr, Nxr, dVel_dx_2_coeffs);
    GenSymm((-2.0*c)/pow(dy,2.0), c/pow(dy,2.0), Nyr, Nyr, dVel_dy_2_coeffs);
    GenTrmm(ax/dx, -ax/dx, Nxr, Nxr, true, dVel_dx_coeffs);
    GenTrmm(ay/dy, -ay/dy, Nyr, Nyr, false, dVel_dy_coeffs);
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
