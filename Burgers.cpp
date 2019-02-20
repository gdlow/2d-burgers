#include <cmath>
#include <fstream>
#include <iomanip>
#include "BLAS_Wrapper.h"
#include "Helpers.h"
#include "Burgers.h"

using namespace std;

bool SELECT_U = true;
bool SELECT_V = false;

Burgers::Burgers(Model &m) {
    model = &m;
}
Burgers::~Burgers() {
    // Get model parameters
    int Nt = model->GetNt();

    // Delete U and V
    for (int k = 1; k < Nt; k++) {
        // U[0] = V[0] = U0 (not dynamically alloc)
        delete[] U[k];
        delete[] V[k];
    }
    delete[] U; delete[] V;
    U = nullptr; V = nullptr;

    // Delete U0
    delete[] U0;
    U0 = nullptr;

    // Delete E
    delete[] E;
    E = nullptr;

    // model is not dynamically alloc
}

/**
 * Sets initial velocity field in x,y for U0 (V0 = U0)
 * */
void Burgers::SetInitialVelocity() {
    // Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double x0 = model->GetX0();
    double y0 = model->GetY0();
    double dx = model->GetDx();
    double dy = model->GetDy();

    // Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    // Compute U0;
    U0 = nullptr;
    U0 = new double[Nyr*Nxr];
    for (int i = 0; i < Nxr; i++) {
        for (int j = 0; j < Nyr; j++) {
            // Assumes x0 and y0 are identifying top LHS of matrix
            double y = y0 - j*dy;
            double x = x0 + i*dx;
            double r = ComputeR(x, y);
            // Store in column-major format
            U0[i*Nyr+j] = (r <= 1.0)? pow(2.0*(1.0-r),4.0) * (4.0*r-1.0) : 0.0;
        }
    }
}

/**
 * Sets velocity field in x,y for U, V
 * */
void Burgers::SetIntegratedVelocity() {
    // Get model parameters
    int Nt = model->GetNt();

    // Generate U, V
    U = nullptr; V = nullptr;
    U = new double*[Nt]; V = new double*[Nt];

    // Set initial velocity field
    U[0] = U0; V[0] = U0;

    // Set Matrix Coefficients
    SetMatrixCoefficients();

    // Compute U, V for every step k
    for (int k = 0; k < Nt-1; k++) {
        U[k+1] = NextVelocityState(U[k], V[k], SELECT_U);
        V[k+1] = NextVelocityState(U[k], V[k], SELECT_V);
    }
}

/**
 * Writes the velocity field for U, V into a file
 * IMPORTANT: Run SetIntegratedVelocity() first
 * */
void Burgers::WriteVelocityFile() {
    // Get model parameters
    int Nt = model->GetNt();
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    int dt = model->GetDt();

    // Write U, V into "data.txt"
    ofstream of;
    of.open("data.txt", ios::out | ios::trunc);
    of.precision(4); // 4 s.f.
    // Write U velocities
    of << "U velocity field:" << endl;
    for (int k = 0; k < Nt; k++) {
        of << "t = " << k*dt << ":" << endl;
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                of << U[k][j][i] << ' ';
            }
            of << endl;
        }
    }
    // Write V velocities
    of << "V velocity field:" << endl;
    for (int k = 0; k < Nt; k++) {
        of << "t = " << k*dt << ":" << endl;
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                of << V[k][j][i] << ' ';
            }
            of << endl;
        }
    }
    of.close();
}

void Burgers::SetEnergy() {
    // Get Model parameters
    int Nt = model->GetNt();
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    // Calculate Energy
    E = nullptr;
    E = new double[Nt];
    for (int k = 0; k < Nt; k++) {
        double energy = 0;
        // Sum Energy Over Domain
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                double USquare = pow(U[k][j][i], 2.0) + pow(V[k][j][i], 2.0);
                energy += USquare;
            }
        }
        // Prefactor by 1/2
        energy *= 0.5;
        E[k] = energy;
    }
}

double* Burgers::NextVelocityState(double* Ui, double* Vi, bool U_OR_V) {
    // Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double dt = model->GetDt();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double b = model->GetB();

    // Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    // Set aliases for computation
    double* Vel = (U_OR_V) ? Ui : Vi;
    double* Other = (U_OR_V)? Vi : Ui;

    // Generate term arrays
    double* NextVel = new double[Nyr*Nxr];
    double* dVel_dx_2 = new double[Nyr*Nxr];
    double* dVel_dy_2 = new double[Nyr*Nxr];
    double* dVel_dx = new double[Nyr*Nxr];
    double* dVel_dy = new double[Nyr*Nxr];
    double* Vel_Vel = nullptr;
    double* Vel_Other = nullptr;
    double* Vel_Vel_Minus_1 = nullptr;
    double* Vel_Other_Minus_1 = nullptr;

    // Compute second derivatives
    F77NAME(dsymm)('R', 'U', Nyr, Nxr, 1.0, dVel_dx_2_coeffs, Nxr, Vel, Nxr, 0.0, dVel_dx_2, Nxr);
    F77NAME(dsymm)('L', 'U', Nyr, Nxr, 1.0, dVel_dy_2_coeffs, Nyr, Vel, Nxr, 0.0, dVel_dy_2, Nxr);

    // Compute first derivatives
    F77NAME(dcopy)(Nyr*Nxr, Vel, 1, dVel_dx, 1);
    F77NAME(dcopy)(Nyr*Nxr, Vel, 1, dVel_dy, 1);
    F77NAME(dtrmm)('R', 'U', 'N', 'N', Nyr, Nxr, -1.0, dVel_dx_coeffs, Nxr, dVel_dx, Nxr);
    F77NAME(dtrmm)('L', 'L', 'N', 'N', Nyr, Nxr, -1.0, dVel_dy_coeffs, Nyr, dVel_dy, Nxr);

    // Compute b terms
    if (U_OR_V == SELECT_U) {
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

    // Matrix addition through all terms
    for (int i = 0; i < Nyr*Nxr; i++) {
        NextVel[i] = dVel_dx_2[i] + dVel_dy_2[i] + dVel_dx[i] + dVel_dy[i] -
                (Vel_Vel[i] + Vel_Other[i] - Vel_Vel_Minus_1[i] - Vel_Other_Minus_1[i]);
        NextVel[i] *= dt;
        NextVel[i] += Vel[i];
    }
    // Delete pointers
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

void Burgers::SetMatrixCoefficients() {
    // Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double ax = model->GetAx();
    double ay = model->GetAy();
    double c = model->GetC();

    // Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    // Generate and set coefficients
    dVel_dx_2_coeffs = GenSymm((-2.0*c)/pow(dx,2.0), c/pow(dx,2.0), Nxr*Nxr);
    dVel_dy_2_coeffs = GenSymm((-2.0*c)/pow(dy,2.0), c/pow(dy,2.0), Nyr*Nyr);
    dVel_dx_coeffs = GenTrmm(ax/dx, -ax/dx, Nxr*Nxr, true);
    dVel_dy_coeffs = GenTrmm(ay/dy, -ay/dy, Nyr*Nyr, false);
}

double Burgers::ComputeR(double x, double y) {
    double r = pow(x*x+y*y, 0.5);
    return r;
}
