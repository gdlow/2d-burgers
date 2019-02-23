#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include "BLAS_Wrapper.h"
#include "Helpers.h"
#include "Burgers.h"
#include <iostream>

using namespace std;

bool SELECT_U = true;
bool SELECT_V = false;

/**
 * Constructor: Accepts a Model instance pointer as input
 * Sets it as an instance variable
 * */
Burgers::Burgers(Model &m) {
    model = &m;
}

/**
 * Destructor: Deletes all allocated pointers in the class instance
 * */
Burgers::~Burgers() {
    // Get model parameters
    int Nt = model->GetNt();

    // Delete matrix coefficients
    delete[] dVel_dx_2_coeffs;
    delete[] dVel_dy_2_coeffs;
    delete[] dVel_dx_coeffs;
    delete[] dVel_dy_coeffs;

    // Delete Cache
    delete[] Vel_c;
    delete[] Other_c;

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
    // MPI Parameters
    int loc_rank = model->GetRank();

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

    // Split x domain into 2
    int loc_Nxr = (Nxr % 2 != 0 && loc_rank == 0) ? (Nxr/2)+1 : Nxr/2;
    double* loc_U0 = new double[Nyr*loc_Nxr];
    int* recvcounts = new int[2];
    int* displs = new int[2];
    recvcounts[0] = Nyr*((Nxr/2)+1); recvcounts[1] = Nyr*(Nxr/2);
    displs[0] = 0; displs[1] = recvcounts[0];

    // Compute loc_U0 for half domain
    for (int i = 0; i < loc_Nxr; i++) {
        for (int j = 0; j < Nyr; j++) {
            double y = y0 - (j+1)*dy;
            // offset = 1 for proc 0; offset = 2 for proc 1 (1+1)
            double x = (loc_rank == 0)? x0 + (i+1)*dx : x0 + (loc_rank*loc_Nxr+i+2)*dx;
            double r = ComputeR(x, y);
            loc_U0[i*Nyr+j] = (r <= 1.0)? pow(2.0*(1.0-r),4.0) * (4.0*r+1.0) : 0.0;
        }
    }

    // Gather result
    U0 = nullptr;
    U0 = new double[Nyr*Nxr];

    MPI_Allgatherv(loc_U0, Nyr*loc_Nxr, MPI_DOUBLE, U0, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    delete[] loc_U0;
    delete[] recvcounts;
    delete[] displs;
}

/**
 * Sets velocity field in x,y for U, V
 * */
void Burgers::SetIntegratedVelocity() {
    // Get model parameters
    int Nt = model->GetNt();
    int Ny = model->GetNy();

    // Reduced parameters
    int Nyr = Ny - 2;

    // Generate U, V
    U = nullptr; V = nullptr;
    U = new double*[Nt]; V = new double*[Nt];

    // Set initial velocity field
    U[0] = U0; V[0] = U0;

    // Set Matrix Coefficients
    SetMatrixCoefficients();

    // Allocate cache memory
    Vel_c = new double[Nyr];
    Other_c = new double[Nyr];

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
    // MPI Parameters
    int loc_rank = model->GetRank();
    if (loc_rank == 0) {
        // Get model parameters
        int Nt = model->GetNt();
        int Ny = model->GetNy();
        int Nx = model->GetNx();
        double dt = model->GetDt();

        // Reduced parameters
        int Nyr = Ny - 2;
        int Nxr = Nx - 2;

        // Alloc 2D pointer
        double** Vel = new double*[Nyr];
        for (int j = 0; j < Nyr; j++) {
            Vel[j] = new double[Nxr];
        }

        // Write U, V into "data.txt"
        ofstream of;
        of.open("data.txt", ios::out | ios::trunc);
        of.precision(4); // 4 s.f.
        // Write U velocities
        of << "U velocity field:" << endl;
        for (int k = 0; k < Nt; k++) {
            of << "t = " << k*dt << ":" << endl;
            wrap(U[k], Nyr, Nxr, Vel);
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
        }
        // Write V velocities
        of << "V velocity field:" << endl;
        for (int k = 0; k < Nt; k++) {
            of << "t = " << k*dt << ":" << endl;
            wrap(V[k], Nyr, Nxr, Vel);
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
        }
        of.close();

        // Delete 2D temp pointer
        for (int j = 0; j < Nyr; j++) {
            delete[] Vel[j];
        }
        delete[] Vel;
    }
}

void Burgers::SetEnergy() {
    // TODO: Consider Alltoall? (TO TEST)
    // MPI Parameters
    int p = model->GetP();
    int loc_rank = model->GetRank();

    // Get Model parameters
    int Nt = model->GetNt();
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    // Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    // Split x domain into 2
    int loc_Nxr = (Nxr % 2 != 0 && loc_rank == 0) ? (Nxr/2)+1 : Nxr/2;
    double** loc_U = new double*[Nt];
    double** loc_V = new double*[Nt];
    double** U_dots = new double*[Nt];
    double** V_dots = new double*[Nt];

    // Calculate Energy
    E = nullptr;
    E = new double[Nt];
    for (int k = 0; k < Nt; k++) {
        loc_U[k] = new double[loc_Nxr*Nyr];
        loc_V[k] = new double[loc_Nxr*Nyr];
        U_dots[k] = new double[p];
        V_dots[k] = new double[p];
        MPI_Scatter(U[k], loc_Nxr*Nyr, MPI_DOUBLE, loc_U[k], loc_Nxr*Nyr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(V[k], loc_Nxr*Nyr, MPI_DOUBLE, loc_V[k], loc_Nxr*Nyr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double loc_ddotU = F77NAME(ddot)(Nyr*loc_Nxr, loc_U[k], 1, loc_U[k], 1);
        double loc_ddotV = F77NAME(ddot)(Nyr*loc_Nxr, loc_V[k], 1, loc_V[k], 1);
        MPI_Gather(&loc_ddotU, 1, MPI_DOUBLE, U_dots[k], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&loc_ddotV, 1, MPI_DOUBLE, V_dots[k], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (loc_rank == 0) {
            // Consolidate dot prods
            E[k] = 0.5 * (U_dots[k][0] + U_dots[k][1] + V_dots[k][0] + V_dots[k][1]);
        }
    }
    for (int k = 0; k < Nt; k++) {
        delete[] loc_U[k];
        delete[] loc_V[k];
        delete[] U_dots[k];
        delete[] V_dots[k];
    }
    delete[] loc_U;
    delete[] loc_V;
    delete[] U_dots;
    delete[] V_dots;
}

/**
 * Computes and returns next velocity state based on previous inputs
 * */
double* Burgers::NextVelocityState(double* Ui, double* Vi, bool U_OR_V) {
    // MPI Parameters
    int loc_rank = model->GetRank();
    int* sendcounts = new int[2];
    int* displs = new int[2];

    // Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double dt = model->GetDt();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double b = model->GetB();
    double ax = model->GetAx();
    double c = model->GetC();

    // Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    // Split x domain into 2
    int loc_Nxr = (Nxr % 2 != 0 && loc_rank == 0) ? (Nxr/2)+1 : Nxr/2;

    sendcounts[0] = Nyr*((Nxr/2)+1); sendcounts[1] = Nyr*(Nxr/2);
    displs[0] = 0; displs[1] = sendcounts[0];

    // Set aliases for computation
    double* Vel = (U_OR_V) ? Ui : Vi;
    double* Other = (U_OR_V)? Vi : Ui;

    // Local MPI Vel arrays
    double* loc_Vel = new double[Nyr*loc_Nxr];
    double* loc_Other = new double[Nyr*loc_Nxr];
    double* loc_NextVel = new double[Nyr*loc_Nxr];

    // Set cache for Vel and Other
    SetCache(Vel, Vel_c);
    SetCache(Other, Other_c);

    MPI_Scatterv(Vel, sendcounts, displs, MPI_DOUBLE, loc_Vel, sendcounts[loc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(Other, sendcounts, displs, MPI_DOUBLE, loc_Other, sendcounts[loc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Generate term arrays
    double* NextVel = new double[Nyr*Nxr];
    double* dVel_dx_2 = new double[Nyr*loc_Nxr];
    double* dVel_dy_2 = new double[Nyr*loc_Nxr];
    double* dVel_dx = new double[Nyr*loc_Nxr];
    double* dVel_dy = new double[Nyr*loc_Nxr];
    double* Vel_Vel = nullptr;
    double* Vel_Other = nullptr;
    double* Vel_Vel_Minus_1 = nullptr;
    double* Vel_Other_Minus_1 = nullptr;

    // Compute second derivatives
    F77NAME(dsymm)('R', 'U', Nyr, loc_Nxr, 1.0, dVel_dx_2_coeffs, loc_Nxr, loc_Vel, Nyr, 0.0, dVel_dx_2, Nyr);
    F77NAME(dsymm)('L', 'U', Nyr, loc_Nxr, 1.0, dVel_dy_2_coeffs, Nyr, loc_Vel, Nyr, 0.0, dVel_dy_2, Nyr);

    // Compute first derivatives
    F77NAME(dcopy)(Nyr*loc_Nxr, loc_Vel, 1, dVel_dx, 1);
    F77NAME(dcopy)(Nyr*loc_Nxr, loc_Vel, 1, dVel_dy, 1);
    F77NAME(dtrmm)('R', 'U', 'N', 'N', Nyr, loc_Nxr, -1.0, dVel_dx_coeffs, loc_Nxr, dVel_dx, Nyr);
    F77NAME(dtrmm)('L', 'L', 'N', 'N', Nyr, loc_Nxr, -1.0, dVel_dy_coeffs, Nyr, dVel_dy, Nyr);

    // Modify arrays based on cache values
    for (int j = 0; j < Nyr; j++) {
        if (loc_rank==0) {
            // LHS: Fix last col
            dVel_dx_2[(loc_Nxr-1)*Nyr+j] += c/pow(dx,2.0)*Vel_c[j];
        }
        else {
            // RHS: Fix first col
            dVel_dx_2[j] += c/pow(dx,2.0)*Vel_c[j];
            dVel_dx[j] -= ax/dx*Vel_c[j];
        }
    }

    // Compute b terms
    if (U_OR_V == SELECT_U) {
        Vel_Vel = MatMul(loc_Vel, loc_Vel, Nyr, loc_Nxr, false, false, b/dx);
        Vel_Other = MatMul(loc_Vel, loc_Other, Nyr, loc_Nxr, false, false, b/dy);
        Vel_Vel_Minus_1 = MatMul(loc_Vel, loc_Vel, Nyr, loc_Nxr, true, false, b/dx);
        Vel_Other_Minus_1 = MatMul(loc_Vel, loc_Other, Nyr, loc_Nxr, false, true, b/dy);

        // Modify Vel_Vel_Minus_1 to include previous cache
        // Fix first col
        if (loc_rank == 1) {
            for (int i = 0; i < Nyr; i++) {
                Vel_Vel_Minus_1[i] = Vel_c[i] * loc_Vel[i];
            }
        }
    }
    else {
        Vel_Vel = MatMul(loc_Vel, loc_Vel, Nyr, loc_Nxr, false, false, b/dy);
        Vel_Other = MatMul(loc_Vel, loc_Other, Nyr, loc_Nxr, false, false, b/dx);
        Vel_Vel_Minus_1 = MatMul(loc_Vel, loc_Vel, Nyr, loc_Nxr, false, true, b/dy);
        Vel_Other_Minus_1 = MatMul(loc_Vel, loc_Other, Nyr, loc_Nxr, true, false, b/dx);

        // Modify Vel_Other_Minus_1 to include previous cache
        // Fix first col
        if (loc_rank == 1) {
            for (int i = 0; i < Nyr; i++) {
                Vel_Other_Minus_1[i] = Vel_c[i] * loc_Other[i];
            }
        }
    }

    // Matrix addition through all terms
    for (int i = 0; i < Nyr*loc_Nxr; i++) {
        loc_NextVel[i] = dVel_dx_2[i] + dVel_dy_2[i] + dVel_dx[i] + dVel_dy[i] -
                (Vel_Vel[i] + Vel_Other[i] - Vel_Vel_Minus_1[i] - Vel_Other_Minus_1[i]);
        loc_NextVel[i] *= dt;
        loc_NextVel[i] += Vel[i];
    }

    MPI_Allgatherv(loc_NextVel, Nyr*loc_Nxr, MPI_DOUBLE, NextVel, sendcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    // Delete pointers
    delete[] sendcounts;
    delete[] displs;
    delete[] loc_Vel;
    delete[] loc_Other;
    delete[] loc_NextVel;
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
 * Sets matrix coefficients for differentials
 * */
void Burgers::SetMatrixCoefficients() {
    // MPI Parameters
    int loc_rank = model->GetRank();

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

    // Split x domain into 2
    int loc_Nxr = (Nxr % 2 != 0 && loc_rank == 0) ? (Nxr/2)+1 : Nxr/2;

    // Generate and set coefficients
    dVel_dx_2_coeffs = GenSymm((-2.0*c)/pow(dx,2.0), c/pow(dx,2.0), loc_Nxr, loc_Nxr);
    dVel_dy_2_coeffs = GenSymm((-2.0*c)/pow(dy,2.0), c/pow(dy,2.0), Nyr, Nyr);
    dVel_dx_coeffs = GenTrmm(ax/dx, -ax/dx, loc_Nxr, loc_Nxr, true);
    dVel_dy_coeffs = GenTrmm(ay/dy, -ay/dy, Nyr, Nyr, false);
}

double Burgers::ComputeR(double x, double y) {
    double r = pow(x*x+y*y, 0.5);
    return r;
}

void Burgers::SetCache(double* Vel, double* Cache) {
    // Cache should be allocated Nyr memory before this is run

    // MPI Parameters
    int loc_rank = model->GetRank();

    // Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    // Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    // Split x domain into 2
    int loc_Nxr = (Nxr % 2 != 0 && loc_rank == 0) ? (Nxr/2)+1 : Nxr/2;

    // Set Cache
    for (int j = 0; j < Nyr; j++) {
        // Picks out first col of RHS for LHS (col = 5)
        // Picks out last col of LHS for RHS (col = 4)
        Cache[j] = Vel[loc_Nxr*Nyr+j];
    }

}
