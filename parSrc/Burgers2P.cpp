#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include "BLAS_Wrapper.h"
#include "Helpers.h"
#include "Burgers2P.h"

using namespace std;


/**
 * Constructor: Accepts a Model instance pointer as input
 * Sets it as an instance variable
 * Allocates memory to all other instance variables
 * */
Burgers2P::Burgers2P(Model &m) {
    // Set model class pointer as instance variable
    model = &m;

    // Get model parameters
    int Nt = model->GetNt();
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();

    /* Allocate memory to instance variables */

    // U0
    U0 = new double[Nyr*Nxr];

    // Matrix coefficients
    dVel_dx_2_coeffs = new double[Nxr*Nxr];
    dVel_dy_2_coeffs = new double[Nyr*Nyr];
    dVel_dx_coeffs = new double[Nxr*Nxr];
    dVel_dy_coeffs = new double[Nyr*Nyr];

    // Caches
    upVel = new double[Nxr];
    downVel = new double[Nxr];
    leftVel = new double[Nyr];
    rightVel = new double[Nyr];

    // U, V
    U = new double*[Nt];
    V = new double*[Nt];
    // U[0] = V[0] = U0
    for (int k = 1; k < Nt; k++) {
        U[k] = new double[Nyr*Nxr];
    }

    // E
    E = new double[Nt];
}

/**
 * Destructor: Deletes all allocated pointers in the class instance
 * */
Burgers2P::~Burgers2P() {
    // Get model parameters
    int Nt = model->GetNt();

    // Delete E
    delete[] E;

    // Delete U and V
    for (int k = 1; k < Nt; k++) {
        // U[0] = V[0] = U0 (not dynamically alloc)
        delete[] U[k];
        delete[] V[k];
    }
    delete[] U;
    delete[] V;

    // Delete Caches
    delete[] upVel;
    delete[] downVel;
    delete[] leftVel;
    delete[] rightVel;

    // Delete matrix coefficients
    delete[] dVel_dx_2_coeffs;
    delete[] dVel_dy_2_coeffs;
    delete[] dVel_dx_coeffs;
    delete[] dVel_dy_coeffs;

    // Delete U0
    delete[] U0;

    // model is not dynamically alloc
}

/**
 * Sets initial velocity field in x,y for U0 (V0 = U0)
 * */
void Burgers2P::SetInitialVelocity() {
    // Get model parameters
    double x0 = model->GetX0();
    double y0 = model->GetY0();
    double dx = model->GetDx();
    double dy = model->GetDy();
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    int displ_x = model->GetDisplX();
    int displ_y = model->GetDisplY();

    // Compute U0
    // Memory layout in column-major format
    double loc_x0 = x0 + (displ_x+1)*dx;
    double loc_y0 = y0 - (displ_y+1)*dy;
    for (int i = 0; i < Nxr; i++) {
        for (int j = 0; j < Nyr; j++) {
            double x = loc_x0 + i*dx;
            double y = loc_y0 - j*dy;
            double r = ComputeR(x, y);
            U0[i*Nyr+j] = (r <= 1.0)? pow(2.0*(1.0-r),4.0) * (4.0*r+1.0) : 0.0;
        }
    }

    // Gathering occurs only in writing to file
    // MPI_Allgatherv(loc_U0, Nyr*loc_Nxr, MPI_DOUBLE, U0, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
}

/**
 * Sets velocity field in x,y for U, V
 * */
void Burgers2P::SetIntegratedVelocity() {
    // Get model parameters
    int Nt = model->GetNt();

    // Set Matrix Coefficients
    SetMatrixCoefficients();


    // Compute U, V for every step k
    U[0] = U0;
    V[0] = U0;
    for (int k = 0; k < Nt-1; k++) {
        U[k+1] = NextVelocityState(U[k], V[k], true);
        V[k+1] = NextVelocityState(U[k], V[k], false);
    }
}

/**
 * Writes the velocity field for U, V into a file
 * IMPORTANT: Run SetIntegratedVelocity() first
 * */
void Burgers2P::WriteVelocityFile() {
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


void Burgers2P::SetEnergy() {
    // Get Model parameters
    int Nt = model->GetNt();

    // Calculate Energy
    for (int k = 0; k < Nt; k++) {
        E[k] = CalculateEnergyState(U[k], V[k]);
    }
}

double Burgers2P::CalculateEnergyState(double* Ui, double* Vi) {
    // MPI Parameters
    int p = model->GetP();
    int loc_rank = model->GetRank();
    int* sendcounts = new int[p];
    int* displs = new int[p];

    // Get Model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    // Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    // Split x domain into 2
    int loc_Nxr = (Nxr % 2 != 0 && loc_rank == 0) ? (Nxr/2)+1 : Nxr/2;
    sendcounts[0] = Nyr*((Nxr/2)+1); sendcounts[1] = Nyr*(Nxr/2);
    displs[0] = 0; displs[1] = sendcounts[0];

    double* loc_U = new double[loc_Nxr*Nyr];
    double* loc_V = new double[loc_Nxr*Nyr];
    double* U_dots = new double[p];
    double* V_dots = new double[p];

    MPI_Scatterv(Ui, sendcounts, displs, MPI_DOUBLE, loc_U, sendcounts[loc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(Vi, sendcounts, displs, MPI_DOUBLE, loc_V, sendcounts[loc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double loc_ddotU = F77NAME(ddot)(Nyr*loc_Nxr, loc_U, 1, loc_U, 1);
    double loc_ddotV = F77NAME(ddot)(Nyr*loc_Nxr, loc_V, 1, loc_V, 1);

    MPI_Allgather(&loc_ddotU, 1, MPI_DOUBLE, U_dots, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&loc_ddotV, 1, MPI_DOUBLE, V_dots, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    double NextEnergyState = 0.5 * (U_dots[0] + U_dots[1] + V_dots[0] + V_dots[1]);

    delete[] loc_U;
    delete[] loc_V;
    delete[] U_dots;
    delete[] V_dots;

    return NextEnergyState;
}

/**
 * Computes and returns next velocity state based on previous inputs
 * */
double* Burgers2P::NextVelocityState(double* Ui, double* Vi, bool SELECT_U) {
    // Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    double dt = model->GetDt();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double b = model->GetB();

    // Set aliases for computation
    double* Vel = (SELECT_U) ? Ui : Vi;
    double* Other = (SELECT_U)? Vi : Ui;

    // Set caches for Vel
    SetCaches(Vel);

    // No need for scatterv because you are working on a sub-part of the matrix
    // MPI_Scatterv(Vel, sendcounts, displs, MPI_DOUBLE, loc_Vel, sendcounts[loc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Scatterv(Other, sendcounts, displs, MPI_DOUBLE, loc_Other, sendcounts[loc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
    F77NAME(dsymm)('R', 'U', Nyr, Nxr, 1.0, dVel_dx_2_coeffs, Nxr, Vel, Nyr, 0.0, dVel_dx_2, Nyr);
    F77NAME(dsymm)('L', 'U', Nyr, Nxr, 1.0, dVel_dy_2_coeffs, Nyr, Vel, Nyr, 0.0, dVel_dy_2, Nyr);

    // Compute first derivatives
    F77NAME(dcopy)(Nyr*Nxr, Vel, 1, dVel_dx, 1);
    F77NAME(dcopy)(Nyr*Nxr, Vel, 1, dVel_dy, 1);
    F77NAME(dtrmm)('R', 'U', 'N', 'N', Nyr, Nxr, 1.0, dVel_dx_coeffs, Nxr, dVel_dx, Nyr);
    F77NAME(dtrmm)('L', 'L', 'N', 'N', Nyr, Nxr, 1.0, dVel_dy_coeffs, Nyr, dVel_dy, Nyr);

    // Compute b terms
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

    UpdateBoundsLinear(dVel_dx_2, dVel_dy_2, dVel_dx, dVel_dy);
    UpdateBoundsNonLinear(Vel, Other, Vel_Vel_Minus_1, Vel_Other_Minus_1, SELECT_U);

    // Matrix addition through all terms
    for (int i = 0; i < Nyr*Nxr; i++) {
        NextVel[i] = dVel_dx_2[i] + dVel_dy_2[i] - dVel_dx[i] - dVel_dy[i] -
                (Vel_Vel[i] + Vel_Other[i] - Vel_Vel_Minus_1[i] - Vel_Other_Minus_1[i]);
        NextVel[i] *= dt;
        NextVel[i] += Vel[i];
    }

    if (model->GetRank() == 1 && SELECT_U) {
        printDebug(Vel, Nyr, Nxr, 'U');
    }

    // Gathering occurs in writing file
    // MPI_Allgatherv(loc_NextVel, Nyr*Nxr, MPI_DOUBLE, NextVel, sendcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    // Delete term array pointers
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
void Burgers2P::SetMatrixCoefficients() {
    // Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double ax = model->GetAx();
    double ay = model->GetAy();
    double c = model->GetC();

    // Set coefficients
    GenSymm((-2.0*c)/pow(dx,2.0), c/pow(dx,2.0), Nxr, Nxr, dVel_dx_2_coeffs);
    GenSymm((-2.0*c)/pow(dy,2.0), c/pow(dy,2.0), Nyr, Nyr, dVel_dy_2_coeffs);
    GenTrmm(ax/dx, -ax/dx, Nxr, Nxr, true, dVel_dx_coeffs);
    GenTrmm(ay/dy, -ay/dy, Nyr, Nyr, false, dVel_dy_coeffs);
}

void Burgers2P::SetCaches(double* Vel) {
    // Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();

    // Get ranks
    int up = model->GetUp();
    int down = model->GetDown();
    int left = model->GetLeft();
    int right = model->GetRight();


    // Get communicator
    MPI_Comm vu = model->GetComm();
    int flag;

    // Allocate memory into Vel for local bounds
    double* myUpVel = new double[Nxr];
    double* myDownVel = new double[Nxr];
    double* myLeftVel = new double[Nyr];
    double* myRightVel = new double[Nyr];

    // Get Vel bounds for this sub-matrix
    for (int i = 0; i < Nxr; i++) {
        myUpVel[i] = Vel[i*Nyr];
        myDownVel[i] = Vel[i*Nyr+(Nyr-1)];
    }
    for (int j = 0; j < Nyr; j++) {
        myLeftVel[j] = Vel[j];
        myRightVel[j] = Vel[(Nxr-1)*Nyr+j];
    }

    // Exchange up/down
    flag = 0;
    /* Send down boundary to down and receive into up boundary */
    MPI_Sendrecv(myDownVel, Nxr, MPI_DOUBLE, down, flag, upVel, Nxr, MPI_DOUBLE, up, flag, vu, MPI_STATUS_IGNORE);
    /* Send up boundary to up and receive into down boundary */
    MPI_Sendrecv(myUpVel, Nxr, MPI_DOUBLE, up, flag, downVel, Nxr, MPI_DOUBLE, down, flag, vu, MPI_STATUS_IGNORE);

    // Exchange left/right
    flag = 1;
    /* Send right boundary to right and receive into left boundary */
    MPI_Sendrecv(myRightVel, Nyr, MPI_DOUBLE, right, flag, leftVel, Nyr, MPI_DOUBLE, left, flag, vu, MPI_STATUS_IGNORE);
    /* Send left boundary to left and receive into right boundary */
    MPI_Sendrecv(myLeftVel, Nyr, MPI_DOUBLE, left, flag, rightVel, Nyr, MPI_DOUBLE, right, flag, vu, MPI_STATUS_IGNORE);

    // Handle boundaries
    if (up < 0) SetZeroes(upVel, Nxr);
    if (down < 0) SetZeroes(downVel, Nxr);
    if (left < 0) SetZeroes(leftVel, Nyr);
    if (right < 0) SetZeroes(rightVel, Nyr);

    delete[] myUpVel;
    delete[] myDownVel;
    delete[] myLeftVel;
    delete[] myRightVel;
}

void Burgers2P::UpdateBoundsLinear(double* dVel_dx_2, double* dVel_dy_2, double* dVel_dx, double* dVel_dy) {
    // Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double ax = model->GetAx();
    double ay = model->GetAy();
    double c = model->GetC();

    // Modify boundaries on this sub-matrix

    /* Fix left and right boundaries */
    for (int j = 0; j < Nyr; j++) {
        // left
        dVel_dx_2[j] += c/pow(dx,2.0)*leftVel[j];
        dVel_dx[j] -= ax/dx*leftVel[j];

        // right
        dVel_dx_2[(Nxr-1)*Nyr+j] += c/pow(dx,2.0)*rightVel[j];
    }

    /* Fix up and down boundaries */
    for (int i = 0; i < Nxr; i++) {
        // up
        dVel_dy_2[i*Nyr] += c/pow(dy,2.0)*upVel[i];
        dVel_dy[i*Nyr] -= ay/dy*upVel[i];

        // down
        dVel_dy_2[i*Nyr+(Nyr-1)] += c/pow(dy,2.0)*downVel[i];
    }
}

void Burgers2P::UpdateBoundsNonLinear(double* Vel, double* Other, double* Vel_Vel_Minus_1, double* Vel_Other_Minus_1, bool SELECT_U) {
    // Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double b = model->GetB();

    if (SELECT_U) {
        //up
        for (int i = 0; i < Nxr; i++) {
            Vel_Other_Minus_1[i*Nyr] = b/dy * upVel[i] * Other[i*Nyr];
        }
        // left
        for (int j = 0; j < Nyr; j++) {
            Vel_Vel_Minus_1[j] = b/dx * leftVel[j] * Vel[j];
        }
    }
    else {
        // up
        for (int i = 0; i < Nxr; i++) {
            Vel_Vel_Minus_1[i*Nyr] = b/dy * upVel[i] * Vel[i*Nyr];
        }
        // left
        for (int j = 0; j < Nyr; j++) {
            Vel_Other_Minus_1[j] = b/dx * leftVel[j] * Other[j];
        }
    }
}

double Burgers2P::ComputeR(double x, double y) {
    double r = pow(x*x+y*y, 0.5);
    return r;
}
