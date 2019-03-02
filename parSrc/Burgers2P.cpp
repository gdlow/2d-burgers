#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include "BLAS_Wrapper.h"
#include "Helpers.h"
#include "Burgers2P.h"
#include <iostream>

using namespace std;

/**
 * @brief Public Constructor: Accepts a Model instance reference as input
 * @brief Allocates memory to all other instance variables
 * @param &m reference to Model instance
 * */
Burgers2P::Burgers2P(Model &m) {
    /// Set model class pointer as instance variable
    model = &m;

    /// Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    int NyrNxr = model->GetLocNyrNxr();

    /// Allocate memory to instance variables
    U = new double[NyrNxr];
    V = new double[NyrNxr];

    /// Term arrays
    dVel_2 = new double[NyrNxr];
    dVel = new double[NyrNxr];

    /// Matrix coefficients
    dVel_dx_2_coeffs = new double[Nxr*Nxr];
    dVel_dy_2_coeffs = new double[Nyr*Nyr];
    dVel_dx_coeffs = new double[Nxr*Nxr]();
    dVel_dy_coeffs = new double[Nyr*Nyr]();

    /// Caches
    upVel = new double[Nxr];
    downVel = new double[Nxr];
    leftVel = new double[Nyr];
    rightVel = new double[Nyr];
}

/**
 * @brief Destructor: Deletes all allocated pointers in the class instance
 * */
Burgers2P::~Burgers2P() {
    /// Delete U and V
    delete[] U;
    delete[] V;

    /// Delete term arrays
    delete[] dVel_2;
    delete[] dVel;

    /// Delete Caches
    delete[] upVel;
    delete[] downVel;
    delete[] leftVel;
    delete[] rightVel;

    /// Delete matrix coefficients
    delete[] dVel_dx_2_coeffs;
    delete[] dVel_dy_2_coeffs;
    delete[] dVel_dx_coeffs;
    delete[] dVel_dy_coeffs;

    /// model is not dynamically alloc
}

/**
 * @brief Sets initial velocity field in x,y for U0 (V0 = U0)
 * */
void Burgers2P::SetInitialVelocity() {
    /// Get model parameters
    double x0 = model->GetX0();
    double y0 = model->GetY0();
    double dx = model->GetDx();
    double dy = model->GetDy();
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    int displ_x = model->GetDisplX();
    int displ_y = model->GetDisplY();

    /// Compute U0
    /// Memory layout in column-major format
    double loc_x0 = x0 + (displ_x+1)*dx;
    double loc_y0 = y0 - (displ_y+1)*dy;
    for (int i = 0; i < Nxr; i++) {
        for (int j = 0; j < Nyr; j++) {
            double x = loc_x0 + i*dx;
            double y = loc_y0 - j*dy;
            double r = ComputeR(x, y);
            U[i*Nyr+j] = (r <= 1.0)? 2.0*pow(1.0-r,4.0) * (4.0*r+1.0) : 0.0;
            V[i*Nyr+j] = (r <= 1.0)? 2.0*pow(1.0-r,4.0) * (4.0*r+1.0) : 0.0;
        }
    }
}

/**
 * @brief Sets velocity field in x,y for U, V
 * */
void Burgers2P::SetIntegratedVelocity() {
    /// Get model parameters
    int Nt = model->GetNt();

    /// Set Matrix Coefficients
    SetMatrixCoefficients();


    /// Compute U, V for every step k
    for (int k = 0; k < Nt-1; k++) {
        double* NextU = NextVelocityState(U, V, true);
        double* NextV = NextVelocityState(U, V, false);
        CopyAndDelete(NextU, NextV);
        if (model->GetRank() == 0) cout << "step: " << k << "\n";
    }
}

/**
 * @brief Writes the velocity field for U, V into a file
 * IMPORTANT: Run SetIntegratedVelocity() first
 * */
void Burgers2P::WriteVelocityFile() {
    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    /// Allocate 2D pointer
    double** M = new double*[Ny-2];
    for (int j = 0; j < Ny-2; j++) {
        M[j] = new double[Nx-2];
    }

    /// Open output file stream to data.txt
    ofstream of;
    of.open("data.txt", ios::out | ios::trunc);
    of.precision(4); // 4 s.f.

    /// Write U velocity
    WriteOf(U, M, of, 'U');

    /// Write V velocity
    WriteOf(V, M, of, 'V');
    of.close();

    /// Delete 2D pointer
    for (int j = 0; j < Ny-2; j++) {
        delete[] M[j];
    }
    delete[] M;
}

/**
 * @brief Private helper function to write to output stream
 * @param Vel pointer to either U or V
 * @param M 2D pointer representing global matrix (should have been allocated memory)
 * @param &of reference to output file stream
 * @param id Supply 'U' or 'V'
 * */
void Burgers2P::WriteOf(double* Vel, double** M, ofstream &of, char id) {
    int loc_rank = model->GetRank();
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    AssembleMatrix(Vel, M);
    if (loc_rank == 0) {
        of << id << " velocity field:" << endl;
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                if (j == 0 || i == 0 || j == Ny - 1 || i == Nx - 1) {
                    of << 0 << ' ';
                } else {
                    of << M[j - 1][i - 1] << ' ';
                }
            }
            of << endl;
        }
    }
}

/**
 * @brief Calculates and sets energy of velocity field
 * */
void Burgers2P::SetEnergy() {
    E = CalculateEnergyState(U, V);
}

/**
 * @brief Private helper function that calculates next energy state per timestamp
 * @param Ui U velocity per timestamp (i.e. supply U[k])
 * @param Vi V velocity per timestamp (i.e. supply V[k])
 * */
double Burgers2P::CalculateEnergyState(double* Ui, double* Vi) {
    /// Get model parameters
    int NyrNxr = model->GetLocNyrNxr();
    double dx = model->GetDx();
    double dy = model->GetDy();
    MPI_Comm vu = model->GetComm();

    /// Blas calls to compute dot products
    double loc_ddotU = F77NAME(ddot)(NyrNxr, Ui, 1, Ui, 1);
    double loc_ddotV = F77NAME(ddot)(NyrNxr, Vi, 1, Vi, 1);

    /// Compute local energy state
    double NextLocalEnergyState = 0.5 * (loc_ddotU + loc_ddotV) * dx*dy;
    double NextGlobalEnergyState;

    /// Sum into global energy state
    MPI_Allreduce(&NextLocalEnergyState, &NextGlobalEnergyState, 1, MPI_DOUBLE, MPI_SUM, vu);
    return NextGlobalEnergyState;
}

/**
 * @brief Private helper function that computes and returns next velocity state based on previous inputs
 * @param Ui U velocity per timestamp (i.e. supply U[k])
 * @param Vi V velocity per timestamp (i.e. supply V[k])
 * @param SELECT_U true if the computation is for U
 * */
double* Burgers2P::NextVelocityState(double* Ui, double* Vi, bool SELECT_U) {
    /// Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    int NyrNxr = model->GetLocNyrNxr();
    double dt = model->GetDt();
    double bdx = model->GetBDx();
    double bdy = model->GetBDy();

    /// Get ranks
    int up = model->GetUp();
    int left = model->GetLeft();

    /// Set aliases for computation
    double* Vel = (SELECT_U) ? Ui : Vi;
    double* Other = (SELECT_U)? Vi : Ui;

    /// Create MPI requests and stats for this sub-matrix
    MPI_Request* reqs = new MPI_Request[8];
    MPI_Status* stats = new MPI_Status[8];

    /// Set caches for Vel (Non-blocking)
    SetCaches(Vel, reqs);

    /// Generate term arrays
    double* NextVel = new double[NyrNxr];

    /// Compute first & second derivatives
    // x
    for (int i = 0; i < Nyr; i++) {
        F77NAME(dgbmv)('N', Nxr, Nxr, 1, 1, 1.0, dVel_dx_2_coeffs, Nxr, &(Vel[i]), Nyr, 0.0, &(dVel_2[i]), Nyr);
        F77NAME(dgbmv)('N', Nxr, Nxr, 1, 0, 1.0, dVel_dx_coeffs, Nxr, &(Vel[i]), Nyr, 0.0, &(dVel[i]), Nyr);
    }
    // y
    for (int i = 0; i < NyrNxr; i += Nyr) {
        F77NAME(dgbmv)('N', Nyr, Nyr, 1, 1, 1.0, dVel_dy_2_coeffs, Nyr, &(Vel[i]), 1, 1.0, &(dVel_2[i]), 1);
        F77NAME(dgbmv)('N', Nyr, Nyr, 1, 0, 1.0, dVel_dy_coeffs, Nyr, &(Vel[i]), 1, 1.0, &(dVel[i]), 1);
    }

    UpdateBoundsLinear(dVel_2, dVel_2, dVel, dVel, reqs, stats);

    /// Matrix addition through all terms
    if (SELECT_U) {
        for (int i = 0; i < NyrNxr; i++) {
            double Vel_Vel = bdx * Vel[i] * Vel[i];
            double Vel_Other = bdy * Vel[i] * Other[i];
            double Vel_Vel_Minus_1 = (i < Nyr)? 0 : bdx * Vel[i-Nyr] * Vel[i];
            double Vel_Other_Minus_1 = (i % Nyr == 0)? 0 : bdy * Vel[i-1] * Other[i];

            // Update non-linear BC
            if (i < Nyr && left >= 0) Vel_Vel_Minus_1 = bdx * leftVel[i] * Vel[i];
            if (i % Nyr == 0 && up >= 0) Vel_Other_Minus_1 = bdy * upVel[i/Nyr] * Other[i];

            NextVel[i] = dVel_2[i] - dVel[i] -
                         (Vel_Vel + Vel_Other - Vel_Vel_Minus_1 - Vel_Other_Minus_1);
            NextVel[i] *= dt;
            NextVel[i] += Vel[i];
        }
    }
    else {
        for (int i = 0; i < NyrNxr; i++) {
            double Vel_Vel = bdy * Vel[i] * Vel[i];
            double Vel_Other = bdx * Vel[i] * Other[i];
            double Vel_Vel_Minus_1 = (i % Nyr == 0)? 0 : bdy * Vel[i-1] * Vel[i];
            double Vel_Other_Minus_1 = (i < Nyr)? 0 : bdx * Vel[i-Nyr] * Other[i];

            // Update non-linear BC
            if (i < Nyr && left >= 0) Vel_Other_Minus_1 = bdx * leftVel[i] * Other[i];
            if (i % Nyr == 0 && up >= 0) Vel_Vel_Minus_1 = bdy * upVel[i/Nyr] * Vel[i];

            NextVel[i] = dVel_2[i] - dVel[i] -
                    (Vel_Vel + Vel_Other - Vel_Vel_Minus_1 - Vel_Other_Minus_1);
            NextVel[i] *= dt;
            NextVel[i] += Vel[i];
        }
    }
    /// Delete term array pointers
    delete[] reqs;
    delete[] stats;

    return NextVel;
}

/**
 * @brief Copies next state velocity into previous state, and deletes the temporary pointer
 * */
void Burgers2P::CopyAndDelete(double* NextU, double* NextV) {
    /// Delete current U and V pointers
    delete[] U;
    delete[] V;
    U = nullptr;
    V = nullptr;

    /// Set U and V to pointers to NextU and NextV
    U = NextU;
    V = NextV;
}
/**
 * @brief Private helper function that sets matrix coefficients for differentials
 * */
void Burgers2P::SetMatrixCoefficients() {
    /// Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
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
    GenTrmmBanded(alpha_dx_1, beta_dx_1, Nxr, false, dVel_dx_coeffs);
    GenTrmmBanded(alpha_dy_1, beta_dy_1, Nyr, false, dVel_dy_coeffs);
}

/**
 * @brief Private helper function that sets the boundary condition velocities
 * @param Vel pointer to U[k] or V[k]
 * */
void Burgers2P::SetCaches(double* Vel, MPI_Request* reqs) {
    /// Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();

    /// Get ranks
    int up = model->GetUp();
    int down = model->GetDown();
    int left = model->GetLeft();
    int right = model->GetRight();

    /// Get communicator
    MPI_Comm vu = model->GetComm();
    int flag;

    /// Allocate memory into Vel for local bounds
    double* myUpVel = new double[Nxr];
    double* myDownVel = new double[Nxr];
    double* myLeftVel = new double[Nyr];
    double* myRightVel = new double[Nyr];

    /// Get Vel bounds for this sub-matrix
    F77NAME(dcopy)(Nxr, Vel, Nyr, myUpVel, 1);
    F77NAME(dcopy)(Nxr, &(Vel[Nyr-1]), Nyr, myDownVel, 1);
    F77NAME(dcopy)(Nyr, Vel, 1, myLeftVel, 1);
    F77NAME(dcopy)(Nyr, &(Vel[(Nxr-1)*Nyr]), 1, myRightVel, 1);

    /// Exchange up/down
    flag = 0;
    /* Send down boundary to down and receive into up boundary */
    MPI_Isend(myDownVel, Nxr, MPI_DOUBLE, down, flag, vu, &reqs[0]);
    MPI_Irecv(upVel, Nxr, MPI_DOUBLE, up, flag, vu, &reqs[1]);
    /* Send up boundary to up and receive into down boundary */
    MPI_Isend(myUpVel, Nxr, MPI_DOUBLE, up, flag, vu, &reqs[2]);
    MPI_Irecv(downVel, Nxr, MPI_DOUBLE, down, flag, vu, &reqs[3]);

    /// Exchange left/right
    flag = 1;
    /* Send right boundary to right and receive into left boundary */
    MPI_Isend(myRightVel, Nyr, MPI_DOUBLE, right, flag, vu, &reqs[4]);
    MPI_Irecv(leftVel, Nyr, MPI_DOUBLE, left, flag, vu, &reqs[5]);
    /* Send left boundary to left and receive into right boundary */
    MPI_Isend(myLeftVel, Nyr, MPI_DOUBLE, left, flag, vu, &reqs[6]);
    MPI_Irecv(rightVel, Nyr, MPI_DOUBLE, right, flag, vu, &reqs[7]);

    delete[] myUpVel;
    delete[] myDownVel;
    delete[] myLeftVel;
    delete[] myRightVel;
}

/**
 * @brief Private helper function that updates the linear boundary conditions for the program's sub-matrix
 * @param dVel_dx_2 pointer to dVel_dx_2 in NextVelocityState()
 * @param dVel_dy_2 pointer to dVel_dy_2 in NextVelocityState()
 * @param dVel_dx pointer to dVel_dx in NextVelocityState()
 * @param dVel_dy pointer to dVel_dy in NextVelocityState()
 * */
void Burgers2P::UpdateBoundsLinear(double* dVel_dx_2, double* dVel_dy_2, double* dVel_dx, double* dVel_dy, MPI_Request* reqs, MPI_Status* stats) {
    /// Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    double beta_dx_2 = model->GetBetaDx_2();
    double beta_dx_1 = model->GetBetaDx_1();
    double beta_dy_2 = model->GetBetaDy_2();
    double beta_dy_1 = model->GetBetaDy_1();

    /// Get ranks
    int up = model->GetUp();
    int down = model->GetDown();
    int left = model->GetLeft();
    int right = model->GetRight();

    /// MPI wait for all comms to finish
    MPI_Waitall(8, reqs, stats);

    /// Modify boundaries on this sub-matrix

    /// Fix left and right boundaries
    for (int j = 0; j < Nyr; j++) {
        if (left >= 0) {
            dVel_dx_2[j] += beta_dx_2*leftVel[j];
            dVel_dx[j] += beta_dx_1*leftVel[j];
        }
        if (right >= 0) {
            dVel_dx_2[(Nxr-1)*Nyr+j] += beta_dx_2*rightVel[j];
        }
    }

    /// Fix up and down boundaries
    for (int i = 0; i < Nxr; i++) {
        if (up >= 0) {
            dVel_dy_2[i*Nyr] += beta_dy_2*upVel[i];
            dVel_dy[i*Nyr] += beta_dy_1*upVel[i];
        }
        if (down >= 0) {
            dVel_dy_2[i*Nyr+(Nyr-1)] += beta_dy_2*downVel[i];
        }
    }
}

/**
 * @brief Private helper function that assembles the global matrix into a preallocated M
 * @brief Arranges data into row-major format from a column-major format in the 1D pointer Vel
 * @param Vel 1D pointer to Vel in column-major format
 * @param M 2D pointer (pre-allocated memory) to be filled in row-major format
 * */
void Burgers2P::AssembleMatrix(double* Vel, double** M) {
    /// Get model parameters
    int loc_rank = model->GetRank();
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    int Px = model->GetPx();
    int Py = model->GetPy();
    MPI_Comm vu = model->GetComm();

    /// Don't delete these pointers (Part of Model object)
    int* displs = model->GetDispls();
    int* recvcount = model->GetRecvCount();
    int* rankNxrMap = model->GetRankNxrMap();
    int* rankNyrMap = model->GetRankNyrMap();
    int* rankDisplsXMap = model->GetRankDisplsXMap();
    int* rankDisplsYMap = model->GetRankDisplsYMap();

    /// Gather into globalVel in root (rank == 0)
    double* globalVel = new double[(Ny-2)*(Nx-2)];
    MPI_Gatherv(Vel, Nyr*Nxr, MPI_DOUBLE, globalVel, recvcount, displs, MPI_DOUBLE, 0, vu);

    /// Build global matrix in root, convert column-major -> row-major format
    if (loc_rank == 0) {
        for (int k = 0; k < Px*Py; k++) {
            for (int i = 0; i < rankNxrMap[k]; i++) {
                for (int j = 0; j < rankNyrMap[k]; j++) {
                    int loc_displ_y = rankDisplsYMap[k];
                    int loc_displ_x = rankDisplsXMap[k];
                    int loc_Nyr = rankNyrMap[k];
                    int global_displ = displs[k];
                    M[loc_displ_y+j][loc_displ_x+i] = globalVel[global_displ+i*loc_Nyr+j];
                }
            }
        }
    }

    delete[] globalVel;
}

/**
 * @brief Private helper function computing R for SetInitialVelocity()
 * @param x real x distance to origin
 * @param y real y distance to origin
 * */
double Burgers2P::ComputeR(double x, double y) {
    double r = pow(x*x+y*y, 0.5);
    return r;
}