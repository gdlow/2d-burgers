#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include "BLAS_Wrapper.h"
#include "Burgers2P.h"

using namespace std;

/**
 * @brief Public Constructor: Accepts a Model instance reference as input
 * Allocates memory to all other instance variables
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
    NextU = new double[NyrNxr];
    NextV = new double[NyrNxr];

    /// Caches
    upV = new double[Nxr];
    downV = new double[Nxr];
    leftV = new double[Nyr];
    rightV = new double[Nyr];
    myUpV = new double[Nxr];
    myDownV = new double[Nxr];
    myLeftV = new double[Nyr];
    myRightV = new double[Nyr];
    upU = new double[Nxr];
    downU = new double[Nxr];
    leftU = new double[Nyr];
    rightU = new double[Nyr];
    myUpU = new double[Nxr];
    myDownU = new double[Nxr];
    myLeftU = new double[Nyr];
    myRightU = new double[Nyr];

    /// Generate new MPI request and stats
    reqs = new MPI_Request[16];
    stats = new MPI_Status[16];
}

/**
 * @brief Destructor: Deletes all allocated pointers in the class instance
 * */
Burgers2P::~Burgers2P() {
    /// Delete U and V
    delete[] U;
    delete[] V;
    delete[] NextU;
    delete[] NextV;

    /// Delete Caches
    delete[] upV;
    delete[] downV;
    delete[] leftV;
    delete[] rightV;
    delete[] myUpV;
    delete[] myDownV;
    delete[] myLeftV;
    delete[] myRightV;
    delete[] upU;
    delete[] downU;
    delete[] leftU;
    delete[] rightU;
    delete[] myUpU;
    delete[] myDownU;
    delete[] myLeftU;
    delete[] myRightU;

    /// Deallocate memory of MPI requests and stats
    delete[] stats;
    delete[] reqs;

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
            double r = pow(x*x+y*y, 0.5);
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

    /// Compute U, V for every step k
    double* temp = nullptr;
    for (int k = 0; k < Nt-1; k++) {
        GetNextVelocities();

        temp = NextU;
        NextU = U;
        U = temp;

        temp = NextV;
        NextV = V;
        V = temp;
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
 * @brief Calculates and sets energy of velocity field
 * */
void Burgers2P::SetEnergy() {
    E = CalculateEnergyState(U, V);
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
 * */
void Burgers2P::GetNextVelocities() {
    int NyrNxr = model->GetLocNyrNxr();
    SetCaches();
    ComputeNextVelocityState();
    MPI_Waitall(16, reqs, stats);
    FixNextVelocityBoundaries();
    for (int k = 0; k < NyrNxr; k++) {
        NextU[k] += U[k];
        NextV[k] += V[k];
    }
}

/**
 * @brief Private helper function that sets the boundary condition velocities
 * */
void Burgers2P::SetCaches() {
    /// Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();
    int NyrNxr = model->GetLocNyrNxr();

    /// Get ranks
    int up = model->GetUp();
    int down = model->GetDown();
    int left = model->GetLeft();
    int right = model->GetRight();

    /// Get communicator
    MPI_Comm vu = model->GetComm();
    int flag;

    /// Get Vel bounds for this sub-matrix
    for (int k = 0, i = 0; k < NyrNxr; k += Nyr, i++) {
        myUpU[i] = U[k];
        myUpV[i] = V[k];
        int didx = k + Nyr-1;
        myDownU[i] = U[didx];
        myDownV[i] = V[didx];
    }
    for (int k = (Nxr-1)*Nyr, i = 0; k < NyrNxr; k++, i++) {
        myLeftU[i] = U[i];
        myLeftV[i] = V[i];
        myRightU[i] = U[k];
        myRightV[i] = V[k];
    }

    /// Exchange up/down
    flag = 0;
    /* Send down boundary to down and receive into up boundary */
    MPI_Isend(myDownU, Nxr, MPI_DOUBLE, down, flag, vu, &reqs[0]);
    MPI_Irecv(upU, Nxr, MPI_DOUBLE, up, flag, vu, &reqs[1]);
    MPI_Isend(myDownV, Nxr, MPI_DOUBLE, down, flag, vu, &reqs[2]);
    MPI_Irecv(upV, Nxr, MPI_DOUBLE, up, flag, vu, &reqs[3]);
    /* Send up boundary to up and receive into down boundary */
    MPI_Isend(myUpU, Nxr, MPI_DOUBLE, up, flag, vu, &reqs[4]);
    MPI_Irecv(downU, Nxr, MPI_DOUBLE, down, flag, vu, &reqs[5]);
    MPI_Isend(myUpV, Nxr, MPI_DOUBLE, up, flag, vu, &reqs[6]);
    MPI_Irecv(downV, Nxr, MPI_DOUBLE, down, flag, vu, &reqs[7]);

    /// Exchange left/right
    flag = 1;
    /* Send right boundary to right and receive into left boundary */
    MPI_Isend(myRightU, Nyr, MPI_DOUBLE, right, flag, vu, &reqs[8]);
    MPI_Irecv(leftU, Nyr, MPI_DOUBLE, left, flag, vu, &reqs[9]);
    MPI_Isend(myRightV, Nyr, MPI_DOUBLE, right, flag, vu, &reqs[10]);
    MPI_Irecv(leftV, Nyr, MPI_DOUBLE, left, flag, vu, &reqs[11]);
    /* Send left boundary to left and receive into right boundary */
    MPI_Isend(myLeftU, Nyr, MPI_DOUBLE, left, flag, vu, &reqs[12]);
    MPI_Irecv(rightU, Nyr, MPI_DOUBLE, right, flag, vu, &reqs[13]);
    MPI_Isend(myLeftV, Nyr, MPI_DOUBLE, left, flag, vu, &reqs[14]);
    MPI_Irecv(rightV, Nyr, MPI_DOUBLE, right, flag, vu, &reqs[15]);
}

/**
 * @brief Computes linear and non-linear terms for U and V
 * */
void Burgers2P::ComputeNextVelocityState() {
    /// Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();

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
}

/**
 * @brief Fixes boundary conditions for U and V
 * */
void Burgers2P::FixNextVelocityBoundaries() {
    /// Get model parameters
    int Nyr = model->GetLocNyr();
    int Nxr = model->GetLocNxr();

    /// Get ranks
    int up = model->GetUp();
    int down = model->GetDown();
    int left = model->GetLeft();
    int right = model->GetRight();
    double beta_dx_sum = model->GetBetaDx_Sum();
    double beta_dy_sum = model->GetBetaDy_Sum();
    double beta_dx_2 = model->GetBetaDx_2();
    double beta_dy_2 = model->GetBetaDy_2();
    double bdx = model->GetBDx();
    double bdy = model->GetBDy();

    /// Fix left and right boundaries
    double bdxU, bdyV;
    for (int j = 0; j < Nyr; j++) {
        if (left >= 0) {
            bdxU = bdx * U[j];
            double bdxU_total = beta_dx_sum + bdxU;
            NextU[j] += bdxU_total * leftU[j];
            NextV[j] += bdxU_total * leftV[j];
        }
        if (right >= 0) {
            int ridx = (Nxr-1)*Nyr+j;
            NextU[ridx] += beta_dx_2 * rightU[j];
            NextV[ridx] += beta_dx_2 * rightV[j];
        }
    }

    /// Fix up and down boundaries
    for (int i = 0; i < Nxr; i++) {
        int upidx = i*Nyr;
        if (up >= 0) {
            bdyV = bdy * V[upidx];
            double bdyV_total = beta_dy_sum + bdyV;
            NextU[upidx] += bdyV_total * upU[i];
            NextV[upidx] += bdyV_total * upV[i];
        }
        if (down >= 0) {
            int didx = upidx + (Nyr-1);
            NextU[didx] += beta_dy_2 * downU[i];
            NextV[didx] += beta_dy_2 * downV[i];

        }
    }
}

/**
 * @brief Private helper function that assembles the global matrix into a pre-allocated M
 * Arranges data into row-major format from a column-major format in the 1D pointer Vel
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
