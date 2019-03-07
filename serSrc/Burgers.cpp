#include <cmath>
#include <fstream>
#include <iomanip>
#include "BLAS_Wrapper.h"
#include "Helpers.h"
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
        double* NextU = NextVelocityState(true);
        double* NextV = NextVelocityState(false);
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
 * @param Ui U velocity per timestamp
 * @param Vi V velocity per timestamp
 * @param SELECT_U true if the computation is for U
 * */
double* Burgers::NextVelocityState(bool SELECT_U) {
    /// Get model parameters
    int Ny = model->GetNy();
    int Nx = model->GetNx();

    /// Reduced parameters
    int Nyr = Ny - 2;
    int Nxr = Nx - 2;

    /// Set aliases for computation
    double* Vel = (SELECT_U) ? U : V;
    double* Other = (SELECT_U)? V : U;

    /// Generate NextVel
    double* NextVel = new double[Nyr*Nxr];

    SetLinearTerms(Vel, NextVel);
    SetNonLinearTerms(Vel, Other, NextVel, SELECT_U);

    return NextVel;
}

void Burgers::SetLinearTerms(double* Vel, double* NextVel) {
    int Nxr = model->GetNx() - 2;
    int Nyr = model->GetNy() - 2;
    double alpha_dx_2 = model->GetAlphaDx_2();
    double beta_dx_2 = model->GetBetaDx_2();
    double alpha_dy_2 = model->GetAlphaDy_2();
    double beta_dy_2 = model->GetBetaDy_2();
    double alpha_dx_1 = model->GetAlphaDx_1();
    double beta_dx_1 = model->GetBetaDx_1();
    double alpha_dy_1 = model->GetAlphaDy_1();
    double beta_dy_1 = model->GetBetaDy_1();

    // loop blocking + pre-fetching previous & next column from memory
    const int blocksize = 8;
    double* Vel_iMinus = nullptr;
    double* Vel_iPlus = nullptr;
    double* VelPtr = nullptr;
    double* NextVelPtr = nullptr;
    int i, j, i2, j2, ii2, jj2;
    for (i = 0; i < Nxr; i+=blocksize) {
        for (j = 0; j < Nyr; j+=blocksize) {
            for (i2 = 0, ii2 = i+i2, VelPtr = &Vel[i*Nyr+j], NextVelPtr = &NextVel[i*Nyr+j]; ii2 < Nxr && i2 < blocksize; ++i2, ++ii2, VelPtr += Nyr, NextVelPtr+=Nyr) {
                if (ii2 > 0) Vel_iMinus = &Vel[(ii2-1)*Nyr+j];
                if (ii2 < Nxr-1) Vel_iPlus = &Vel[(ii2+1)*Nyr+j];
                for (j2 = 0, jj2 = j+j2; jj2 < Nyr && j2 < blocksize; ++j2, ++jj2) {
                    // Update x
                    NextVelPtr[j2] = (ii2 > 0) ? alpha_dx_1 * VelPtr[j2] + beta_dx_1 * Vel_iMinus[j2] :alpha_dx_1 * VelPtr[j2];
                    NextVelPtr[j2] = (ii2 > 0) ? NextVelPtr[j2] + alpha_dx_2 * VelPtr[j2] + beta_dx_2 * Vel_iMinus[j2] :NextVelPtr[j2] + alpha_dx_2 * VelPtr[j2];
                    NextVelPtr[j2] = (ii2 < Nxr-1) ? NextVelPtr[j2] + beta_dx_2 * Vel_iPlus[j2] : NextVelPtr[j2];
                    // Update y
                    NextVelPtr[j2] = (jj2 > 0) ? NextVelPtr[j2] + alpha_dy_1 * VelPtr[j2] + beta_dy_1 * VelPtr[j2-1] : NextVelPtr[j2] + alpha_dy_1 * VelPtr[j2];
                    NextVelPtr[j2] = (jj2 > 0) ? NextVelPtr[j2] + alpha_dy_2 * VelPtr[j2] + beta_dy_2 * VelPtr[j2-1] : NextVelPtr[j2] + alpha_dy_2 * VelPtr[j2];
                    NextVelPtr[j2] = (jj2 < Nyr-1) ? NextVelPtr[j2] + beta_dy_2 * VelPtr[j2+1] : NextVelPtr[j2];
                }
            }
        }
    }
//    for (int i = 0; i < Nxr; i++) {
//        if (i > 0) Vel_iMinus = &(Vel[(i-1)*Nyr]);
//        if (i < Nxr-1) Vel_iPlus = &(Vel[(i+1)*Nyr]);
//        int start = i*Nyr;
//        for (int j = 0; j < Nyr; j+=blocksize) {
//            for (int k = j; k < Nyr && k < j + blocksize; k++) {
//                int curr = start + k;
//                // Update x
//                NextVel[curr] = (i > 0) ? alpha_dx_1 * Vel[curr] + beta_dx_1 * Vel_iMinus[k] :alpha_dx_1 * Vel[curr];
//                NextVel[curr] = (i > 0) ? NextVel[curr] + alpha_dx_2 * Vel[curr] + beta_dx_2 * Vel_iMinus[k] :NextVel[curr] + alpha_dx_2 * Vel[curr];
//                NextVel[curr] = (i < Nxr-1) ? NextVel[curr] + beta_dx_2 * Vel_iPlus[k] : NextVel[curr];
//                // Update y
//                NextVel[curr] = (k > 0) ? NextVel[curr] + alpha_dy_1 * Vel[curr] + beta_dy_1 * Vel[curr-1] : NextVel[curr] + alpha_dy_1 * Vel[curr];
//                NextVel[curr] = (k > 0) ? NextVel[curr] + alpha_dy_2 * Vel[curr] + beta_dy_2 * Vel[curr-1] : NextVel[curr] + alpha_dy_2 * Vel[curr];
//                NextVel[curr] = (k < Nyr-1) ? NextVel[curr] + beta_dy_2 * Vel[curr+1] : NextVel[curr];
//            }
//        }
//    }
}

void Burgers::SetNonLinearTerms(double* Vel, double* Other, double* NextVel, bool SELECT_U) {
    /// Get model parameters
    int Nyr = model->GetNy() - 2;
    int Nxr = model->GetNx() - 2;
    double dt = model->GetDt();
    double bdx = model->GetBDx();
    double bdy = model->GetBDy();
    const int blocksize = 4;

    double* Vel_iMinus = nullptr;
    double Vel_Vel, Vel_Other, Vel_Vel_Minus_1, Vel_Other_Minus_1;
    if (SELECT_U) {
        for (int i = 0; i < Nxr; i++) {
            if (i > 0) Vel_iMinus = &(Vel[(i-1)*Nyr]);
            int start = i*Nyr;
            for (int j = 0; j < Nyr; j+=blocksize) {
                for (int k = j; k < Nyr && k < j + blocksize; k++) {
                    int curr = start + k;
                    Vel_Vel = bdx * Vel[curr] * Vel[curr];
                    Vel_Other = bdy * Vel[curr] * Other[curr];
                    Vel_Vel_Minus_1 = (i == 0) ? 0 : bdx * Vel_iMinus[k] * Vel[curr];
                    Vel_Other_Minus_1 = (k == 0) ? 0 : bdy * Vel[curr - 1] * Other[curr];
                    NextVel[curr] -= (Vel_Vel + Vel_Other - Vel_Vel_Minus_1 - Vel_Other_Minus_1);
                    NextVel[curr] *= dt;
                    NextVel[curr] += Vel[curr];
                }
            }
        }
    }
    else {
        for (int i = 0; i < Nxr; i++) {
            if (i > 0) Vel_iMinus = &(Vel[(i-1)*Nyr]);
            int start = i*Nyr;
            for (int j = 0; j < Nyr; j+=blocksize) {
                for (int k = j; k < Nyr && k < j + blocksize; k++) {
                    int curr = start + k;
                    Vel_Vel = bdy * Vel[curr] * Vel[curr];
                    Vel_Other = bdx * Vel[curr] * Other[curr];
                    Vel_Vel_Minus_1 = (k == 0) ? 0 : bdy * Vel[curr - 1] * Vel[curr];
                    Vel_Other_Minus_1 = (i == 0) ? 0 : bdx * Vel_iMinus[k] * Other[curr];
                    NextVel[curr] -= (Vel_Vel + Vel_Other - Vel_Vel_Minus_1 - Vel_Other_Minus_1);
                    NextVel[curr] *= dt;
                    NextVel[curr] += Vel[curr];
                }
            }
        }
    }
}
