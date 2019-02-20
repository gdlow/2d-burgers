#include <cmath>
#include <fstream>
#include <iomanip>
#include "Burgers.h"

using namespace std;

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
    int Ny = model->GetNy(); int Nx = model->GetNx();
    double x0 = model->GetX0(); double y0 = model->GetY0();
    double dx = model->GetDx(); double dy = model->GetDy();

    // Compute U0;
    U0 = nullptr;
    U0 = new double[Ny*Nx];
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // Assumes x0 and y0 are identifying top LHS of matrix
            double y = y0 - j*dy;
            double x = x0 + i*dx;
            double r = ComputeR(x, y);
            // Store in column-major format
            U0[i*Ny+j] = (r <= 1.0)? pow(2.0*(1.0-r),4.0) * (4.0*r-1.0) : 0.0;
        }
    }
}

/**
 * Sets velocity field in x,y for U, V
 * */
void Burgers::SetIntegratedVelocity() {
    // Get model parameters
    int Nt = model->GetNt();
    int Ny = model->GetNy();
    int Nx = model->GetNx();
    double dx = model->GetDx();
    double dy = model->GetDy();
    double dt = model->GetDt();
    double ax = model->GetAx();
    double ay = model->GetAy();
    double b = model->GetB();
    double c = model->GetC();

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

double Burgers::ComputeR(double x, double y) {
    double r = pow(x*x+y*y, 0.5);
    return r;
}
