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
    int Ny = model->GetNy();

    // Delete U and V
    for (int k = 1; k < Nt; k++) {
        // U[0] = V[0] = U0 (not dynamically alloc)
        for (int j = 0; j < Ny; j++) {
            delete[] U[k][j];
            delete[] V[k][j];
        }
        delete[] U[k];
        delete[] V[k];
    }
    delete[] U; delete[] V;
    U = nullptr; V = nullptr;

    // Delete U0
    for (int j = 0; j < Ny; j++) {
        delete[] U0[j];
    }
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
    U0 = new double*[Ny];
    for (int j = 0; j < Ny; j++) {
        U0[j] = new double[Nx];
        for (int i = 0; i < Nx; i++) {
            // Assumes x0 and y0 are identifying bottom LHS of matrix
            double y = y0 + j*dy;
            double x = x0 + i*dx;
            double r = ComputeR(x, y);
            U0[j][i] = (r <= 1.0)? pow(2.0*(1.0-r),4.0) * (4.0*r-1.0) : 0.0;
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

    // Compute U, V (coupled problem)
    U = nullptr; V = nullptr;
    U = new double**[Nt]; V = new double**[Nt];

    // Set initial velocity field
    U[0] = U0; V[0] = U0;

    // Compute velocity field states for k+1 = 1..Nt-1 => k = 0..Nt-2
    for (int k = 0; k < Nt-1; k++) {
        U[k+1] = new double*[Ny]; V[k+1] = new double*[Ny];
        // Compute velocity field states in (x,y)
        for (int j = 0; j < Ny; j++) {
            U[k+1][j] = new double[Nx]; V[k+1][j] = new double[Nx];
            for (int i = 0; i < Nx; i++) {
                // Set boundary conditions to == 0
                if (j == 0 || i == 0 || j == Ny-1 || i == Nx-1) {
                    U[k+1][j][i] = 0; V[k+1][j][i] = 0;
                }
                else {
                    // Compute differentials using FDS
                    double dudx = (U[k][j][i] - U[k][j][i-1]) / dx;
                    double dudy = (U[k][j][i] - U[k][j-1][i]) / dy;
                    double dudx_2 = (U[k][j][i+1] - 2.0*U[k][j][i] + U[k][j][i-1]) / pow(dx, 2.0);
                    double dudy_2 = (U[k][j+1][i] - 2.0*U[k][j][i] + U[k][j-1][i]) / pow(dy, 2.0);

                    double dvdx = (V[k][j][i] - V[k][j][i-1]) / dx;
                    double dvdy = (V[k][j][i] - V[k][j-1][i]) / dy;
                    double dvdx_2 = (V[k][j][i+1] - 2.0*V[k][j][i] + V[k][j][i-1]) / pow(dx, 2.0);
                    double dvdy_2 = (V[k][j+1][i] - 2.0*V[k][j][i] + V[k][j-1][i]) / pow(dy, 2.0);

                    // Compute U[k+1][j][i]
                    double temp_u = c*(dudx_2 + dudy_2) - (ax + b*U[k][j][i])*dudx - (ay + b*V[k][j][i])*dudy;
                    U[k+1][j][i] = temp_u * dt + U[k][j][i];

                    // Compute V[k+1][j][i]
                    double temp_v = c*(dvdx_2 + dvdy_2) - (ax + b*U[k][j][i])*dvdx - (ay + b*V[k][j][i])*dvdy;
                    V[k+1][j][i] = temp_v * dt + V[k][j][i];
                }
            }
        }
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

double Burgers::ComputeR(double x, double y) {
    double r = pow(x*x+y*y, 0.5);
    return r;
}
