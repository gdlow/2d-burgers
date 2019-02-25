#include <iostream>
#include <mpi.h>
#include "Model.h"
#include "ParseException.h"

/**
 * CPP source file for Model class
 * Constructors, public and private member functions defined here
 * */

using namespace std;

/**
 * Constructor: sets constants from arg parameters
 * */
Model::Model(int argc, char** argv) {
    try {
        ParseParameters(argc, argv);
        if (Is_Parallel) {
            MPI_Init(&argc, &argv);
            MPI_Comm_rank(MPI_COMM_WORLD, &loc_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &p);
            SetGridParameters();
        }
    } catch (IllegalArgumentException &e) {
        cout << e.what() << endl;
    }
    // validate then SetNumerics
    ValidateParameters();
}

Model::~Model() {
    if (Is_Parallel) {
        delete[] loc_Nxr;
        delete[] loc_Nyr;
        delete[] displs_x;
        delete[] displs_y;
        loc_Nxr = nullptr;
        loc_Nyr = nullptr;
        displs_x = nullptr;
        displs_y = nullptr;
        MPI_Finalize();
    }
}

/**
 * Parses parameters from command line into program
 * */
void Model::ParseParameters(int argc, char **argv) {
    if (argc == 9) {
        ax = atof(argv[1]);
        ay = atof(argv[2]);
        b = atof(argv[3]);
        c = atof(argv[4]);
        Lx = atof(argv[5]);
        Ly = atof(argv[6]);
        T = atof(argv[7]);
        Is_Parallel = atoi(argv[8]);
        // Last parameter for choosing Serial (0) or Parallel (1)
        cout << "Parameters saved successfully." << endl;
    }
    else throw illegalArgumentException;
}

/**
 * Prints model parameters
 * */
void Model::PrintParameters() {
    if (loc_rank == 0) {
        cout << "ax: " << ax << endl;
        cout << "ay: " << ay << endl;
        cout << "b: " << b << endl;
        cout << "c: " << c << endl;
        cout << "Lx: " << Lx << endl;
        cout << "Ly: " << Ly << endl;
        cout << "T: " << T << endl;
    }
}

/**
 * Checks if parameters supplied are valid
 * */
bool Model::IsValid() {
    return ax >= 0 && ay >= 0 && b >= 0 && c >= 0 && Lx >= 0 && Ly >= 0 && T >= 0;
}

void Model::ValidateParameters() {
    if (!IsValid()) cout << "WARN: Parameter values have to be (>=0)" << endl;
    else SetNumerics();
}

/**
 * Set appropriate values for various members
 * */
void Model::SetNumerics() {
    Nx = 11;
    Ny = 11;
    Nt = 11;
    // dx,dy and dt are dependent on L,T and Nx,Ny,Nt:
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    dt = T / (Nt-1);
    // x0 and y0 represent the top LHS of the matrix:
    x0 = -Lx/2.0;
    y0 = Ly/2.0;
}

void Model::SetGridParameters() {
    int Ny_f = 201;
    int Nx_f = 201;
    // Px, Py to be redefined as instance variables
    // parsed into command line
    int Px = 10;
    int Py = 10;

    loc_Nxr = new int[Px];
    loc_Nyr = new int[Py];
    displs_x = new int[Px];
    displs_y = new int[Py];

    // Reduced parameters
    int Nxr = Nx_f - 2;
    int Nyr = Ny_f - 2;

    // Define subcols and displs along x
    int rem_x = Nxr % Px;
    int sum_x = 0;
    for (int i = 0; i < Px; i++) {
        loc_Nxr[i] = Nxr / Px;
        if (rem_x > 0) {
            loc_Nxr[i]++;
            rem_x--;
        }
        displs_x[i] = sum_x;
        sum_x += loc_Nxr[i];
    }

    // Define subcols and displs along y
    int rem_y = Nyr % Py;
    int sum_y = 0;
    for (int i = 0; i < Py; i++) {
        loc_Nyr[i] = Nyr / Py;
        if (rem_y > 0) {
            loc_Nyr[i]++;
            rem_y--;
        }
        displs_y[i] = sum_y;
        sum_y += loc_Nyr[i];
    }

    // Print result
    if (loc_rank == 0) {
        for (int j = 0; j < Py; j++) {
            for (int i = 0; i < Px; i++) {
                cout << "(" << loc_Nxr[i] << "," << loc_Nyr[j] << ")" << ' ';
            }
            cout << endl;
        }
    }
}