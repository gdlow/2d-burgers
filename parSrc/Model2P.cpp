#include <iostream>
#include <mpi.h>
#include "Model2P.h"
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
    } catch (IllegalArgumentException &e) {
        cout << e.what() << endl;
    }
    ValidateParameters();

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &loc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    SetGridParameters();
    SetCartesianGrid();
}

Model::~Model() {
    delete[] loc_coord;
    delete[] loc_Nxr;
    delete[] loc_Nyr;
    delete[] displs_x;
    delete[] displs_y;
    loc_coord = nullptr;
    loc_Nxr = nullptr;
    loc_Nyr = nullptr;
    displs_x = nullptr;
    displs_y = nullptr;
    MPI_Finalize();
}

/**
 * Parses parameters from command line into program
 * */
void Model::ParseParameters(int argc, char **argv) {
    if (argc == 10) {
        ax = atof(argv[1]);
        ay = atof(argv[2]);
        b = atof(argv[3]);
        c = atof(argv[4]);
        Lx = atof(argv[5]);
        Ly = atof(argv[6]);
        T = atof(argv[7]);
        Px = atoi(argv[8]);
        Py = atoi(argv[9]);
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
        cout << "Px: " << Px << endl;
        cout << "Py: " << Py << endl;
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
    // Reduced parameters
    int Nxr = Nx - 2;
    int Nyr = Ny - 2;

    loc_Nxr = new int[Px];
    loc_Nyr = new int[Py];
    displs_x = new int[Px];
    displs_y = new int[Py];
    int rem, sum;

    // Define subcols and displs along x
    rem = Nxr % Px;
    sum = 0;
    for (int i = 0; i < Px; i++) {
        loc_Nxr[i] = Nxr / Px;
        if (rem > 0) {
            loc_Nxr[i]++;
            rem--;
        }
        displs_x[i] = sum;
        sum += loc_Nxr[i];
    }

    // Define subcols and displs along y
    rem = Nyr % Py;
    sum = 0;
    for (int i = 0; i < Py; i++) {
        loc_Nyr[i] = Nyr / Py;
        if (rem > 0) {
            loc_Nyr[i]++;
            rem--;
        }
        displs_y[i] = sum;
        sum += loc_Nyr[i];
    }

    // Print result
    if (loc_rank == 0) {
        for (int j = 0; j < Py; j++) {
            for (int i = 0; i < Px; i++) {
                cout << "(" << loc_Nyr[j] << "," << loc_Nxr[i] << ")" << ' ';
            }
            cout << endl;
        }
    }
}

void Model::SetCartesianGrid() {
    int dim[2] = {Py, Px};
    int period[2] = {0,0};
    int reorder = 1;
    loc_coord = new int[2];

    // Create cartesian grid of processes
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &vu);

    // Set local coordinates
    MPI_Cart_coords(vu, loc_rank, 2, loc_coord);

    // Set neighbours
    SetNeighbours();

    // Print loc_rank, coordinates, local Nxr and Nyr
    cout << "Rank: " << loc_rank << endl;
    cout << "Coordinates: (" << loc_coord[0] << "," << loc_coord[1] << ")" << endl;
    cout << "Nyr, Nxr: (" << GetLocNyr() << "," << GetLocNxr() << ")" << endl;
}

void Model::SetNeighbours() {
    MPI_Cart_shift(vu, 0, 1, &up, &down);
    MPI_Cart_shift(vu, 1, 1, &left, &right);
}

int Model::GetCoordX() {
    return loc_coord[1];
}

int Model::GetCoordY() {
    return loc_coord[0];
}

int Model::GetLocNxr() {
    return loc_Nxr[loc_coord[1]];
}

int Model::GetLocNyr() {
    return loc_Nyr[loc_coord[0]];
}

int Model::GetDisplX() {
    return displs_x[loc_coord[1]];
}

int Model::GetDisplY() {
    return displs_y[loc_coord[0]];
}

