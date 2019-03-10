#include <iostream>
#include <mpi.h>
#include <cmath>
#include "Model2P.h"
#include "ParseException.h"

using namespace std;

/**
 * @brief Constructor: sets constants from arg parameters
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

/**
 * @brief Destructor: deallocates memory, finalizes MPI program and destroys Model instance
 * */
Model::~Model() {
    delete[] loc_coord;
    delete[] loc_Nxr;
    delete[] loc_Nyr;
    delete[] displs_x;
    delete[] displs_y;
    delete[] displs;
    delete[] recvcount;
    delete[] rankNxrMap;
    delete[] rankNyrMap;
    delete[] rankDisplsXMap;
    delete[] rankDisplsYMap;
    MPI_Finalize();
}

/**
 * @brief Parses parameters from command line into program
 * Throws an exception if invalid number of arguments are supplied
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
    }
    else throw illegalArgumentException;
}

/**
 * @brief Prints model parameters
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
 * @brief Checks if parameters supplied are valid
 * */
bool Model::IsValid() {
    return ax >= 0 && ay >= 0 && b >= 0 && c >= 0 && Lx >= 0 && Ly >= 0 && T >= 0;
}

/**
 * @brief Validates the parameters. If parameters supplied are valid, set them as instance vars
 * */
void Model::ValidateParameters() {
    if (!IsValid()) cout << "WARN: Parameter values have to be (>=0)" << endl;
    else SetNumerics();
}

/**
 * @brief Set appropriate values for various members
 * */
void Model::SetNumerics() {
    Nx = 2001;
    Ny = 2001;
    Nt = 4001;
    /// dx,dy and dt are dependent on L,T and Nx,Ny,Nt:
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    dt = T / (Nt-1);
    /// x0 and y0 represent the top LHS of the matrix:
    x0 = -Lx/2.0;
    y0 = Ly/2.0;
    /// b/dx and b/dy saves computation time in the future
    bdx = b/dx;
    bdy = b/dy;
    /// constants used in SetIntegratedVelocity()
    double alpha_dx_2 = (-2.0*c)/pow(dx,2.0);
    double alpha_dy_2 = (-2.0*c)/pow(dy,2.0);
    double alpha_dx_1 = -ax/dx;
    double alpha_dy_1 = -ay/dy;
    double beta_dx_1 = ax/dx;
    double beta_dy_1 = ay/dy;
    beta_dx_2 = c/pow(dx,2.0);
    beta_dy_2 = c/pow(dy,2.0);
    alpha_sum = alpha_dx_1 + alpha_dx_2 + alpha_dy_1 + alpha_dy_2;
    beta_dx_sum = beta_dx_1 + beta_dx_2;
    beta_dy_sum = beta_dy_1 + beta_dy_2;
    /// multiply by dt for pre-computational purposes
    bdx *= dt;
    bdy *= dt;
    alpha_sum *= dt;
    beta_dx_sum *= dt;
    beta_dy_sum *= dt;
    beta_dx_2 *= dt;
    beta_dy_2 *= dt;
}

/**
 * @brief Sets the local and global displacements and sizes of each sub-matrix
 * */
void Model::SetGridParameters() {
    /// Reduced parameters
    int Nxr = Nx - 2;
    int Nyr = Ny - 2;

    loc_Nxr = new int[Px];
    loc_Nyr = new int[Py];
    displs_x = new int[Px];
    displs_y = new int[Py];
    int rem, sum;

    /// Define subcols and displs along x
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

    /// Define subcols and displs along y
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

    /// Define parameters for assembling matrix
    recvcount = new int[Px*Py];
    displs = new int[Px*Py];
    rankNxrMap = new int[Px*Py];
    rankNyrMap = new int[Px*Py];
    rankDisplsXMap = new int[Px*Py];
    rankDisplsYMap = new int[Px*Py];

    /* Layout of cartesian grid in row-major format */
    sum = 0;
    for (int j = 0; j < Py; j++) {
        for (int i = 0; i < Px; i++) {
            rankNxrMap[j*Px+i] = loc_Nxr[i];
            rankNyrMap[j*Px+i] = loc_Nyr[j];
            rankDisplsXMap[j*Px+i] = displs_x[i];
            rankDisplsYMap[j*Px+i] = displs_y[j];
            displs[j*Px+i] = sum;
            recvcount[j*Px+i] = loc_Nyr[j] * loc_Nxr[i];
            sum += recvcount[j*Px+i];
        }
    }

    /// Print result
    if (loc_rank == 0) {
        for (int j = 0; j < Py; j++) {
            for (int i = 0; i < Px; i++) {
                cout << "(" << loc_Nyr[j] << "," << loc_Nxr[i] << ")" << ' ';
            }
            cout << endl;
        }
    }
}

/**
 * @brief Sets up a cartesian grid of Px * Py processors and identifies local neighbours
 * */
void Model::SetCartesianGrid() {
    int dim[2] = {Py, Px};
    int period[2] = {0,0};
    int reorder = 1;
    loc_coord = new int[2];

    /// Create cartesian grid of processes
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &vu);

    /// Recast loc_rank and p wrt vu
    MPI_Comm_rank(vu, &loc_rank);
    MPI_Comm_size(vu, &p);

    /// Set local coordinates
    MPI_Cart_coords(vu, loc_rank, 2, loc_coord);

    /// Set neighbours
    SetNeighbours();

    /// Print loc_rank, coordinates, local Nxr and Nyr
    cout << "Rank: " << loc_rank << endl;
    cout << "Coordinates: (" << loc_coord[0] << "," << loc_coord[1] << ")" << endl;
    cout << "Nyr, Nxr: (" << GetLocNyr() << "," << GetLocNxr() << ")" << endl;
}

/**
 * @brief Sets processor ids of neighbours to the local sub-matrix
 * */
void Model::SetNeighbours() {
    MPI_Cart_shift(vu, 0, 1, &up, &down);
    MPI_Cart_shift(vu, 1, 1, &left, &right);
}

/**
 * @brief Get local Nxr
 * */
int Model::GetLocNxr() const {
    return loc_Nxr[loc_coord[1]];
}

/**
 * @brief Get local Nyr
 * */
int Model::GetLocNyr() const {
    return loc_Nyr[loc_coord[0]];
}

/**
 * @brief Get local x displacement from global (0,0)
 * */
int Model::GetDisplX() const {
    return displs_x[loc_coord[1]];
}

/**
 * @brief Get local y displacement from global (0,0)
 * */
int Model::GetDisplY() const {
    return displs_y[loc_coord[0]];
}

/**
 * @brief Gets submatrix size
 * */
int Model::GetLocNyrNxr() const {
    return GetLocNxr() * GetLocNyr();
}

