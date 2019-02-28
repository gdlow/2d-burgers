#include <iostream>
#include "Model.h"
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
}

/**
 * @brief Destructor: deallocates memory, finalizes MPI program and destroys Model instance
 * */
Model::~Model() {

}

/**
 * @brief Parses parameters from command line into program
 * @brief Throws an exception if invalid number of arguments are supplied
 * */
void Model::ParseParameters(int argc, char **argv) {
    if (argc == 8) {
        ax = atof(argv[1]);
        ay = atof(argv[2]);
        b = atof(argv[3]);
        c = atof(argv[4]);
        Lx = atof(argv[5]);
        Ly = atof(argv[6]);
        T = atof(argv[7]);
        cout << "Parameters saved successfully." << endl;
    }
    else throw illegalArgumentException;
}

/**
 * @brief Prints model parameters
 * */
void Model::PrintParameters() {
    cout << "ax: " << ax << endl;
    cout << "ay: " << ay << endl;
    cout << "b: " << b << endl;
    cout << "c: " << c << endl;
    cout << "Lx: " << Lx << endl;
    cout << "Ly: " << Ly << endl;
    cout << "T: " << T << endl;
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
    Nx = 101;
    Ny = 101;
    Nt = 101;
    // dx,dy and dt are dependent on L,T and Nx,Ny,Nt:
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    dt = T / (Nt-1);
    // x0 and y0 represent the top LHS of the matrix:
    x0 = -Lx/2.0;
    y0 = Ly/2.0;
}
