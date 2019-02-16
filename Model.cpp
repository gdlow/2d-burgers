#include <iostream>
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
        SetNumerics();
    } catch (ParseException &e) {
        cout << e.what() << endl;
    }
}

/**
 * Public member functions
 * */

/**
 * Parses parameters from command line into program
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
        cout << "Parameters have been saved." << endl;
    }
    else throw parseException;
}

void Model::SetU0(double u_value) {
    u0 = u_value;
}

void Model::SetV0(double v_value) {
    u0 = v_value;
}

/**
 * Private member functions
 * */

/**
 * Set appropriate values for various members
 * */
void Model::SetNumerics() {
    x0 = 0.0;
    y0 = 0.0;
    Nx = 100;
    Ny = 100;
    Nt = 100;
    // dx,dy and dt are dependent on L,T; and discretization Nx,Ny,Nt:

    dx = 0.01;
    dy = 0.01;
    dt = 0.01;
}