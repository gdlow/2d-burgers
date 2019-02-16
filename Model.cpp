#include <iostream>
#include "Model.h"

/**
 * CPP source file for Model class
 * Constructors, public and private member functions defined here
 * */

using namespace std;

/**
 * Constructor: sets constants from arg parameters
 * */

Model::Model(int argc, char** argv) {
    if (argc == 5) {
        ax = atof(argv[1]);
        ay = atof(argv[2]);
        b = atof(argv[3]);
        c = atof(argv[4]);
        cout << "Parameters have been saved." << endl;
    }
    else cout << "Wrong number of arguments. Expected 5" << endl;
}