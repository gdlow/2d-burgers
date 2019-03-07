#include "Helpers.h"
#include <iostream>

/* Helper functions with 2D matrices */

/**
 * @brief Wraps a column-major 1D pointer into a generated row-major 2D pointer
 * @param A 1D pointer in column-major format
 * @param Nxr Nxr
 * @param Nyr Nyr
 * */
double** wrap(double* A, int Nyr, int Nxr) {
    double ** res = new double*[Nyr];
    for (int i = 0; i < Nyr; i++) {
        res[i] = new double[Nxr];
    }

    for (int i = 0; i < Nyr*Nxr; i++) {
        int col = i / Nyr; // divisor result
        int row = i % Nyr; // remainder after division
        res[row][col] = A[i];
    }
    return res;
}

/**
 * @brief Wraps a column-major 1D pointer into a pre-allocated row-major 2D pointer
 * @param A 1D pointer in column-major format
 * @param Nxr Nxr
 * @param Nyr Nyr
 * @param res 2D pointer pre-allocated with memory
 * */
void wrap(double* A, int Nyr, int Nxr, double** res) {
    for (int i = 0; i < Nyr*Nxr; i++) {
        int col = i / Nyr; // divisor result
        int row = i % Nyr; // remainder after division
        res[row][col] = A[i];
    }
}