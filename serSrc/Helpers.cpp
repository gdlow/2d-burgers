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

/**
 * @brief Generates a symmetric, banded matrix stored in banded-format
 * @param alpha constant along the leading diagonal
 * @param beta constant along the sub/super-diagonal
 * @param N number of rows/columns (should be a square matrix)
 * @param M pre-allocated with N*N memory
 * */
void GenSymmBanded(double alpha, double beta, int N, double* M) {
    /// Generate first row <=> upper diagonal
    for (int i = 1; i < N; i++) {
        M[i*N] = beta;
    }
    /// Generate second row <=> leading diagonal
    for (int i = 0; i < N; i++) {
        M[i*N+1] = alpha;
    }
    /// Generate third row <=> lower diagonal
    for (int i = 0; i < N-1; i++) {
        M[i*N+2] = beta;
    }
}

/**
 * @brief Generates a lower triangular, banded matrix stored in banded-format
 * @param alpha constant along the leading diagonal
 * @param beta constant along the sub-diagonal
 * @param N number of rows/columns (should be a square matrix)
 * @param M pre-allocated with N*N memory
 * */
void GenTrmmBanded(double alpha, double beta, int N, double* M) {
    /// Generate first row <=> leading diagonal
    for (int i = 0; i < N; i++) {
        M[i*N] = alpha;
    }
    /// Generate second row  <=> lower diagonal
    for (int i = 0; i < N-1; i++) {
        M[i*N+1] = beta;
    }
}