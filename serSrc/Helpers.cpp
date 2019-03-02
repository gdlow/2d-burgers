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
 * @brief Print the 2D representation of a column-major 1D pointer for debugging
 * @param A 1D pointer in column-major format
 * @param Nyr Nyr
 * @param Nxr Nxr
 * @param c ID for printing
 * */
void printDebug(double* A, int Nyr, int Nxr, char c) {
    double ** res = wrap(A, Nyr, Nxr);
    std::cout << "PrintDebug " << c << ":" << std::endl;
    for (int j = 0; j < Nyr; j++) {
        for (int i = 0; i < Nxr; i++) {
            std::cout << res[j][i] << ' ';
        }
        std::cout << std::endl;
    }

    // destroy res
    for (int j = 0; j < Nyr; j++) {
        delete[] res[j];
    }
    delete[] res;
}

void GenSymmBanded(double alpha, double beta, int N, double* M) {
    /// M should be of size (N*N)
    /// IMPORTANT: This fills it in column-major format!

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