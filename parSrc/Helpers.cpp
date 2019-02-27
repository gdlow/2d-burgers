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

/**
 * @brief Generates a symmetrical, square matrix for getting matrix coefficients
 * @brief Stored in column-major format
 * @param alpha prefactor along the leading diagonal
 * @param beta prefactor along the banded rows above and below
 * @param Nyr Nyr
 * @param Nxr Nxr
 * @param M pre-allocated matrix to be filled in symmetrically
 * */
void GenSymm(double alpha, double beta, int Nyr, int Nxr, double* M) {
    for (int i = 0; i < Nyr*Nxr; i++) {
        M[i] = 0;
    }
    for (int i = 0; i < Nxr; i++) {
        if (i>0) M[i*Nyr+(i-1)] = beta;
        M[i*Nyr+i] = alpha;
        if (i<Nxr-1) M[i*Nyr+(i+1)] = beta;
    }
}

/**
 * @brief Generates a symmetrical, square matrix for getting matrix coefficients
 * @brief Stored in column-major format
 * @param alpha prefactor along the leading diagonal
 * @param beta prefactor along the banded row either above OR below
 * @param Nyr Nyr
 * @param Nxr Nxr
 * @param UPPER specifies whether the matrix is upper triangular
 * @param M pre-allocated matrix to be filled in symmetrically
 * */
void GenTrmm(double alpha, double beta, int Nyr, int Nxr, bool UPPER, double* M) {
    for (int i = 0; i < Nyr*Nxr; i++) {
        M[i] = 0;
    }
    for (int i = 0; i < Nxr; i++) {
        if (UPPER && i>0) M[i*Nyr+(i-1)] = beta;
        M[i*Nyr+i] = alpha;
        if (!UPPER && i<Nxr-1) M[i*Nyr+i+1] = beta;
    }
}

/**
 * @brief Performs element-wise multiplication of 2 matrices
 * @param Ui pointer to column-major, Velocity matrix. Always the one offset
 * @param Vi pointer to column-major, Other matrix. Is the one NOT offset
 * @param Nyr Nyr
 * @param Nxr Nxr
 * @param offset_i is there an i offset?
 * @param offset_j is there a j offset?
 * @param p prefactor
 * */
double* MatMul(double* Ui, double* Vi, int Nyr, int Nxr, bool offset_i, bool offset_j, double p) {
    double* M = new double[Nyr*Nxr];
    if (offset_i) {
        for (int i = 0; i < Nyr*Nxr; i++) {
            if (i < Nyr) M[i] = 0;
            else M[i] = p * Ui[i-Nyr] * Vi[i];
        }
    }
    else if (offset_j) {
        for (int i = 0; i < Nyr*Nxr; i++) {
            if (i % Nyr == 0) M[i] = 0;
            else M[i] = p * Ui[i-1] * Vi[i];
        }
    }
    else {
        for (int i = 0; i <Nyr*Nxr; i++) {
            M[i] = p * Ui[i] * Vi[i];
        }
    }
    return M;
}

/**
 * @brief Sets the entire 1D array to 0
 * */
void SetZeroes(double* arr, int sz) {
    for (int i = 0; i < sz; i++) {
        arr[i] = 0;
    }
}