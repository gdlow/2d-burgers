#include "Helpers.h"
#include <iostream>


/* Helper functions with 2D matrices */

/**
 * Wraps a column-major 1D pointer into a row-major 2D array
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

void wrap(double* A, int Nyr, int Nxr, double** res) {
    for (int i = 0; i < Nyr*Nxr; i++) {
        int col = i / Nyr; // divisor result
        int row = i % Nyr; // remainder after division
        res[row][col] = A[i];
    }
}

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

double** transpose(double** A, int N, int M) {
    double** T = new double*[M];
    for (int i = 0; i < M; i++) {
        T[i] = new double[N];
        for (int j = 0; j < N; j++) {
            T[i][j] = A[j][i];
        }
    }
    return T;
}


/**
 * Generates a symmetrical, square matrix
 * Stored in column-major format
 * alpha takes on the leading diagonal
 * beta takes on the banded row above and below
 * */
double* GenSymm(double alpha, double beta, int Nyr, int Nxr) {
    double* M = new double[Nyr*Nxr]; // Column-major format (symmetric anyway)
    for (int i = 0; i < Nyr*Nxr; i++) {
        M[i] = 0;
    }
    for (int i = 0; i < Nxr; i++) {
        if (i>0) M[i*Nyr+(i-1)] = beta;
        M[i*Nyr+i] = alpha;
        if (i<Nxr-1) M[i*Nyr+(i+1)] = beta;
    }
    return M;
}

/**
 * Generates a triangular, square matrix
 * Stored in column-major format
 * alpha takes on the leading diagonal
 * beta takes on the banded row above (UPPER) or below
 * */
double* GenTrmm(double alpha, double beta, int Nyr, int Nxr, bool UPPER) {
    double* M = new double[Nyr*Nxr];
    for (int i = 0; i < Nyr*Nxr; i++) {
        M[i] = 0;
    }
    for (int i = 0; i < Nxr; i++) {
        if (UPPER && i>0) M[i*Nyr+(i-1)] = beta;
        M[i*Nyr+i] = alpha;
        if (!UPPER && i<Nxr-1) M[i*Nyr+i+1] = beta;
    }
    return M;
}

/**
 * Performs element-wise multiplication of 2 matrices
 * Ui is always the one offset (Must -> Vel)
 * p = prefactor
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

