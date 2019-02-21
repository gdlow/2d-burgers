#include "Helpers.h"

/* Helper functions with 2D matrices */
double* unwrap(double** M, int n, int m) {
    double* res = new double[n*m];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            res[i*m+j] = M[i][j];
        }
    }
    return res;
}

double** wrap(double* A, int n, int m) {
    double ** res = new double*[n];
    for (int i = 0; i < n; i++) {
        res[i] = new double[m];
    }

    for (int i = 0; i < n*m; i++) {
        int row = i / m; // divisor result
        int col = i % m; // remainder after division
        res[row][col] = A[i];
    }
    return res;
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
double* GenSymm(double alpha, double beta, int N) {
    double* M = new double[N*N]; // Column-major format (symmetric anyway)
    for (int i = 0; i < N*N; i++) {
        M[i] = 0;
    }
    for (int i = 0; i < N; i++) {
        if (i>0) M[i*N+(i-1)] = beta;
        M[i*N+i] = alpha;
        if (i<N-1) M[i*N+i+1] = beta;
    }
    return M;
}

/**
 * Generates a triangular, square matrix
 * Stored in column-major format
 * alpha takes on the leading diagonal
 * beta takes on the banded row above (UPPER) or below
 * */
double* GenTrmm(double alpha, double beta, int N, bool UPPER) {
    double* M = new double[N*N];
    for (int i = 0; i < N*N; i++) {
        M[i] = 0;
    }
    for (int i = 0; i < N; i++) {
        if (UPPER && i>0) M[i*N+(i-1)] = beta;
        M[i*N+i] = alpha;
        if (!UPPER && i<N-1) M[i*N+i+1] = beta;
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
