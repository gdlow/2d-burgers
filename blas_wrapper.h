#ifndef BLAS_WRAPPER_H
#define BLAS_WRAPPER_H


// define LAPACK and BLAS libraries in file
// link in your build file

#define F77NAME(x) x##_
extern "C" {
/* Level 1 functions */
double F77NAME(ddot)(
        const int& n,
        const double *x, const int& incx,
        const double *y, const int& incy
);

double F77NAME(dasum)(
        const int& n,
        const double *x, const int& incx
);

double F77NAME(idamax)(
        const int& n,
        const double *x, const int& incx
);

double F77NAME(dcopy)(
        const int& n,
        const double *x, const int& incx,
        const double *y, const int& incy
);
double F77NAME(daxpy)(
        const int& n, const double& alpha,
        const double *x,const int& incx,
        const double *y, const int& incy
);

/* Level 2 functions */
void F77NAME(dgemv)(
        const char& trans, const int& m,
        const int& n, const double& alpha,
        const double* A, const int& lda,
        const double* x, const int& incx,
        const double& beta, double* y,
        const int& incy
);

void F77NAME(dsymv)(
        const char& uplo, const int& n,
        const double& alpha, const double* A,
        const int& lda, const double* x,
        const int& incx, const double& beta,
        double* y, const int& incy
);

void F77NAME(dgbmv)(
        const char& trans, const int& m,
        const int& n, const int& kl,
        const int& ku, const double& alpha,
        const double* A, const int& lda,
        const double* x, const int& incx,
        const double& beta, double* y,
        const int& incy
);

/* Level 3 functions */

// C = A * B
// ldx = no. of cols in x (M in a NxM matrix)
void F77NAME(dgemm) (
        const char& trans, const char& transb,
        const int& m, const int& n,
        const int& k, const double& alpha,
        const double* A, const int& lda,
        const double* B, const int& ldb,
        const double& beta, double* C,
        const int& ldc
);

// A is the symmetric matrix
// control whether it appears on LHS or RHS below
void F77NAME(dsymm) (
        const char& side, const char& uplo,
        const int& m, const int& n,
        const double& alpha, const double* A,
        const int& lda, const double* B,
        const int& ldb, const double& beta,
        double* C, const int& ldc
);

/* LAPACK functions */

// solving a system of linear equations
void F77NAME(dgesv)(
        const int& n, const int& nrhs, const double * A,
        const int& lda, int * ipiv, double * B,
        const int& ldb, int& info
);
}

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

void matMul(double* A, int N_A, int M_A, double* B, int N_B, int M_B, double* C) {
    // result is [N_A x M_B]
    // N_B == M_A
    // supply result as A^T and B^T because you are storing it in columnar format!
    F77NAME(dgemm)('T','T', N_A, M_B, M_A, 1.0, A, M_A, B, M_B, 0.0, C, M_B);
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

#endif //BLAS_WRAPPER_H
