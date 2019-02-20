#ifndef BLAS_WRAPPER_H
#define BLAS_WRAPPER_H

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

void F77NAME(dtrmm) (
        const char& side, const char& uplo,
        const char& transa, const char& diag,
        const int& m, const int& n,
        const double& alpha, const double* A,
        const int& lda, const double* B,
        const int& ldb
);

}


#endif //BLAS_WRAPPER_H
