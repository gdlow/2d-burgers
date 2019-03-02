#ifndef F77NAME
#define F77NAME(x) x##_
extern "C" {
/* Level 1 functions */
double F77NAME(ddot)(
        const int& n,
        const double *x, const int& incx,
        const double *y, const int& incy
);

double F77NAME(dcopy)(
        const int& n,
        const double *x, const int& incx,
        const double *y, const int& incy
);

/* Level 2 functions */
void F77NAME(dgbmv)(
        const char& trans, const int& m,
        const int& n, const int& kl,
        const int& ku, const double& alpha,
        const double* A, const int& lda,
        const double* x, const int& incx,
        const double& beta, double* y,
        const int& incy
);
}
#endif //F77NAME
