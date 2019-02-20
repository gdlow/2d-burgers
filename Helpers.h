#ifndef HELPERS_H
#define HELPERS_H
double* unwrap(double** M, int n, int m);
double** wrap(double* A, int n, int m);
double** transpose(double** A, int N, int M);
double* GenSymm(double alpha, double beta, int N);
double* GenTrmm(double alpha, double beta, int N, bool UPPER);
double* MatMul(double* Ui, double* Vi, int Ny, int Nx, bool offset_i, bool offset_j, double p);
#endif //HELPERS_H
