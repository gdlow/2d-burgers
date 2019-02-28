#ifndef HELPERS_H
#define HELPERS_H
double** wrap(double* A, int Nyr, int Nxr);
void wrap(double* A, int Nyr, int Nxr, double** res);
void printDebug(double* A, int Nyr, int Nxr, char c);
void GenSymm(double alpha, double beta, int Nyr, int Nxr, double* M);
void GenTrmm(double alpha, double beta, int Nyr, int Nxr, bool UPPER, double* M);
double* MatMul(double* Ui, double* Vi, int Ny, int Nx, bool offset_i, bool offset_j, double p);
#endif //HELPERS_H
