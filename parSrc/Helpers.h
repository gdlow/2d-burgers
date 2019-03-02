#ifndef HELPERS_H
#define HELPERS_H
double** wrap(double* A, int Nyr, int Nxr);
void wrap(double* A, int Nyr, int Nxr, double** res);
void printDebug(double* A, int Nyr, int Nxr, char c);
void GenSymm(double alpha, double beta, int Nyr, int Nxr, double* M);
void GenTrmm(double alpha, double beta, int Nyr, int Nxr, bool UPPER, double* M);
void GenSymmBanded(double alpha, double beta, int N, double* M);
void GenTrmmBanded(double alpha, double beta, int N, bool IS_UPPER, double* M);
#endif //HELPERS_H
