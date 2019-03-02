#ifndef CLASS_BURGERS2P
#define CLASS_BURGERS2P

#include "Model2P.h"
#include <fstream>

/**
 * @class Burgers2P
 * @brief Creates a Burgers instance that does computations on Burger's equation
 * */
class Burgers2P {
public:
    Burgers2P(Model &m);
    ~Burgers2P();

    void SetInitialVelocity();
    void SetIntegratedVelocity();
    void WriteVelocityFile();
    void SetEnergy();
    double GetE()     const { return E; }
private:
    void WriteOf(double* Vel, double** M, std::ofstream &of, char id);
    double ComputeR(double x, double y);
    double CalculateEnergyState(double* Ui, double* Vi);
    double* NextVelocityState(double* Ui, double* Vi, bool U_OR_V);
    void SetMatrixCoefficients();
    void SetCaches(double* Vel, MPI_Request* reqs);
    void AssembleMatrix(double* Vel, double** M);
    void UpdateBoundsLinear(double* dVel_dx_2, double* dVel_dy_2, double* dVel_dx, double* dVel_dy, MPI_Request* reqs, MPI_Status* stats);
    void CopyAndDelete(double* NextU, double* NextV);

    /// Burger parameters
    Model* model;
    double* U;
    double* V;
    double* dVel_dx_2_coeffs;
    double* dVel_dy_2_coeffs;
    double* dVel_dx_coeffs;
    double* dVel_dy_coeffs;
    double E;

    /// Term arrays
    double* dVel_2;
    double* dVel;

    /// Caches for partitioning matrix
    double* upVel;
    double* downVel;
    double* leftVel;
    double* rightVel;
    double* myUpVel;
    double* myDownVel;
    double* myLeftVel;
    double* myRightVel;
};
#endif //CLASS_BURGERS2P
