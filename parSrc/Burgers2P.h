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
    void SetMatrixCoefficients();
    double* NextVelocityState(double* Ui, double* Vi, bool U_OR_V);
    void SetCaches(double* Vel);
    void UpdateBoundsLinear(double* dVel_2, double* dVel);
    double CalculateEnergyState(double* Ui, double* Vi);
    void AssembleMatrix(double* Vel, double** M);
    void WriteOf(double* Vel, double** M, std::ofstream &of, char id);
    double ComputeR(double x, double y);

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

    /// MPI Requests and Statuses
    MPI_Request* reqs;
    MPI_Status* stats;
};
#endif //CLASS_BURGERS2P
