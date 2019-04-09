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
    explicit Burgers2P(Model &m);
    ~Burgers2P();

    void SetInitialVelocity();
    void SetIntegratedVelocity();
    void WriteVelocityFile();
    void SetEnergy();
    double GetE()     const { return E; }
private:

    void GetNextVelocities();
    void ComputeNextVelocityState();
    void FixNextVelocityBoundaries();
    void SetCaches();
    double CalculateEnergyState(double* Ui, double* Vi);
    void AssembleMatrix(double* Vel, double** M);
    void WriteOf(double* Vel, double** M, std::ofstream &of, char id);

    /// Burger parameters
    Model* model;
    double* U;
    double* V;
    double* NextU;
    double* NextV;
    double E;

    /// Caches for partitioning matrix
    double* upU;
    double* downU;
    double* leftU;
    double* rightU;
    double* myUpU;
    double* myDownU;
    double* myLeftU;
    double* myRightU;

    double* upV;
    double* downV;
    double* leftV;
    double* rightV;
    double* myUpV;
    double* myDownV;
    double* myLeftV;
    double* myRightV;

    /// MPI Requests and Statuses
    MPI_Request* reqs;
    MPI_Status* stats;
};
#endif //CLASS_BURGERS2P
