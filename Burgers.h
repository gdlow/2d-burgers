#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include "Model.h"

class Burgers {
public:
    Burgers(Model &m);
    ~Burgers();

    void SetInitialVelocity();
    void SetIntegratedVelocity();
    void WriteVelocityFile();
    void SetEnergy();
private:
    double ComputeR(double x, double y);

    Model* model;
    double** U;
    double** V;
    double* U0;
    double* E;
};
#endif //CLASS_BURGERS
