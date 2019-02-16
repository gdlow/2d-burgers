#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include "Model.h"

class Burgers {
public:
    Burgers(Model &m);
    ~Burgers();

    void SetInitialVelocity();
    double IntegrateVelocity();
    void WriteVelocityFile();
    double CalculateEnergy();
private:
    double ComputeR();
    Model* model;
};
#endif //CLASS_BURGERS
