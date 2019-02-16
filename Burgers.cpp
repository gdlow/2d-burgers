#include <cmath>

#include "Burgers.h"

Burgers::Burgers(Model &m) {
    model = &m;
}

void Burgers::SetInitialVelocity() {
    double r = ComputeR();
    double res;
    if (r <= 1) {
        res = pow(2.0*(1.0-r),4.0) * (4*r-1);
    }
    else res = 0;
    model->SetU0(res);
    model->SetV0(res);
}

double Burgers::IntegrateVelocity() {

}

void Burgers::WriteVelocityFile() {

}

double Burgers::CalculateEnergy() {

}

double Burgers::ComputeR() {
    double x = model->GetX0();
    double y = model->GetY0();
    double r = pow(x*x+y*y, 0.5);
    return r;
}