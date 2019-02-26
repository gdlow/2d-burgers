#include <chrono>
#include "Model2P.h"
#include "Burgers2P.h"

int main(int argc, char* argv[]) {
    Model m(argc, argv);

    // typedef std::chrono::high_resolution_clock hrc;
    // typedef std::chrono::milliseconds ms;
    Burgers2P b(m);
    // Call code to initialise the problem here;
    m.PrintParameters();

    // hrc::time_point start = hrc::now();

    // Call code to perform time integration here
    b.SetInitialVelocity();
    b.SetIntegratedVelocity();

    // hrc::time_point end = hrc::now();

    // Calculate final energy and write output
    // b.SetEnergy();
    // b.WriteVelocityFile();

    return 0;
}
