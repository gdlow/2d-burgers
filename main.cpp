#include <chrono>

#include "Model.h"
#include "Burgers.h"

int main(int argc, char* argv[]) {
    Model m(argc, argv);
    Burgers b(m);

    // Call code to initialise the problem here;
    m.PrintParameters();

    // typedef std::chrono::high_resolution_clock hrc;
    // typedef std::chrono::milliseconds ms;
    // hrc::time_point start = hrc::now();

    // Call code to perform time integration here
    b.SetInitialVelocity();
    b.SetIntegratedVelocity();

    // hrc::time_point end = hrc::now();

    // Calculate final energy and write output
    // b.SetEnergy(); // Suppress to only test velocity settings
    b.WriteVelocityFile();

    return 0;
}
