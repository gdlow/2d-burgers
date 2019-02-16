#include <chrono>

#include "Model.h"
#include "Burgers.h"

int main(int argc, char* argv[]) {

    // Call code to initialise the problem here;
    Model m(argc, argv);
    Burgers b(m);

    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

    // Call code to perform time integration here
    b.SetIntegratedVelocity();
    b.WriteVelocityFile();

    hrc::time_point end = hrc::now();

    // Calculate final energy and write output

    return 0;
}
