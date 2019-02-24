#include <chrono>
#include <cstdlib>
#include "Model.h"
#include "Burgers2P.h"
#include "Burgers.h"

int main(int argc, char* argv[]) {
    Model m(argc, argv);

    // typedef std::chrono::high_resolution_clock hrc;
    // typedef std::chrono::milliseconds ms;

    if (argc == 9 && atoi(argv[8]) == 1) {
        Burgers2P b(m);
        // Call code to initialise the problem here;
        m.PrintParameters();

        // hrc::time_point start = hrc::now();

        // Call code to perform time integration here
        b.SetInitialVelocity();
        b.SetIntegratedVelocity();

        // hrc::time_point end = hrc::now();

        // Calculate final energy and write output
        b.SetEnergy();
        b.WriteVelocityFile();
    } else {
        Burgers b(m);
        // Call code to initialise the problem here;
        m.PrintParameters();

        // hrc::time_point start = hrc::now();

        // Call code to perform time integration here
        b.SetInitialVelocity();
        b.SetIntegratedVelocity();

        // hrc::time_point end = hrc::now();

        // Calculate final energy and write output
        b.SetEnergy();
        b.WriteVelocityFile();
    }

    return 0;
}
