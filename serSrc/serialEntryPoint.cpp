#include <chrono>
#include "Model.h"
#include "Burgers.h"
#include <iostream>

int main(int argc, char* argv[]) {
    Model m(argc, argv);

    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    typedef std::chrono::duration<double> fsec;

    Burgers b(m);
    // Call code to initialise the problem here;
    m.PrintParameters();

    hrc::time_point start = hrc::now();

    // Call code to perform time integration here
    b.SetInitialVelocity();
    b.SetIntegratedVelocity();

    hrc::time_point end = hrc::now();
    fsec elapsed_seconds = end-start;
    ms elapsed_ms = std::chrono::duration_cast<ms>(elapsed_seconds);
    std::cout << "Time elapsed: " << elapsed_seconds.count() << " s" << std::endl;
    std::cout << "Time elapsed: " << elapsed_ms.count() << " ms" << std::endl;

    // Calculate final energy and write output
    b.SetEnergy();
    b.WriteVelocityFile();

    return 0;
}
