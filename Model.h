#ifndef CLASS_MODEL
#define CLASS_MODEL

// #include "mpi.h"

class Model {
public:
    Model(int argc, char* argv[]);
    ~Model();

    void PrintParameters();

    bool IsValid();

    // Getters
    bool   IsVerbose() const { return verbose; }
    bool   IsHelp()    const { return help; }
    double GetX0()     const { return x0; }
    double GetY0()     const { return y0; }
    double GetLx()     const { return Lx; }
    double GetLy()     const { return Ly; }
    double GetT()      const { return T; }
    int    GetNx()     const { return Nx; }
    int    GetNy()     const { return Ny; }
    int    GetNt()     const { return Nt; }
    double GetDx()     const { return dx; }
    double GetDy()     const { return dy; }
    double GetDt()     const { return dt; }
    double GetAx()     const { return ax; }
    double GetAy()     const { return ay; }
    double GetB()      const { return b; }
    double GetC()      const { return c; }

    // Add any other getters here...
    double GetU0()      const { return u0; }
    double GetV0()      const { return v0; }

    // Public Setters
    void SetU0(double u_value);
    void SetV0(double v_value);
private:
    void ParseParameters(int argc, char* argv[]);
    void ValidateParameters();

    // Private Setters
    void SetNumerics();
    bool verbose;
    bool help;

    // Numerics: Everything here has to be predefined
    // ParseParameters:
    double Lx;
    double Ly;
    double T;
    // SetNumerics:
    double x0;
    double y0;
    int    Nx;
    int    Ny;
    int    Nt;
    double dx;
    double dy;
    double dt;

    // Physics
    double ax;
    double ay;
    double b;
    double c;

    // Add any additional parameters here...

    // Added Numerics
    double u0;
    double v0;
};

#endif //CLASS_MODEL